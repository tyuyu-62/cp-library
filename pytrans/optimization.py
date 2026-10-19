from __future__ import annotations

import ast
from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class BindingOptimizations:
    callable_aliases: frozenset[str]
    private_lists: frozenset[str]


class OptimizationPlanner:
    def __init__(self, analyzer: Any) -> None:
        self.analyzer = analyzer
        self.constant_globals = self._constant_globals()
        self.immutable_callables = self._immutable_callables()

    def _constant_globals(self) -> dict[str, ast.stmt]:
        tree = self.analyzer.tree
        if any(
            isinstance(child, ast.Call)
            and isinstance(child.func, ast.Name)
            and child.func.id in {"globals", "locals", "vars", "eval", "exec"}
            for child in ast.walk(tree)
        ):
            return {}
        stores: dict[str, list[ast.Name]] = {}
        for child in ast.walk(tree):
            if isinstance(child, ast.Name) and not isinstance(child.ctx, ast.Load):
                stores.setdefault(child.id, []).append(child)
        result = {}
        prior_loads: set[str] = set()
        prefix_safe = True
        for statement in tree.body:
            target = None
            if isinstance(statement, ast.Assign) and len(statement.targets) == 1:
                target = statement.targets[0]
            elif isinstance(statement, ast.AnnAssign):
                target = statement.target
            value = getattr(statement, "value", None)
            literal = object()
            try:
                literal = ast.literal_eval(value)
            except (ValueError, TypeError, SyntaxError):
                pass
            if (
                prefix_safe
                and isinstance(target, ast.Name)
                and len(stores.get(target.id, [])) == 1
                and target.id not in prior_loads
                and isinstance(literal, int)
                and -(1 << 63) < literal < (1 << 63)
                and getattr(self.analyzer.global_types.get(target.id), "kind", "unknown") in {"i64", "bool"}
            ):
                result[target.id] = statement
            prior_loads.update(child.id for child in ast.walk(statement) if isinstance(child, ast.Name) and isinstance(child.ctx, ast.Load))
            if isinstance(statement, (ast.Import, ast.ImportFrom)):
                continue
            if (
                isinstance(statement, ast.FunctionDef)
                and not statement.decorator_list
                and all(
                    isinstance(default, ast.Constant)
                    for default in statement.args.defaults + [item for item in statement.args.kw_defaults if item is not None]
                )
                and not statement.args.vararg
                and not statement.args.kwarg
            ):
                continue
            if isinstance(target, ast.Name) and isinstance(literal, (int, float, str, bool, type(None))) and value is not None:
                continue
            prefix_safe = False
        return result

    def _immutable_callables(self) -> dict[str, frozenset[str]]:
        tree = self.analyzer.tree
        if any(
            isinstance(child, ast.Call)
            and (
                (isinstance(child.func, ast.Attribute) and child.func.attr == "__init__")
                or (isinstance(child.func, ast.Name) and child.func.id in {"setattr", "delattr", "globals", "eval", "exec"})
            )
            for child in ast.walk(tree)
        ):
            return {}
        parents = {child: parent for parent in ast.walk(tree) for child in ast.iter_child_nodes(parent)}
        result = {}
        for name, attributes in self.analyzer.class_attributes.items():
            candidates = {attribute for attribute, kind in attributes.items() if kind.kind == "function"}
            for child in ast.walk(tree):
                if not isinstance(child, ast.Attribute) or isinstance(child.ctx, ast.Load) or child.attr not in candidates:
                    continue
                scope = parents.get(child)
                while scope is not None and not isinstance(scope, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                    scope = parents.get(scope)
                owner = parents.get(scope)
                if not (
                    isinstance(scope, ast.FunctionDef)
                    and scope.name == "__init__"
                    and isinstance(owner, ast.ClassDef)
                    and owner.name == name
                    and isinstance(child.value, ast.Name)
                    and child.value.id == "self"
                ):
                    candidates.discard(child.attr)
            result[name] = frozenset(candidates)
        return result

    def local_bindings(self, node: ast.FunctionDef, env: dict[str, Any]) -> BindingOptimizations:
        parents = {child: parent for parent in ast.walk(node) for child in ast.iter_child_nodes(parent)}
        names: dict[str, list[ast.Name]] = {}
        for child in ast.walk(node):
            if child is not node and isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda, ast.ClassDef)):
                return BindingOptimizations(frozenset(), frozenset())
            if isinstance(child, (ast.Global, ast.Nonlocal)):
                return BindingOptimizations(frozenset(), frozenset())
            if isinstance(child, ast.Name):
                if self.analyzer.hoisted_captures.get(child.id):
                    return BindingOptimizations(frozenset(), frozenset())
                names.setdefault(child.id, []).append(child)
        aliases: set[str] = set()
        private: set[str] = set()
        parameters = {
            argument.arg
            for argument in node.args.posonlyargs + node.args.args + node.args.kwonlyargs
        }
        for statement in node.body:
            if not isinstance(statement, ast.Assign) or len(statement.targets) != 1 or not isinstance(statement.targets[0], ast.Name):
                continue
            name = statement.targets[0].id
            if name in parameters:
                continue
            uses = names.get(name, [])
            if sum(not isinstance(use.ctx, ast.Load) for use in uses) != 1:
                continue
            loads = [use for use in uses if isinstance(use.ctx, ast.Load)]
            if any(use.lineno <= statement.lineno for use in loads):
                continue
            value = statement.value
            value_type = env.get(name)
            if (
                getattr(value_type, "kind", "unknown") == "function"
                and isinstance(value, ast.Attribute)
                and isinstance(value.value, ast.Name)
                and value.value.id == "self"
                and getattr(env.get("self"), "kind", "unknown") == "class"
                and (
                    value.attr in self.analyzer.class_methods(env["self"].name)
                    or (node.name != "__init__" and value.attr in self.immutable_callables.get(env["self"].name, ()))
                )
                and all(isinstance(parents[use], ast.Call) and parents[use].func is use for use in loads)
            ):
                aliases.add(name)
            item_type = None
            if value_type is not None and getattr(value_type, "args", ()):
                item_type = value_type.args[0]
            if (
                not isinstance(value, ast.List)
                or getattr(value_type, "kind", "unknown") != "list"
                or getattr(item_type, "kind", "unknown") not in {"i64", "f64", "bool", "char", "string", "bytes"}
            ):
                continue
            safe = True
            for use in loads:
                parent = parents[use]
                if isinstance(parent, (ast.While, ast.If)) and parent.test is use:
                    continue
                if isinstance(parent, ast.Attribute) and parent.value is use and parent.attr in {"append", "pop"}:
                    call = parents.get(parent)
                    if isinstance(call, ast.Call) and call.func is parent and not call.keywords and len(call.args) == (1 if parent.attr == "append" else 0):
                        continue
                if isinstance(parent, ast.Call) and isinstance(parent.func, ast.Name) and parent.func.id == "reversed" and len(parent.args) == 1 and not parent.keywords:
                    loop = parents.get(parent)
                    if isinstance(loop, ast.For) and loop.iter is parent and not any(
                        isinstance(child, ast.Name) and child.id == name
                        for statement in loop.body for child in ast.walk(statement)
                    ):
                        continue
                safe = False
                break
            if safe:
                private.add(name)
        return BindingOptimizations(frozenset(aliases), frozenset(private))

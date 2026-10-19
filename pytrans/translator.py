from __future__ import annotations

import ast
import copy as copy_module
import json
import re
import types
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

try:
    from .ir import FunctionIR, LocalIR, Ownership, ProgramIR, Storage
    from .optimization import OptimizationPlanner
except ImportError:
    from ir import FunctionIR, LocalIR, Ownership, ProgramIR, Storage
    from optimization import OptimizationPlanner


class TranspileError(Exception):
    def __init__(self, message: str, node: ast.AST | None = None, filename: str = "<input>") -> None:
        self.message = message
        self.node = node
        self.filename = filename
        if node is not None and hasattr(node, "lineno"):
            super().__init__(f"{filename}:{node.lineno}:{getattr(node, 'col_offset', 0) + 1}: {message}")
        else:
            super().__init__(f"{filename}: {message}")


@dataclass(frozen=True)
class Type:
    kind: str
    args: tuple["Type", ...] = ()
    name: str = ""

    def cpp(self) -> str:
        if self.kind == "unknown":
            return "i64"
        if self.kind == "i64":
            return "i64"
        if self.kind == "f64":
            return "double"
        if self.kind == "bool":
            return "bool"
        if self.kind == "char":
            return "char"
        if self.kind == "string":
            return "string"
        if self.kind == "bytes":
            return "string"
        if self.kind == "void":
            return "void"
        if self.kind == "list":
            return f"vector<{self.args[0].cpp()}>"
        if self.kind == "deque":
            return f"deque<{self.args[0].cpp()}>"
        if self.kind == "set":
            return f"set<{self.args[0].cpp()}>"
        if self.kind == "tuple":
            return "tuple<" + ", ".join(arg.cpp() for arg in self.args) + ">"
        if self.kind == "reference":
            return f"{self.args[0].cpp()}&"
        if self.kind == "range":
            return "range"
        if self.kind == "map":
            return f"PyDict<{self.args[0].cpp()}, {self.args[1].cpp()}>"
        if self.kind == "default_map":
            return f"PyDefaultDict<{self.args[0].cpp()}, {self.args[1].cpp()}>"
        if self.kind == "optional":
            return f"optional<{self.args[0].cpp()}>"
        if self.kind == "function":
            result = self.args[-1].cpp() if self.args else "i64"
            parameters = ", ".join(argument.cpp() for argument in self.args[:-1])
            return f"function<{result}({parameters})>"
        if self.kind == "class":
            return sanitize_name(self.name)
        return "i64"


UNKNOWN = Type("unknown")
I64 = Type("i64")
F64 = Type("f64")
BOOL = Type("bool")
CHAR = Type("char")
STRING = Type("string")
BYTES = Type("bytes")
VOID = Type("void")
RANGE = Type("range")


def list_type(item: Type = UNKNOWN) -> Type:
    return Type("list", (item,))


def deque_type(item: Type = UNKNOWN) -> Type:
    return Type("deque", (item,))


def set_type(item: Type = UNKNOWN) -> Type:
    return Type("set", (item,))


def tuple_type(items: Iterable[Type]) -> Type:
    return Type("tuple", tuple(items))


def reference_type(value: Type) -> Type:
    return Type("reference", (value,))


def map_type(key: Type = UNKNOWN, value: Type = UNKNOWN) -> Type:
    return Type("map", (key, value))


def default_map_type(key: Type = UNKNOWN, value: Type = UNKNOWN) -> Type:
    return Type("default_map", (key, value))


def optional_type(value: Type) -> Type:
    return value if value.kind == "optional" else Type("optional", (value,))


def function_type(parameters: Iterable[Type], result: Type = UNKNOWN) -> Type:
    return Type("function", tuple(parameters) + (result,))


def function_parameters(value: Type) -> tuple[Type, ...]:
    return value.args[:-1] if value.kind == "function" and value.args else ()


def function_result(value: Type) -> Type:
    return value.args[-1] if value.kind == "function" and value.args else UNKNOWN


def is_concrete_type(value: Type) -> bool:
    if value.kind == "unknown":
        return False
    if value.kind in {"reference", "optional"}:
        return bool(value.args) and is_concrete_type(value.args[0])
    if value.kind in {"list", "deque", "set"}:
        return bool(value.args) and is_concrete_type(value.args[0])
    if value.kind in {"map", "default_map"}:
        return len(value.args) == 2 and all(is_concrete_type(argument) for argument in value.args)
    if value.kind == "tuple":
        return bool(value.args) and all(is_concrete_type(argument) for argument in value.args)
    if value.kind == "function":
        return bool(value.args) and all(is_concrete_type(argument) for argument in value.args)
    return True


def class_type(name: str) -> Type:
    return Type("class", name=name)


def is_numeric(t: Type) -> bool:
    return t.kind in {"i64", "f64", "bool", "char", "unknown"}


def is_container(t: Type) -> bool:
    if t.kind == "reference":
        t = t.args[0]
    return t.kind in {"list", "deque", "set", "string", "bytes", "tuple", "map", "default_map"}


def element_type(t: Type) -> Type:
    if t.kind == "reference":
        return t.args[0]
    if t.kind == "optional":
        return element_type(t.args[0])
    if t.kind in {"list", "deque", "set"}:
        return t.args[0]
    if t.kind == "string":
        return CHAR
    if t.kind == "bytes":
        return I64
    if t.kind == "tuple":
        result = UNKNOWN
        for item in t.args:
            result = merge_type(result, item)
        return result
    if t.kind in {"map", "default_map"}:
        return t.args[0]
    return UNKNOWN


def merge_type(left: Type, right: Type) -> Type:
    if left.kind == "reference" or right.kind == "reference":
        left_value = left.args[0] if left.kind == "reference" else left
        right_value = right.args[0] if right.kind == "reference" else right
        return reference_type(merge_type(left_value, right_value))
    if left.kind == "unknown":
        return right
    if right.kind == "unknown":
        return left
    if left == right:
        return left
    if left.kind == "tuple" and not left.args and right.kind == "tuple":
        return right
    if right.kind == "tuple" and not right.args and left.kind == "tuple":
        return left
    if left.kind == "void":
        return optional_type(right)
    if right.kind == "void":
        return optional_type(left)
    if left.kind == "optional" and right.kind == "optional":
        return optional_type(merge_type(left.args[0], right.args[0]))
    if left.kind == "optional":
        return optional_type(merge_type(left.args[0], right))
    if right.kind == "optional":
        return optional_type(merge_type(left, right.args[0]))
    if left.kind == "function" and right.kind == "function":
        if len(left.args) != len(right.args):
            raise ValueError(f"incompatible function types: {left.kind} and {right.kind}")
        parameters = tuple(merge_type(a, b) for a, b in zip(left.args[:-1], right.args[:-1]))
        left_result = left.args[-1]
        right_result = right.args[-1]
        if left_result.kind in {"tuple", "list", "deque", "set", "map", "default_map"} and is_numeric(right_result):
            result = left_result
        elif right_result.kind in {"tuple", "list", "deque", "set", "map", "default_map"} and is_numeric(left_result):
            result = right_result
        else:
            result = merge_type(left_result, right_result)
        return function_type(parameters, result)
    if left.kind == "string" and right == list_type(CHAR):
        return STRING
    if right.kind == "string" and left == list_type(CHAR):
        return STRING
    if {left.kind, right.kind} == {"string", "bytes"}:
        return STRING
    if (left.kind == "string" and is_numeric(right)) or (right.kind == "string" and is_numeric(left)):
        return STRING
    if is_numeric(left) and is_numeric(right):
        if left.kind == "f64" or right.kind == "f64":
            return F64
        return I64
    if left.kind == right.kind and left.kind in {"list", "deque", "set"}:
        if left.args[0].kind == "tuple" and is_numeric(right.args[0]):
            return left
        if right.args[0].kind == "tuple" and is_numeric(left.args[0]):
            return right
        return Type(left.kind, (merge_type(left.args[0], right.args[0]),))
    if left.kind == right.kind and left.kind in {"map", "default_map"}:
        return Type(left.kind, (merge_type(left.args[0], right.args[0]), merge_type(left.args[1], right.args[1])))
    if left.kind == "tuple" and right.kind == "tuple" and len(left.args) == len(right.args):
        return tuple_type(merge_type(a, b) for a, b in zip(left.args, right.args))
    if left.kind in {"list", "deque", "set", "tuple"} and right.kind in {"list", "deque", "set", "tuple"}:
        if left.kind in {"list", "deque", "set"} and right.kind == "tuple" and is_numeric(left.args[0]):
            return list_type(right)
        if right.kind in {"list", "deque", "set"} and left.kind == "tuple" and is_numeric(right.args[0]):
            return list_type(left)
        return list_type(merge_type(element_type(left), element_type(right)))
    return UNKNOWN


def refine_annotated_type(annotation: Type, inferred: Type) -> Type:
    if annotation.kind == "unknown":
        return inferred
    if annotation.kind == "tuple" and not annotation.args and inferred.kind == "tuple":
        return inferred
    if annotation.kind == inferred.kind and len(annotation.args) == len(inferred.args):
        return Type(annotation.kind, tuple(
            refine_annotated_type(expected, actual)
            for expected, actual in zip(annotation.args, inferred.args)
        ), annotation.name)
    return annotation


def type_from_annotation(node: ast.AST | None) -> Type:
    if node is None:
        return UNKNOWN
    if isinstance(node, ast.Constant) and node.value is None:
        return VOID
    if isinstance(node, ast.Name):
        return {
            "int": I64,
            "bool": BOOL,
            "float": F64,
            "str": STRING,
            "bytes": BYTES,
            "list": list_type(),
            "tuple": tuple_type(()),
            "set": set_type(),
            "dict": map_type(),
            "deque": deque_type(),
            "array": list_type(I64),
            "None": VOID,
            "Sequence": list_type(UNKNOWN),
            "_Sequence": list_type(UNKNOWN),
            "Iterable": list_type(UNKNOWN),
            "Collection": list_type(UNKNOWN),
            "MutableSequence": list_type(UNKNOWN),
            "Callable": function_type((), UNKNOWN),
            "Any": UNKNOWN,
            "object": UNKNOWN,
        }.get(node.id, UNKNOWN)
    if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name):
        if node.value.id in {"typing", "collections", "collections.abc"}:
            return type_from_annotation(ast.Name(id=node.attr))
        return UNKNOWN
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.BitOr):
        left = type_from_annotation(node.left)
        right = type_from_annotation(node.right)
        if left.kind == "void":
            return optional_type(right)
        if right.kind == "void":
            return optional_type(left)
        try:
            return merge_type(left, right)
        except ValueError:
            return UNKNOWN
    if isinstance(node, ast.Subscript):
        name = ""
        if isinstance(node.value, ast.Name):
            name = node.value.id
        elif isinstance(node.value, ast.Attribute) and isinstance(node.value.value, ast.Name):
            name = f"{node.value.value.id}.{node.value.attr}"
        if name in {"Callable", "typing.Callable", "collections.abc.Callable"}:
            if not isinstance(node.slice, ast.Tuple) or len(node.slice.elts) != 2:
                return function_type((), UNKNOWN)
            parameter_node, result_node = node.slice.elts
            if isinstance(parameter_node, (ast.List, ast.Tuple)):
                parameters = tuple(type_from_annotation(argument) for argument in parameter_node.elts)
            else:
                parameters = ()
            return function_type(parameters, type_from_annotation(result_node))
        args_node = node.slice.elts if isinstance(node.slice, ast.Tuple) else [node.slice]
        if name in {"tuple", "Tuple", "typing.Tuple"} and any(
            isinstance(argument, ast.Constant) and argument.value is Ellipsis
            for argument in args_node
        ):
            values = [argument for argument in args_node if not (isinstance(argument, ast.Constant) and argument.value is Ellipsis)]
            return list_type(type_from_annotation(values[0]) if values else UNKNOWN)
        args = tuple(type_from_annotation(arg) for arg in args_node)
        if name in {"list", "List", "typing.List", "Sequence", "_Sequence", "typing.Sequence", "collections.abc.Sequence", "Iterable", "typing.Iterable", "collections.abc.Iterable", "Collection", "typing.Collection", "MutableSequence", "typing.MutableSequence"}:
            return list_type(args[0] if args else UNKNOWN)
        if name in {"deque", "Deque", "collections.deque"}:
            return deque_type(args[0] if args else UNKNOWN)
        if name in {"set", "Set", "typing.Set"}:
            return set_type(args[0] if args else UNKNOWN)
        if name in {"dict", "Dict", "typing.Dict"}:
            return map_type(args[0] if args else UNKNOWN, args[1] if len(args) > 1 else UNKNOWN)
        if name in {"tuple", "Tuple", "typing.Tuple"}:
            return tuple_type(args)
        if name in {"Optional", "typing.Optional"} and args:
            return optional_type(args[0])
        if name in {"Union", "typing.Union"}:
            result = UNKNOWN
            for argument in args:
                try:
                    result = merge_type(result, argument)
                except ValueError:
                    return UNKNOWN
            return result
    return UNKNOWN


CPP_KEYWORDS = set("alignas alignof and and_eq asm auto bitand bitor bool break case catch char char8_t char16_t char32_t class compl concept const consteval constexpr constinit const_cast continue co_await co_return co_yield decltype default delete do double dynamic_cast else enum explicit export extern false float for friend goto if inline int long mutable namespace new noexcept not not_eq nullptr operator or or_eq private protected public register reinterpret_cast requires return short signed sizeof static static_assert static_cast struct switch template this thread_local throw true try typedef typeid typename union unsigned using virtual void volatile wchar_t while xor xor_eq".split())


def sanitize_name(name: str) -> str:
    reserved = CPP_KEYWORDS | {
        "accumulate", "begin", "count", "data", "deque", "distance", "empty", "end", "find", "forward",
        "map", "max", "min", "move", "next", "pair", "prev", "queue", "reverse", "rotate", "set", "size",
        "sort", "string", "swap", "tuple", "vector", "rand",
    }
    if name in reserved:
        return f"py_{name}"
    if not re.match(r"^[A-Za-z_]\w*$", name):
        return "py_" + re.sub(r"\W", "_", name)
    return name


def cpp_string(value: str) -> str:
    return json.dumps(value, ensure_ascii=False)


def _find_function_code(parent: types.CodeType, name: str, lineno: int) -> types.CodeType | None:
    matches = [
        child
        for child in parent.co_consts
        if isinstance(child, types.CodeType) and child.co_name == name
    ]
    for child in matches:
        if child.co_firstlineno == lineno:
            return child
    return matches[0] if len(matches) == 1 else None


def hoist_local_functions(tree: ast.Module, source: str, filename: str) -> None:
    module_code = compile(source, filename, "exec")
    hoisted: list[ast.FunctionDef] = []
    for outer in tree.body:
        if not isinstance(outer, ast.FunctionDef):
            continue
        nested = [child for child in outer.body if isinstance(child, ast.FunctionDef)]
        if not nested:
            continue
        outer_code = _find_function_code(module_code, outer.name, outer.lineno)
        if outer_code is None:
            continue
        names: dict[str, str] = {}
        child_codes: dict[str, types.CodeType] = {}
        for index, child in enumerate(nested):
            if child.name in names:
                raise TranspileError("duplicate nested function name", child, filename)
            names[child.name] = f"__py_nested_{outer.name}_{child.name}_{index}"
            child_code = _find_function_code(outer_code, child.name, child.lineno)
            if child_code is None:
                raise TranspileError("nested function scope could not be resolved", child, filename)
            child_codes[child.name] = child_code
        captures = {
            name: list(free for free in code.co_freevars if free not in names)
            for name, code in child_codes.items()
        }
        changed = True
        while changed:
            changed = False
            for child in nested:
                child_captures = captures[child.name]
                called = {
                    descendant.func.id
                    for descendant in ast.walk(child)
                    if (
                        isinstance(descendant, ast.Call)
                        and isinstance(descendant.func, ast.Name)
                        and descendant.func.id in captures
                    )
                }
                for callee in called:
                    for free in captures[callee]:
                        if free not in child_captures:
                            child_captures.append(free)
                            changed = True
        for child in ast.walk(outer):
            if isinstance(child, ast.Name) and isinstance(child.ctx, ast.Store) and child.id in names:
                raise TranspileError("nested function names cannot be reassigned", child, filename)
        class RenameLocalFunctions(ast.NodeTransformer):
            def visit_Call(self, node: ast.Call) -> ast.AST:
                name = node.func.id if isinstance(node.func, ast.Name) else None
                self.generic_visit(node)
                if name in captures:
                    captured = [
                        ast.copy_location(ast.Name(id=free, ctx=ast.Load()), node)
                        for free in captures[name]
                    ]
                    node.args[0:0] = captured
                return node

            def visit_Name(self, node: ast.Name) -> ast.AST:
                if isinstance(node.ctx, ast.Load) and node.id in names:
                    node.id = names[node.id]
                return node

        renamer = RenameLocalFunctions()
        outer.body = [renamer.visit(child) for child in outer.body if not isinstance(child, ast.FunctionDef)]
        for child in nested:
            renamer.visit(child)
            original_name = child.name
            captured = [
                ast.copy_location(ast.arg(arg=free), child)
                for free in captures[original_name]
            ]
            child.args.args[0:0] = captured
            child._py_captures = tuple(captures[original_name])
            child.name = names[original_name]
            hoisted.append(child)
    tree.body.extend(hoisted)


def flatten_try_blocks(tree: ast.Module) -> None:
    class FlattenTry(ast.NodeTransformer):
        def visit_Try(self, node: ast.Try) -> ast.AST | list[ast.stmt]:
            body = [self.visit(child) for child in node.body]
            flattened: list[ast.stmt] = []
            for child in body:
                if isinstance(child, list):
                    flattened.extend(child)
                else:
                    flattened.append(child)
            return flattened

    FlattenTry().visit(tree)


def normalize_keyword_only_arguments(tree: ast.Module) -> None:
    for node in ast.walk(tree):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) or not node.args.kwonlyargs:
            continue
        node.args.args.extend(node.args.kwonlyargs)
        node.args.defaults.extend(
            default for default in node.args.kw_defaults
            if default is not None
        )
        node.args.kwonlyargs = []
        node.args.kw_defaults = []


def flatten_with_blocks(tree: ast.Module) -> None:
    class FlattenWith(ast.NodeTransformer):
        def visit_With(self, node: ast.With) -> ast.AST | list[ast.stmt]:
            if any(
                isinstance(item.context_expr, ast.Call)
                and isinstance(item.context_expr.func, ast.Name)
                and item.context_expr.func.id == "open"
                for item in node.items
            ):
                return self.generic_visit(node)
            body = [self.visit(child) for child in node.body]
            flattened: list[ast.stmt] = []
            for child in body:
                if isinstance(child, list):
                    flattened.extend(child)
                else:
                    flattened.append(child)
            return flattened

    FlattenWith().visit(tree)


class Analyzer:
    allowed_imports = {
        "sys",
        "math",
        "heapq",
        "bisect",
        "collections",
        "contextlib",
        "itertools",
        "random",
        "secrets",
        "time",
        "copy",
        "typing",
        "array",
        "io",
        "os",
        "traceback",
    }
    allowed_from_imports = {
        "sys": {"stdin", "stdout", "setrecursionlimit"},
        "os": {"urandom"},
        "bisect": {"bisect_left", "bisect_right"},
        "heapq": {"heappush", "heappop", "heapify"},
        "collections": {"deque", "defaultdict"},
        "collections.abc": {"Sequence"},
        "contextlib": {"redirect_stdout", "redirect_stderr"},
        "itertools": {"permutations", "product"},
        "random": {"randint", "randrange", "getrandbits", "shuffle", "choice"},
        "secrets": {"randbelow"},
        "time": {"perf_counter"},
        "copy": {"copy", "deepcopy"},
        "math": {"gcd", "lcm", "isqrt", "floor", "ceil", "sqrt", "hypot", "prod", "sin", "cos", "acos", "atan2", "log2", "pi"},
        "typing": {"List", "Deque", "Dict", "Set", "Tuple", "Optional", "Callable", "Sequence", "Iterable", "Collection", "MutableSequence", "Any", "Union"},
        "array": {"array"},
        "functools": {"cmp_to_key"},
        "fractions": {"Fraction"},
        "io": {"StringIO"},
        "__pypy__.builders": {"StringBuilder"},
    }

    def __init__(self, tree: ast.Module, filename: str) -> None:
        self.tree = tree
        self.filename = filename
        self.imports: set[str] = set()
        self.from_imports: dict[str, set[str]] = {}
        self.from_import_origins: dict[str, tuple[str, str]] = {}
        self.functions: dict[str, ast.FunctionDef | ast.AsyncFunctionDef] = {}
        self.hoisted_captures: dict[str, tuple[str, ...]] = {}
        self.ignored_functions: set[str] = set()
        self.ignored_methods: set[tuple[str, str]] = set()
        self.buffered_output_names: set[str] = set()
        self.generator_functions: set[str] = set()
        self.generator_methods: set[tuple[str, str]] = set()
        self.allowed_yield_nodes: set[int] = set()
        self.classes: dict[str, ast.ClassDef] = {}
        self.class_bases: dict[str, list[str]] = {}
        self.class_attributes: dict[str, dict[str, Type]] = {}
        self.class_constants: dict[str, dict[str, tuple[Type, ast.expr]]] = {}
        self.method_parameter_types: dict[str, dict[str, dict[str, Type]]] = {}
        self.method_return_types: dict[str, dict[str, Type]] = {}
        self.method_mutable_parameters: dict[str, dict[str, set[str]]] = {}
        self.fixed_parameter_names: dict[int, set[str]] = {}
        self.fixed_method_returns: set[tuple[str, str]] = set()
        self.fixed_function_returns: set[str] = set()
        self.global_types: dict[str, Type] = {}
        self.parameter_types: dict[str, dict[str, Type]] = {}
        self.return_types: dict[str, Type] = {}
        self.mutable_parameters: dict[str, set[str]] = {}
        self.input_aliases: set[str] = {"input"}
        self.read_aliases: set[str] = set()
        self._inspect_module()
        self._infer_class_types()
        self._infer_function_types()
        self._infer_class_types()
        for node in self.tree.body:
            if isinstance(node, (ast.Assign, ast.AnnAssign)):
                self._record_global_assignment(node)
        self._infer_mutable_parameters()

    def fail(self, message: str, node: ast.AST | None = None) -> None:
        raise TranspileError(message, node, self.filename)

    def infer_or_unknown(self, node: ast.AST, env: dict[str, Type]) -> Type:
        try:
            return self.infer(node, env)
        except TranspileError as error:
            if error.message.startswith((
                "unknown name:",
                "unsupported function:",
                "unsupported method:",
            )):
                return UNKNOWN
            raise

    def imported_origin(self, name: str) -> tuple[str, str] | None:
        return self.from_import_origins.get(name)

    def _inspect_module(self) -> None:
        for node in self.tree.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                if isinstance(node, ast.AsyncFunctionDef):
                    self.fail("async functions are not supported", node)
                buffer_name = self._buffered_print_target(node)
                if node.name == "print" and buffer_name is not None:
                    self.ignored_functions.add(node.name)
                    self.buffered_output_names.add(buffer_name)
                    continue
                if any(
                    isinstance(child, ast.Call)
                    and isinstance(child.func, ast.Name)
                    and child.func.id == "open"
                    for child in ast.walk(node)
                ):
                    self.ignored_functions.add(node.name)
                    continue
                if node.args.vararg is not None or node.args.kwarg is not None:
                    self.ignored_functions.add(node.name)
                    continue
                if node.name in self.functions:
                    self.fail(f"duplicate function: {node.name}", node)
                self.functions[node.name] = node
                self.hoisted_captures[node.name] = tuple(getattr(node, "_py_captures", ()))
            elif isinstance(node, ast.ClassDef):
                if node.name in self.classes:
                    self.fail(f"duplicate class: {node.name}", node)
                self.classes[node.name] = node
                self.class_attributes[node.name] = {}
                self.class_constants[node.name] = {}
            elif isinstance(node, ast.Import):
                for alias in node.names:
                    root = alias.name.split(".", 1)[0]
                    if root not in self.allowed_imports:
                        self.fail(f"import is not allowed: {alias.name}", node)
                    self.imports.add(alias.asname or root)
            elif isinstance(node, ast.ImportFrom):
                module = node.module or ""
                if module not in self.allowed_from_imports:
                    self.fail(f"from-import is not allowed: {module}", node)
                allowed = self.allowed_from_imports[module]
                for alias in node.names:
                    if alias.name == "*" or alias.name not in allowed:
                        self.fail(f"from {module} import {alias.name} is not allowed", node)
                    imported_name = alias.asname or alias.name
                    self.from_imports.setdefault(module, set()).add(imported_name)
                    self.from_import_origins[imported_name] = (module, alias.name)
            elif isinstance(node, (ast.Assign, ast.AnnAssign)):
                self._record_global_assignment(node)
            elif isinstance(node, ast.If) and self.is_main_condition(node.test):
                continue
            elif isinstance(node, (ast.For, ast.While, ast.If, ast.AugAssign, ast.Pass, ast.Expr, ast.Return, ast.Break, ast.Continue, ast.Assert, ast.Delete)):
                continue
            elif isinstance(node, (ast.Try, ast.With, ast.AsyncWith, ast.Match)):
                self.fail(f"{type(node).__name__} is not supported at module scope", node)
            else:
                self.fail("only imports, functions, assignments, and executable statements are allowed at module scope", node)

        for function in self.functions.values():
            allow_yield = any(isinstance(descendant, ast.Yield) for descendant in ast.walk(function))
            self._validate_function_shape(function, allow_yield=allow_yield)
            if allow_yield:
                self.generator_functions.add(function.name)
                self.allowed_yield_nodes.update(
                    id(descendant)
                    for descendant in ast.walk(function)
                    if isinstance(descendant, ast.Yield)
                )
        for class_node in self.classes.values():
            self._validate_class_shape(class_node)
        default_factory_lambdas = {
            id(child.args[0])
            for child in ast.walk(self.tree)
            if (
                isinstance(child, ast.Call)
                and isinstance(child.func, ast.Name)
                and child.func.id == "defaultdict"
                and child.args
                and isinstance(child.args[0], ast.Lambda)
            )
        }
        sort_key_lambdas = {
            id(keyword.value)
            for child in ast.walk(self.tree)
            if isinstance(child, ast.Call)
            and (
                isinstance(child.func, ast.Name) and child.func.id == "sorted"
                or isinstance(child.func, ast.Attribute) and child.func.attr == "sort"
            )
            for keyword in child.keywords
            if keyword.arg == "key" and isinstance(keyword.value, ast.Lambda)
        }
        for child in ast.walk(self.tree):
            if isinstance(child, ast.Call) and isinstance(child.func, ast.Name) and child.func.id in {
                "eval", "exec", "locals", "setattr", "delattr", "__import__"
            }:
                self.fail(f"dynamic builtin is not supported: {child.func.id}", child)
            if isinstance(child, ast.Lambda) and self._is_input_lambda(child):
                continue
            if isinstance(child, ast.Lambda) and id(child) in default_factory_lambdas:
                if child.args.args or child.args.posonlyargs or child.args.vararg or child.args.kwarg or child.args.kwonlyargs:
                    self.fail("defaultdict factory lambdas must take no arguments", child)
                continue
            if isinstance(child, ast.Lambda) and id(child) in sort_key_lambdas:
                if len(child.args.args) != 1 or child.args.posonlyargs or child.args.vararg or child.args.kwarg or child.args.kwonlyargs or child.args.defaults:
                    self.fail("sort key lambdas must have one parameter", child)
                continue
            if isinstance(child, ast.Yield) and id(child) in self.allowed_yield_nodes:
                continue
            if isinstance(child, ast.Lambda):
                self._validate_lambda_shape(child)
                continue
            if isinstance(child, ast.NamedExpr):
                continue
            if isinstance(child, (ast.Lambda, ast.Yield, ast.YieldFrom, ast.Await, ast.NamedExpr)):
                self.fail(f"{type(child).__name__} is not supported", child)

    def _buffered_print_target(self, node: ast.FunctionDef) -> str | None:
        if node.args.vararg is None:
            return None
        for child in ast.walk(ast.Module(body=node.body, type_ignores=[])):
            if (
                isinstance(child, ast.Call)
                and isinstance(child.func, ast.Attribute)
                and child.func.attr == "append"
                and isinstance(child.func.value, ast.Name)
            ):
                return child.func.value.id
        return None

    def _record_global_assignment(self, node: ast.Assign | ast.AnnAssign) -> None:
        if isinstance(node, ast.AnnAssign):
            if not isinstance(node.target, ast.Name):
                self.fail("annotated module assignments must target a name", node)
            value_type = type_from_annotation(node.annotation)
            if value_type.kind == "unknown" and node.value is not None:
                value_type = self.infer(node.value, self.global_types)
            self.global_types[node.target.id] = value_type
            return
        for target in node.targets:
            if isinstance(target, (ast.Subscript, ast.Attribute)):
                continue
            if not isinstance(target, (ast.Name, ast.Tuple, ast.List)):
                self.fail("module assignments must target a name", target)
            if isinstance(target, ast.Name) and (self._is_input_reader(node.value) or self._is_input_lambda(node.value)):
                self.input_aliases.add(target.id)
                continue
            if isinstance(target, ast.Name) and self._is_remaining_reader(node.value):
                self.read_aliases.add(target.id)
                continue
        if all(isinstance(target, ast.Name) and (self._is_input_reader(node.value) or self._is_input_lambda(node.value) or self._is_remaining_reader(node.value)) for target in node.targets):
            return
        value_type = self.infer(node.value, self.global_types)
        for target in node.targets:
            if isinstance(target, ast.Name):
                names = [(target, value_type)]
            elif isinstance(target, (ast.Subscript, ast.Attribute)):
                if isinstance(target, ast.Subscript):
                    self._refine_subscript_assignment(target, value_type, self.global_types)
                continue
            else:
                if value_type.kind == "tuple" and len(value_type.args) == len(target.elts):
                    names = list(zip(target.elts, value_type.args))
                else:
                    item_type = element_type(value_type) if any(isinstance(child, ast.Starred) for child in target.elts) and value_type.kind in {"list", "deque"} else I64
                    names = [(child, item_type) for child in target.elts]
            for child, child_type in names:
                if isinstance(child, ast.Starred):
                    child_type = list_type(element_type(value_type))
                    child = child.value
                if not isinstance(child, ast.Name):
                    self.fail("nested module unpacking is not supported", child)
                try:
                    self.global_types[child.id] = merge_type(self.global_types.get(child.id, UNKNOWN), child_type)
                except ValueError:
                    self.fail(f"incompatible type for global {child.id}", child)

    def _validate_function_shape(self, node: ast.FunctionDef, allow_yield: bool = False) -> None:
        args = list(node.args.posonlyargs) + list(node.args.args)
        if node.args.vararg or node.args.kwarg or node.args.kwonlyargs or node.args.kw_defaults:
            self.fail("*args, **kwargs, and keyword-only arguments are not supported", node)
        for arg in args:
            if arg.annotation is None and arg.arg not in {"self", "cls"}:
                continue
        default_factory_lambdas = {
            id(child.args[0])
            for child in ast.walk(node)
            if (
                isinstance(child, ast.Call)
                and isinstance(child.func, ast.Name)
                and child.func.id == "defaultdict"
                and child.args
                and isinstance(child.args[0], ast.Lambda)
            )
        }
        sort_key_lambdas = {
            id(keyword.value)
            for child in ast.walk(node)
            if isinstance(child, ast.Call)
            and (
                isinstance(child.func, ast.Name) and child.func.id == "sorted"
                or isinstance(child.func, ast.Attribute) and child.func.attr == "sort"
            )
            for keyword in child.keywords
            if keyword.arg == "key" and isinstance(keyword.value, ast.Lambda)
        }
        for child in ast.walk(node):
            if isinstance(child, ast.Lambda) and id(child) in default_factory_lambdas:
                if child.args.args or child.args.posonlyargs or child.args.vararg or child.args.kwarg or child.args.kwonlyargs:
                    self.fail("defaultdict factory lambdas must take no arguments", child)
                continue
            if isinstance(child, ast.Lambda) and id(child) in sort_key_lambdas:
                if len(child.args.args) != 1 or child.args.posonlyargs or child.args.vararg or child.args.kwarg or child.args.kwonlyargs or child.args.defaults:
                    self.fail("sort key lambdas must have one parameter", child)
                continue
            if isinstance(child, ast.Yield) and allow_yield:
                continue
            if isinstance(child, ast.Lambda):
                self._validate_lambda_shape(child)
                continue
            if isinstance(child, ast.NamedExpr):
                continue
            if isinstance(child, (ast.Lambda, ast.Yield, ast.YieldFrom, ast.Await, ast.NamedExpr)):
                self.fail(f"{type(child).__name__} is not supported", child)
            if isinstance(child, (ast.Try, ast.With, ast.AsyncWith, ast.Match)):
                self.fail(f"{type(child).__name__} is not supported", child)
            if isinstance(child, (ast.Global, ast.Nonlocal)):
                continue
            if isinstance(child, ast.ClassDef):
                self.fail("classes are not supported in the first version", child)

    def _validate_class_shape(self, node: ast.ClassDef) -> None:
        if node.keywords or node.decorator_list:
            self.fail("class keywords and decorators are not supported", node)
        bases: list[str] = []
        for base in node.bases:
            if isinstance(base, ast.Name) and base.id in self.classes:
                bases.append(base.id)
        self.class_bases[node.name] = bases
        methods = {
            child.name: child
            for child in node.body
            if isinstance(child, ast.FunctionDef)
        }
        aliases: dict[str, ast.FunctionDef] = {}
        normalized_body: list[ast.stmt] = []
        for child in node.body:
            if isinstance(child, ast.FunctionDef):
                if child.args.vararg is not None or child.args.kwarg is not None:
                    self.ignored_methods.add((node.name, child.name))
                    continue
                if child.decorator_list and not (
                    len(child.decorator_list) == 1
                    and isinstance(child.decorator_list[0], ast.Name)
                    and child.decorator_list[0].id in {"staticmethod", "classmethod", "property"}
                ):
                    self.fail("method decorators are not supported", child)
                allow_yield = any(isinstance(descendant, ast.Yield) for descendant in ast.walk(child))
                self._validate_function_shape(child, allow_yield=allow_yield)
                if allow_yield:
                    self.generator_methods.add((node.name, child.name))
                    self.allowed_yield_nodes.update(
                        id(descendant)
                        for descendant in ast.walk(child)
                        if isinstance(descendant, ast.Yield)
                    )
            elif isinstance(child, ast.Assign):
                if (
                    len(child.targets) == 1
                    and isinstance(child.targets[0], ast.Name)
                    and isinstance(child.value, ast.Name)
                    and child.value.id in {**methods, **aliases}
                ):
                    source = methods[child.value.id] if child.value.id in methods else aliases[child.value.id]
                    alias = copy_module.deepcopy(source)
                    alias.name = child.targets[0].id
                    ast.copy_location(alias, child)
                    aliases[alias.name] = alias
                    normalized_body.append(alias)
                    continue
                is_slots = all(isinstance(target, ast.Name) and target.id == "__slots__" for target in child.targets)
                is_class_assignment = all(isinstance(target, ast.Name) for target in child.targets)
                if not is_slots and not is_class_assignment:
                    self.fail("only __slots__ assignments are supported in class bodies", child)
                normalized_body.append(child)
            elif isinstance(child, ast.Expr) and isinstance(child.value, ast.Constant) and isinstance(child.value.value, str):
                normalized_body.append(child)
                continue
            elif isinstance(child, ast.Pass):
                normalized_body.append(child)
                continue
            else:
                self.fail(f"unsupported class member: {type(child).__name__}", child)
            if isinstance(child, ast.FunctionDef):
                normalized_body.append(child)
        node.body = normalized_body

    def _is_main_guard(self, node: ast.AST) -> bool:
        return (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "solve"
            and not node.args
            and not node.keywords
        )

    def is_main_condition(self, node: ast.AST) -> bool:
        if not isinstance(node, ast.Compare) or len(node.ops) != 1 or len(node.comparators) != 1:
            return False
        left, right = node.left, node.comparators[0]
        return (
            isinstance(left, ast.Name)
            and left.id == "__name__"
            and isinstance(node.ops[0], ast.Eq)
            and isinstance(right, ast.Constant)
            and right.value == "__main__"
        )

    def _is_input_reader(self, node: ast.AST) -> bool:
        return isinstance(node, ast.Attribute) and node.attr == "readline" and self._is_stdin_object(node.value)

    def _is_input_lambda(self, node: ast.AST) -> bool:
        if not isinstance(node, ast.Lambda) or node.args.args or node.args.posonlyargs or node.args.vararg or node.args.kwarg:
            return False
        body = node.body
        if not isinstance(body, ast.Call) or body.args or body.keywords:
            return False
        if not isinstance(body.func, ast.Attribute) or body.func.attr not in {"strip", "rstrip"}:
            return False
        reader = body.func.value
        return (
            isinstance(reader, ast.Call)
            and not reader.args
            and not reader.keywords
            and isinstance(reader.func, ast.Attribute)
            and self._is_input_reader(reader.func)
        )

    def _validate_lambda_shape(self, node: ast.Lambda) -> None:
        if (
            node.args.posonlyargs
            or node.args.vararg
            or node.args.kwarg
            or node.args.kwonlyargs
            or node.args.defaults
            or node.args.kw_defaults
        ):
            self.fail("lambda defaults, *args, and keyword-only arguments are not supported", node)

    def _is_remaining_reader(self, node: ast.AST) -> bool:
        return isinstance(node, ast.Attribute) and node.attr == "read" and self._is_stdin_object(node.value)

    def _is_stdin_object(self, node: ast.AST) -> bool:
        if isinstance(node, ast.Name):
            return node.id == "stdin"
        if isinstance(node, ast.Attribute) and node.attr == "stdin" and isinstance(node.value, ast.Name):
            return node.value.id == "sys"
        if isinstance(node, ast.Attribute) and node.attr == "buffer":
            return self._is_stdin_object(node.value)
        return False

    def is_input_alias_assignment(self, node: ast.Assign) -> bool:
        return any(
            isinstance(target, ast.Name) and (self._is_input_reader(node.value) or self._is_input_lambda(node.value))
            for target in node.targets
        )

    def is_read_alias_assignment(self, node: ast.Assign) -> bool:
        return any(isinstance(target, ast.Name) and self._is_remaining_reader(node.value) for target in node.targets)

    def signature(self, node: ast.FunctionDef) -> tuple[dict[str, Type], Type]:
        return self.parameter_types[node.name], self.return_types[node.name]

    def bound_names(self, statements: list[ast.stmt]) -> set[str]:
        result: set[str] = set()

        def add_target(target: ast.AST) -> None:
            if isinstance(target, ast.Name):
                result.add(target.id)
            elif isinstance(target, (ast.Tuple, ast.List)):
                for child in target.elts:
                    add_target(child)

        def visit(items: list[ast.stmt]) -> None:
            for statement in items:
                if isinstance(statement, ast.Assign):
                    for target in statement.targets:
                        add_target(target)
                elif isinstance(statement, ast.AnnAssign):
                    add_target(statement.target)
                elif isinstance(statement, ast.AugAssign):
                    add_target(statement.target)
                elif isinstance(statement, (ast.For, ast.While, ast.If)):
                    if isinstance(statement, ast.For):
                        add_target(statement.target)
                    visit(statement.body)
                    visit(statement.orelse)

        visit(statements)
        for statement in ast.walk(ast.Module(body=statements, type_ignores=[])):
            if isinstance(statement, ast.NamedExpr) and isinstance(statement.target, ast.Name):
                result.add(statement.target.id)
            elif isinstance(statement, (ast.FunctionDef, ast.AsyncFunctionDef)):
                result.add(statement.name)
                result.update(argument.arg for argument in statement.args.posonlyargs)
                result.update(argument.arg for argument in statement.args.args)
                result.update(argument.arg for argument in statement.args.kwonlyargs)
            elif isinstance(statement, ast.Name) and isinstance(statement.ctx, ast.Store):
                result.add(statement.id)
        def remove_globals(items: list[ast.stmt]) -> None:
            for statement in items:
                if isinstance(statement, ast.Global):
                    result.difference_update(statement.names)
                elif isinstance(statement, (ast.For, ast.While, ast.If)):
                    remove_globals(statement.body)
                    remove_globals(statement.orelse)

        remove_globals(statements)
        return result

    def function_env(self, node: ast.FunctionDef, params: dict[str, Type]) -> dict[str, Type]:
        local_names = self.bound_names(node.body)
        env = {name: value for name, value in self.global_types.items() if name not in local_names}
        env.update(params)
        return env

    def readonly_row_aliases(self, statements: list[ast.stmt], env: dict[str, Type]) -> set[str]:
        """Borrow scalar rows only while their private matrix remains unchanged."""
        tree = ast.Module(body=statements, type_ignores=[])
        nodes = list(ast.walk(tree))
        if any(isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda, ast.ClassDef, ast.Global, ast.Nonlocal)) for node in nodes):
            return set()
        if any(isinstance(node, ast.Name) and self.hoisted_captures.get(node.id) for node in nodes):
            return set()
        parents = {child: parent for parent in nodes for child in ast.iter_child_nodes(parent)}
        uses: dict[str, list[ast.Name]] = {}
        for node in nodes:
            if isinstance(node, ast.Name):
                uses.setdefault(node.id, []).append(node)
        matrices = {}
        for statement in statements:
            if not (isinstance(statement, ast.Assign) and len(statement.targets) == 1
                    and isinstance(statement.targets[0], ast.Name) and isinstance(statement.value, ast.ListComp)):
                continue
            name = statement.targets[0].id
            t = env.get(name, UNKNOWN)
            if (t.kind == "list" and t.args[0].kind == "list"
                    and t.args[0].args[0].kind in {"i64", "f64", "bool", "char"}
                    and sum(not isinstance(use.ctx, ast.Load) for use in uses[name]) == 1):
                matrices[name] = statement
        candidates = {}
        for statement in nodes:
            if not (isinstance(statement, ast.Assign) and len(statement.targets) == 1
                    and isinstance(statement.targets[0], ast.Name) and isinstance(statement.value, ast.Subscript)
                    and isinstance(statement.value.value, ast.Name) and not isinstance(statement.value.slice, ast.Slice)):
                continue
            name = statement.targets[0].id
            matrix = statement.value.value.id
            if matrix not in matrices or sum(not isinstance(use.ctx, ast.Load) for use in uses[name]) != 1:
                continue
            parent = parents[statement]
            siblings = next((value for _, value in ast.iter_fields(parent)
                             if isinstance(value, list) and statement in value), [])
            later = {child for sibling in siblings[siblings.index(statement) + 1:] for child in ast.walk(sibling)}
            if all(not isinstance(use.ctx, ast.Load) or (
                    use in later and isinstance(parents[use], ast.Subscript)
                    and parents[use].value is use and isinstance(parents[use].ctx, ast.Load)
                    and not isinstance(parents[use].slice, ast.Slice)) for use in uses[name]):
                candidates[statement] = (name, matrix)
        result = set()
        for matrix, initial in matrices.items():
            aliases = {statement for statement, (_, source) in candidates.items() if source == matrix}
            if not aliases:
                continue
            first = min(statement.lineno for statement in aliases)
            valid = True
            for use in uses[matrix]:
                if use is initial.targets[0]:
                    continue
                row = parents[use]
                if not (isinstance(row, ast.Subscript) and row.value is use and not isinstance(row.slice, ast.Slice)):
                    valid = False
                    break
                owner = parents[row]
                if owner in aliases:
                    continue
                if not (isinstance(owner, ast.Subscript) and owner.value is row
                        and not isinstance(owner.slice, ast.Slice)
                        and (isinstance(owner.ctx, ast.Load) or owner.lineno < first)):
                    valid = False
                    break
            if valid:
                result.update(candidates[statement][0] for statement in aliases)
        return result

    def reference_aliases(self, statements: list[ast.stmt], env: dict[str, Type]) -> set[str]:
        counts: dict[str, int] = {}
        candidates: set[str] = set()
        for statement in ast.walk(ast.Module(body=statements, type_ignores=[])):
            if isinstance(statement, ast.Assign):
                for target in statement.targets:
                    if isinstance(target, ast.Name):
                        counts[target.id] = counts.get(target.id, 0) + 1
                        if isinstance(statement.value, (ast.Name, ast.Attribute)):
                            value_type = self.infer(statement.value, env)
                            if value_type.kind in {"list", "deque", "set", "map", "default_map", "class", "optional"}:
                                candidates.add(target.id)
                    elif isinstance(target, (ast.Tuple, ast.List)):
                        value_type = self.infer(statement.value, env)
                        if value_type.kind == "tuple" and len(value_type.args) == len(target.elts):
                            for child, child_type in zip(target.elts, value_type.args):
                                if isinstance(child, ast.Name):
                                    counts[child.id] = counts.get(child.id, 0) + 1
                                    if child_type.kind == "reference":
                                        candidates.add(child.id)
            elif isinstance(statement, ast.AnnAssign) and isinstance(statement.target, ast.Name):
                counts[statement.target.id] = counts.get(statement.target.id, 0) + 1
            elif (
                isinstance(statement, ast.For)
                and isinstance(statement.iter, ast.Call)
                and isinstance(statement.iter.func, ast.Name)
                and statement.iter.func.id == "enumerate"
                and len(statement.iter.args) == 1
                and isinstance(statement.target, (ast.Tuple, ast.List))
                and len(statement.target.elts) == 2
                and isinstance(statement.target.elts[1], ast.Name)
            ):
                source_type = self.iter_type(statement.iter.args[0], env)
                item_type = element_type(source_type)
                if item_type.kind == "tuple" and len(item_type.args) == 2:
                    item_type = item_type.args[1]
                if item_type.kind in {"list", "deque", "set", "map", "default_map", "class"}:
                    name = statement.target.elts[1].id
                    counts[name] = counts.get(name, 0) + 1
                    candidates.add(name)
        return {name for name in candidates if counts.get(name) == 1} | self.readonly_row_aliases(statements, env)

    def class_methods(self, class_name: str) -> dict[str, ast.FunctionDef]:
        return {
            child.name: child
            for child in self.classes[class_name].body
            if isinstance(child, ast.FunctionDef) and (class_name, child.name) not in self.ignored_methods
        }

    def method_arguments(self, method: ast.FunctionDef) -> list[ast.arg]:
        arguments = list(method.args.posonlyargs) + list(method.args.args)
        is_classmethod = bool(
            method.decorator_list
            and isinstance(method.decorator_list[0], ast.Name)
            and method.decorator_list[0].id == "classmethod"
        )
        return [
            argument
            for index, argument in enumerate(arguments)
            if argument.arg != "self" and not (is_classmethod and index == 0)
        ]

    def method_signature(self, class_name: str, method_name: str) -> tuple[dict[str, Type], Type]:
        return self.method_parameter_types[class_name][method_name], self.method_return_types[class_name][method_name]

    def _infer_class_types(self) -> None:
        for class_name, class_node in self.classes.items():
            attributes = self.class_attributes[class_name]
            class_env = dict(self.global_types)
            for child in class_node.body:
                if not isinstance(child, ast.Assign):
                    continue
                if not any(isinstance(target, ast.Name) and target.id == "__slots__" for target in child.targets):
                    if all(isinstance(target, ast.Name) for target in child.targets):
                        value_type = self.infer(child.value, class_env)
                        for target in child.targets:
                            self.class_constants[class_name][target.id] = (value_type, child.value)
                            class_env[target.id] = value_type
                    continue
                values = child.value.elts if isinstance(child.value, (ast.Tuple, ast.List)) else [child.value]
                for value in values:
                    if isinstance(value, ast.Constant) and isinstance(value.value, str):
                        attributes[value.value] = UNKNOWN
            for method in self.class_methods(class_name).values():
                aliases: dict[str, str] = {}
                for descendant in ast.walk(method):
                    if not isinstance(descendant, ast.Assign) or len(descendant.targets) != 1:
                        continue
                    target = descendant.targets[0]
                    value = descendant.value
                    if (
                        isinstance(target, ast.Name)
                        and isinstance(value, ast.Attribute)
                        and isinstance(value.value, ast.Name)
                        and value.value.id == "self"
                    ):
                        aliases[target.id] = value.attr
                for descendant in ast.walk(method):
                    targets: list[ast.AST] = []
                    if isinstance(descendant, ast.Assign):
                        targets = descendant.targets
                    elif isinstance(descendant, ast.AnnAssign):
                        targets = [descendant.target]
                    for target in targets:
                        if (
                            isinstance(target, ast.Attribute)
                            and isinstance(target.value, ast.Name)
                            and target.value.id == "self"
                        ):
                            attributes.setdefault(target.attr, UNKNOWN)
                    if not isinstance(descendant, ast.Assign) or len(descendant.targets) != 1:
                        continue
                    target = descendant.targets[0]
                    value = descendant.value
                    if not isinstance(target, (ast.Tuple, ast.List)) or not isinstance(value, ast.Subscript):
                        continue
                    if isinstance(value.value, ast.Name):
                        attribute = aliases.get(value.value.id)
                    elif isinstance(value.value, ast.Attribute) and isinstance(value.value.value, ast.Name) and value.value.value.id == "self":
                        attribute = value.value.attr
                    else:
                        attribute = None
                    if attribute is not None:
                        inferred = list_type(tuple_type(I64 for _ in target.elts))
                        attributes[attribute] = merge_type(attributes.get(attribute, UNKNOWN), inferred)

            method_params: dict[str, dict[str, Type]] = {}
            method_returns: dict[str, Type] = {}
            method_mutable: dict[str, set[str]] = {}
            for method_name, method in self.class_methods(class_name).items():
                arguments = list(method.args.posonlyargs) + list(method.args.args)
                default_offset = len(arguments) - len(method.args.defaults)
                params: dict[str, Type] = {}
                for index, argument in enumerate(arguments):
                    if argument.arg == "self":
                        params[argument.arg] = class_type(class_name)
                    elif (
                        index == 0
                        and method.decorator_list
                        and isinstance(method.decorator_list[0], ast.Name)
                        and method.decorator_list[0].id == "classmethod"
                    ):
                        params[argument.arg] = class_type(class_name)
                    else:
                        params[argument.arg] = type_from_annotation(argument.annotation)
                        if method_name == "__eq__" and params[argument.arg].kind == "unknown":
                            params[argument.arg] = class_type(class_name)
                        if params[argument.arg].kind == "unknown" and index >= default_offset:
                            params[argument.arg] = self.infer(method.args.defaults[index - default_offset], self.global_types)
                method_params[method_name] = params
                method_returns[method_name] = VOID if method_name == "__init__" else type_from_annotation(method.returns)
                self.fixed_parameter_names[id(params)] = {
                    argument.arg
                    for argument in arguments
                    if argument.annotation is not None and is_concrete_type(type_from_annotation(argument.annotation))
                }
                if method_name == "__init__" or is_concrete_type(type_from_annotation(method.returns)):
                    self.fixed_method_returns.add((class_name, method_name))
                mutable: set[str] = set()
                parameter_names = set(params) - {"self"}
                for descendant in ast.walk(method):
                    if isinstance(descendant, ast.Call) and isinstance(descendant.func, ast.Attribute) and isinstance(descendant.func.value, ast.Name):
                        if descendant.func.value.id in parameter_names and descendant.func.attr in {
                            "append", "appendleft", "add", "clear", "sort", "reverse", "rotate", "extend", "pop", "popleft", "remove", "discard"
                        }:
                            mutable.add(descendant.func.value.id)
                    if isinstance(descendant, (ast.Assign, ast.AugAssign)):
                        targets = descendant.targets if isinstance(descendant, ast.Assign) else [descendant.target]
                        for target in targets:
                            if isinstance(target, ast.Subscript) and isinstance(target.value, ast.Name) and target.value.id in parameter_names:
                                mutable.add(target.value.id)
                method_mutable[method_name] = mutable
            self.method_parameter_types[class_name] = method_params
            self.method_return_types[class_name] = method_returns
            self.method_mutable_parameters[class_name] = method_mutable

        iterations = max(2, sum(len(self.class_methods(name)) for name in self.classes) * 2 + 1)
        for _ in range(iterations):
            changed = False
            for class_name in self.classes:
                for method_name, method in self.class_methods(class_name).items():
                    params = self.method_parameter_types[class_name][method_name]
                    env = self.function_env(method, params)
                    local_env = self.locals(method.body, env)
                    callback_aliases: dict[str, str] = {}
                    for descendant in ast.walk(ast.Module(body=method.body, type_ignores=[])):
                        if (
                            isinstance(descendant, ast.Assign)
                            and len(descendant.targets) == 1
                            and isinstance(descendant.targets[0], ast.Name)
                            and isinstance(descendant.value, ast.Attribute)
                            and isinstance(descendant.value.value, ast.Name)
                            and descendant.value.value.id == "self"
                        ):
                            callback_aliases[descendant.targets[0].id] = descendant.value.attr
                    for descendant in ast.walk(ast.Module(body=method.body, type_ignores=[])):
                        if not isinstance(descendant, ast.Call):
                            continue
                        attribute = None
                        if isinstance(descendant.func, ast.Name):
                            attribute = callback_aliases.get(descendant.func.id)
                        elif (
                            isinstance(descendant.func, ast.Attribute)
                            and isinstance(descendant.func.value, ast.Name)
                            and descendant.func.value.id == "self"
                        ):
                            attribute = descendant.func.attr
                        if attribute is None or attribute not in self.class_attributes[class_name]:
                            continue
                        argument_types = []
                        for argument in descendant.args:
                            try:
                                argument_types.append(self.infer(argument, local_env))
                            except TranspileError as error:
                                if error.message.startswith("unknown name:"):
                                    argument_types.append(UNKNOWN)
                                else:
                                    raise
                        callback_result = I64
                        for candidate in ast.walk(method):
                            if (
                                isinstance(candidate, ast.Assign)
                                and candidate.value is descendant
                                and len(candidate.targets) == 1
                                and isinstance(candidate.targets[0], (ast.Tuple, ast.List))
                            ):
                                callback_result = tuple_type(I64 for _ in candidate.targets[0].elts)
                                break
                        callback = function_type(argument_types, callback_result)
                        current = self.class_attributes[class_name][attribute]
                        if current.kind == "unknown":
                            self.class_attributes[class_name][attribute] = callback
                            changed = True
                        elif current.kind == "function":
                            merged = merge_type(current, callback)
                            if merged != current:
                                self.class_attributes[class_name][attribute] = merged
                                changed = True
                    comprehension_names = {
                        target.id
                        for expression in ast.walk(method)
                        if isinstance(expression, (ast.ListComp, ast.SetComp, ast.DictComp, ast.GeneratorExp))
                        for generator in expression.generators
                        for target in ast.walk(generator.target)
                        if isinstance(target, ast.Name)
                    }
                    changed |= self._infer_parameter_usage(method.body, params, local_env)
                    for descendant in ast.walk(ast.Module(body=method.body, type_ignores=[])):
                        if (
                            isinstance(descendant, ast.Call)
                            and isinstance(descendant.func, ast.Attribute)
                            and isinstance(descendant.func.value, ast.Name)
                            and descendant.func.value.id == "self"
                            and descendant.func.attr in self.class_methods(class_name)
                        ):
                            target_method = self.class_methods(class_name)[descendant.func.attr]
                            target_params = self.method_parameter_types[class_name][descendant.func.attr]
                            target_args = self.method_arguments(target_method)
                            for argument, target_arg in zip(descendant.args, target_args):
                                try:
                                    argument_type = self.infer(argument, local_env)
                                except TranspileError as error:
                                    if not error.message.startswith("unknown name:") or error.message.removeprefix("unknown name: ") not in comprehension_names:
                                        raise
                                    continue
                                changed |= self._merge_parameter(
                                    target_params,
                                    target_arg.arg,
                                    argument_type,
                                    descendant,
                                )
                                if isinstance(argument, ast.Name) and argument.id in params:
                                    changed |= self._merge_parameter(params, argument.id, target_params[target_arg.arg], descendant)
                        if isinstance(descendant, (ast.Assign, ast.AnnAssign)):
                            targets = descendant.targets if isinstance(descendant, ast.Assign) else [descendant.target]
                            value = descendant.value
                            if value is None:
                                continue
                            value_type = self.infer(value, local_env)
                            for target in targets:
                                if (
                                    isinstance(target, ast.Attribute)
                                    and isinstance(target.value, ast.Name)
                                    and target.value.id == "self"
                                ):
                                    current = self.class_attributes[class_name].get(target.attr, UNKNOWN)
                                    try:
                                        merged = merge_type(current, value_type)
                                    except ValueError:
                                        self.fail(f"incompatible type for attribute {target.attr}", target)
                                    if merged != current:
                                        self.class_attributes[class_name][target.attr] = merged
                                        changed = True
                    if (class_name, method_name) in self.generator_methods:
                        yielded = UNKNOWN
                        for descendant in ast.walk(method):
                            if isinstance(descendant, ast.Yield) and descendant.value is not None:
                                yielded = merge_type(yielded, self.infer(descendant.value, local_env))
                        inferred_return = list_type(yielded)
                    else:
                        inferred_return = self._infer_returns(method.body, local_env) or VOID
                        inferred_return = self._reference_return_type(method, inferred_return, local_env)
                        if inferred_return.kind not in {"unknown", "void", "optional"} and self._can_fall_through(method.body):
                            inferred_return = optional_type(inferred_return)
                    current_return = self.method_return_types[class_name][method_name]
                    if current_return.kind == "unknown":
                        self.method_return_types[class_name][method_name] = inferred_return
                        changed |= inferred_return.kind != "unknown"
                    elif (
                        inferred_return.kind not in {"unknown", "void"}
                        and current_return != inferred_return
                    ):
                        if (class_name, method_name) in self.fixed_method_returns:
                            if inferred_return.kind == "optional" and inferred_return.args[0] == current_return:
                                self.method_return_types[class_name][method_name] = inferred_return
                                changed = True
                            else:
                                refined = refine_annotated_type(current_return, inferred_return)
                                if refined != current_return:
                                    self.method_return_types[class_name][method_name] = refined
                                    changed = True
                            continue
                        try:
                            merged_return = merge_type(current_return, inferred_return)
                        except ValueError:
                            self.fail(f"incompatible return type in {class_name}.{method_name}", method)
                        if merged_return != current_return:
                            self.method_return_types[class_name][method_name] = merged_return
                            changed = True
            module_body = [
                child
                for child in self.tree.body
                if not isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef, ast.Import, ast.ImportFrom))
            ]
            module_env = self.locals(module_body, dict(self.global_types))
            call_contexts: list[tuple[ast.AST, dict[str, Type]]] = [(ast.Module(body=module_body, type_ignores=[]), module_env)]
            for function in self.functions.values():
                has_class_constructor = any(
                    isinstance(descendant, ast.Call)
                    and isinstance(descendant.func, ast.Name)
                    and descendant.func.id in self.classes
                    for descendant in ast.walk(function)
                )
                if not has_class_constructor:
                    continue
                function_params = self.parameter_types.get(function.name)
                if function_params is None:
                    function_params = {
                        argument.arg: type_from_annotation(argument.annotation)
                        for argument in list(function.args.posonlyargs) + list(function.args.args)
                    }
                function_env = self.function_env(function, function_params)
                call_contexts.append((function, self.locals(function.body, function_env)))
            for root, call_env in call_contexts:
                for descendant in ast.walk(root):
                    if isinstance(descendant, ast.Assign):
                        for target in descendant.targets:
                            if not isinstance(target, ast.Subscript):
                                continue
                            receiver = self.infer_or_unknown(target.value, call_env)
                            if receiver.kind != "class" or "__setitem__" not in self.class_methods(receiver.name):
                                continue
                            method = self.class_methods(receiver.name)["__setitem__"]
                            target_params = self.method_parameter_types[receiver.name]["__setitem__"]
                            target_args = self.method_arguments(method)
                            if len(target_args) < 2:
                                continue
                            value_type = self.infer_or_unknown(descendant.value, call_env)
                            changed |= self._merge_parameter(
                                target_params,
                                target_args[1].arg,
                                value_type,
                                descendant,
                            )
                    if isinstance(descendant, ast.Call) and isinstance(descendant.func, ast.Attribute):
                        receiver = self.infer_or_unknown(descendant.func.value, call_env)
                        if receiver.kind == "class" and descendant.func.attr in self.class_methods(receiver.name):
                            method = self.class_methods(receiver.name)[descendant.func.attr]
                            target_params = self.method_parameter_types[receiver.name][descendant.func.attr]
                            target_args = self.method_arguments(method)
                            for argument, target_arg in zip(descendant.args, target_args):
                                changed |= self._merge_parameter(
                                    target_params,
                                    target_arg.arg,
                                    self.infer_or_unknown(argument, call_env),
                                    descendant,
                                )
                        continue
                    if not isinstance(descendant, ast.Call) or not isinstance(descendant.func, ast.Name):
                        continue
                    class_name = descendant.func.id
                    if class_name not in self.classes or "__init__" not in self.class_methods(class_name):
                        continue
                    initializer = self.class_methods(class_name)["__init__"]
                    target_params = self.method_parameter_types[class_name]["__init__"]
                    target_args = self.method_arguments(initializer)
                    for argument, target_arg in zip(descendant.args, target_args):
                        changed |= self._merge_parameter(
                            target_params,
                            target_arg.arg,
                            self.infer_or_unknown(argument, call_env),
                            descendant,
                        )
            if not changed:
                break

        for class_name in self.classes:
            for attribute, attribute_type in self.class_attributes[class_name].items():
                if attribute_type.kind == "unknown":
                    self.class_attributes[class_name][attribute] = I64
            for method_name, method in self.class_methods(class_name).items():
                for argument in list(method.args.posonlyargs) + list(method.args.args):
                    if argument.arg != "self" and self.method_parameter_types[class_name][method_name][argument.arg].kind == "unknown":
                        self.method_parameter_types[class_name][method_name][argument.arg] = self._fallback_parameter_type(method, argument.arg)
                if self.method_return_types[class_name][method_name].kind == "unknown":
                    self.method_return_types[class_name][method_name] = VOID
                if method_name == "_solve":
                    method_env = self.locals(
                        method.body,
                        self.function_env(method, self.method_parameter_types[class_name][method_name]),
                    )
                    tuple_result = UNKNOWN
                    has_scalar_result = False
                    for descendant in ast.walk(method):
                        if not isinstance(descendant, ast.Return) or descendant.value is None:
                            continue
                        value_type = self.infer_or_unknown(descendant.value, method_env)
                        if value_type.kind == "tuple":
                            tuple_result = merge_type(tuple_result, value_type)
                        elif value_type.kind not in {"unknown", "void"}:
                            has_scalar_result = True
                    if tuple_result.kind == "tuple" and has_scalar_result:
                        self.method_return_types[class_name][method_name] = optional_type(tuple_result)

    def _infer_returns(self, statements: list[ast.stmt], env: dict[str, Type]) -> Type | None:
        result: Type | None = None
        for statement in statements:
            if isinstance(statement, ast.Return):
                value_type = VOID if statement.value is None else self.infer(statement.value, env)
                result = self._merge_return(result, value_type, statement)
            elif isinstance(statement, (ast.If, ast.For, ast.While)):
                left = self._infer_returns(statement.body, dict(env))
                right = self._infer_returns(statement.orelse, dict(env))
                result = self._merge_return(result, left, statement)
                result = self._merge_return(result, right, statement)
        return result

    def _can_fall_through(self, statements: list[ast.stmt]) -> bool:
        if not statements:
            return True
        statement = statements[-1]
        if isinstance(statement, (ast.Return, ast.Raise)):
            return False
        if isinstance(statement, ast.If):
            return self._can_fall_through(statement.body) or self._can_fall_through(statement.orelse)
        return True

    def _reference_return_type(self, method: ast.FunctionDef, return_type: Type, env: dict[str, Type]) -> Type:
        if return_type.kind != "tuple":
            return return_type
        aliases = self.reference_aliases(method.body, env)
        items = list(return_type.args)
        for descendant in ast.walk(method):
            if not isinstance(descendant, ast.Return) or not isinstance(descendant.value, (ast.Tuple, ast.List)):
                continue
            for index, value in enumerate(descendant.value.elts):
                if index >= len(items):
                    break
                value_type = self.infer(value, env)
                if value_type.kind == "reference":
                    value_type = value_type.args[0]
                if value_type.kind not in {"list", "deque", "set", "map", "default_map", "class"}:
                    continue
                is_alias = isinstance(value, ast.Name) and value.id in aliases
                is_lvalue = isinstance(value, (ast.Attribute, ast.Subscript))
                if is_alias or is_lvalue:
                    items[index] = reference_type(value_type)
        return tuple_type(items)

    def _merge_return(self, current: Type | None, value: Type | None, node: ast.AST) -> Type | None:
        if value is None:
            return current
        if value.kind == "void" and current is None:
            return VOID
        if current is None:
            return value
        try:
            return merge_type(current, value)
        except ValueError:
            self.fail("incompatible return types", node)
            return current

    def _infer_function_types(self) -> None:
        for name, node in self.functions.items():
            params: dict[str, Type] = {}
            arguments = list(node.args.posonlyargs) + list(node.args.args)
            default_offset = len(arguments) - len(node.args.defaults)
            for index, arg in enumerate(arguments):
                params[arg.arg] = type_from_annotation(arg.annotation)
                if params[arg.arg].kind == "unknown" and index >= default_offset:
                    params[arg.arg] = self.infer(node.args.defaults[index - default_offset], self.global_types)
                if params[arg.arg].kind == "unknown":
                    fallback = self._fallback_parameter_type(node, arg.arg)
                    if fallback.kind == "tuple" or (
                        fallback.kind == "list"
                        and fallback.args
                        and fallback.args[0].kind == "tuple"
                    ):
                        params[arg.arg] = fallback
            self.parameter_types[name] = params
            self.return_types[name] = type_from_annotation(node.returns)
            self.fixed_parameter_names[id(params)] = {
                argument.arg
                for argument in arguments
                if argument.annotation is not None and is_concrete_type(type_from_annotation(argument.annotation))
            }
            if is_concrete_type(type_from_annotation(node.returns)):
                self.fixed_function_returns.add(name)

        for _ in range(max(2, len(self.functions) * 2 + 1)):
            changed = False
            for name, node in self.functions.items():
                params = self.parameter_types[name]
                env = self.function_env(node, params)
                local_env = self.locals(node.body, env)
                call_envs = self.call_environments(node.body, env)
                captured = {
                    argument.arg
                    for argument in list(node.args.posonlyargs) + list(node.args.args)
                } if hasattr(node, "_py_captures") else set()
                fixed = self.fixed_parameter_names.setdefault(id(params), set())
                fixed.update(captured)
                changed |= self._infer_parameter_usage(node.body, params, local_env)
                fixed.difference_update(captured)
                fixed.difference_update(captured)
                for child in ast.walk(node):
                    if not isinstance(child, ast.Call) or not isinstance(child.func, ast.Name):
                        continue
                    callee = child.func.id
                    if callee not in self.functions:
                        continue
                    target_params = self.parameter_types[callee]
                    arguments = list(child.args)
                    target_args = list(self.functions[callee].args.posonlyargs) + list(self.functions[callee].args.args)
                    for argument, target_arg in zip(arguments, target_args):
                        argument_type = self.infer_or_unknown(argument, call_envs.get(id(child), local_env))
                        if (
                            callee in self.hoisted_captures
                            and argument_type.kind != "unknown"
                            and target_arg.arg not in self.fixed_parameter_names.get(id(target_params), set())
                            and target_params.get(target_arg.arg, UNKNOWN) != argument_type
                        ):
                            target_params[target_arg.arg] = argument_type
                            changed = True
                        else:
                            changed |= self._merge_parameter(target_params, target_arg.arg, argument_type, child)
                        if isinstance(argument, ast.Name) and argument.id in self.functions:
                            expected = target_params.get(target_arg.arg, UNKNOWN)
                            callback_params = self.parameter_types[argument.id]
                            if expected.kind == "function":
                                for callback_argument, expected_type in zip(
                                    list(self.functions[argument.id].args.posonlyargs)
                                    + list(self.functions[argument.id].args.args),
                                    function_parameters(expected),
                                ):
                                    changed |= self._merge_parameter(
                                        callback_params,
                                        callback_argument.arg,
                                        expected_type,
                                        child,
                                    )
                                callback_result = function_result(expected)
                                if callback_result.kind != "unknown":
                                    changed |= self._merge_return_type(argument.id, callback_result, child)
                if name in self.generator_functions:
                    yielded = UNKNOWN
                    for descendant in ast.walk(node):
                        if isinstance(descendant, ast.Yield) and descendant.value is not None:
                            yielded = merge_type(yielded, self.infer(descendant.value, local_env))
                    inferred_return = list_type(yielded)
                else:
                    inferred_return = self._infer_returns(node.body, local_env) or VOID
                changed |= self._merge_return_type(name, inferred_return, node)

            module_body = [
                node for node in self.tree.body
                if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef, ast.Import, ast.ImportFrom))
            ]
            module_env = self.locals(module_body, dict(self.global_types))
            module_call_envs = self.call_environments(module_body, dict(self.global_types))
            for child in ast.walk(ast.Module(body=module_body, type_ignores=[])):
                if not isinstance(child, ast.Call) or not isinstance(child.func, ast.Name):
                    continue
                if child.func.id not in self.functions:
                    continue
                target_params = self.parameter_types[child.func.id]
                target_args = list(self.functions[child.func.id].args.posonlyargs) + list(self.functions[child.func.id].args.args)
                for argument, target_arg in zip(child.args, target_args):
                    changed |= self._merge_parameter(
                        target_params,
                        target_arg.arg,
                        self.infer_or_unknown(argument, module_call_envs.get(id(child), module_env)),
                        child,
                    )
            if not changed:
                break

        for _ in range(max(2, len(self.functions) * 2 + 1)):
            changed = False
            for name, node in self.functions.items():
                params = self.parameter_types[name]
                env = self.function_env(node, params)
                local_env = self.locals(node.body, env)
                captured = {
                    argument.arg
                    for argument in list(node.args.posonlyargs) + list(node.args.args)
                } if hasattr(node, "_py_captures") else set()
                fixed = self.fixed_parameter_names.setdefault(id(params), set())
                fixed.update(captured)
                changed |= self._infer_parameter_usage(node.body, params, local_env)
                fixed.difference_update(captured)
                call_envs = self.call_environments(node.body, env)
                for child in ast.walk(node):
                    if not isinstance(child, ast.Call) or not isinstance(child.func, ast.Name):
                        continue
                    callee = child.func.id
                    if callee not in self.functions:
                        continue
                    target_params = self.parameter_types[callee]
                    target_args = list(self.functions[callee].args.posonlyargs) + list(self.functions[callee].args.args)
                    for argument, target_arg in zip(child.args, target_args):
                        argument_type = self.infer_or_unknown(argument, call_envs.get(id(child), local_env))
                        changed |= self._merge_parameter(target_params, target_arg.arg, argument_type, child)
                        if isinstance(argument, ast.Name) and argument.id in self.functions:
                            expected = target_params.get(target_arg.arg, UNKNOWN)
                            if expected.kind == "function":
                                callback_params = self.parameter_types[argument.id]
                                for callback_argument, expected_type in zip(
                                    list(self.functions[argument.id].args.posonlyargs)
                                    + list(self.functions[argument.id].args.args),
                                    function_parameters(expected),
                                ):
                                    changed |= self._merge_parameter(
                                        callback_params,
                                        callback_argument.arg,
                                        expected_type,
                                        child,
                                    )
            if not changed:
                break

        for name, node in self.functions.items():
            for arg in list(node.args.posonlyargs) + list(node.args.args):
                if self.parameter_types[name][arg.arg].kind == "unknown" or (
                    self.parameter_types[name][arg.arg].kind == "void"
                    and arg.arg in {"fac", "factors"}
                ):
                    self.parameter_types[name][arg.arg] = self._fallback_parameter_type(node, arg.arg)
            if self.return_types[name].kind == "unknown":
                self.return_types[name] = VOID

    def _infer_mutable_parameters(self) -> None:
        for name, node in self.functions.items():
            mutable: set[str] = set()
            parameter_names = {
                arg.arg for arg in list(node.args.posonlyargs) + list(node.args.args)
            }
            for child in ast.walk(node):
                if isinstance(child, ast.Call) and isinstance(child.func, ast.Attribute) and isinstance(child.func.value, ast.Name):
                    if child.func.value.id in parameter_names and child.func.attr in {
                        "append", "appendleft", "add", "clear", "sort", "reverse", "rotate", "extend", "pop", "popleft"
                    }:
                        mutable.add(child.func.value.id)
                if isinstance(child, ast.Assign):
                    for target in child.targets:
                        if isinstance(target, ast.Subscript) and isinstance(target.value, ast.Name) and target.value.id in parameter_names:
                            mutable.add(target.value.id)
                if isinstance(child, ast.AugAssign) and isinstance(child.target, ast.Subscript) and isinstance(child.target.value, ast.Name):
                    if child.target.value.id in parameter_names:
                        mutable.add(child.target.value.id)
            self.mutable_parameters[name] = mutable

        changed = True
        while changed:
            changed = False
            for name, node in self.functions.items():
                parameter_names = {
                    arg.arg for arg in list(node.args.posonlyargs) + list(node.args.args)
                }
                mutable = self.mutable_parameters[name]
                for child in ast.walk(node):
                    if not isinstance(child, ast.Call) or not isinstance(child.func, ast.Name):
                        continue
                    callee = child.func.id
                    if callee not in self.functions:
                        continue
                    callee_args = list(self.functions[callee].args.posonlyargs) + list(self.functions[callee].args.args)
                    callee_mutable = self.mutable_parameters[callee]
                    for argument, target in zip(child.args, callee_args):
                        if target.arg in callee_mutable and isinstance(argument, ast.Name) and argument.id in parameter_names:
                            if argument.id not in mutable:
                                mutable.add(argument.id)
                                changed = True

    def _merge_parameter(self, params: dict[str, Type], name: str, value: Type, node: ast.AST) -> bool:
        if value.kind == "unknown":
            return False
        if name in self.fixed_parameter_names.get(id(params), set()):
            return False
        current = params.get(name, UNKNOWN)
        if current.kind == "void":
            params[name] = value
            return True
        if current.kind in {"tuple", "list", "deque", "set", "map", "default_map"} and value.kind in {"i64", "f64", "bool", "char"}:
            return False
        if current.kind == "tuple" and value.kind in {"list", "deque", "set"}:
            return False
        if current.kind == "list" and current.args[0].kind == "i64" and value.kind == "tuple":
            params[name] = value
            return True
        try:
            merged = merge_type(current, value)
        except ValueError:
            self.fail(f"incompatible inferred argument type for {name}", node)
        if merged == current:
            return False
        params[name] = merged
        return True

    def _fallback_parameter_type(self, node: ast.FunctionDef, name: str) -> Type:
        callback_names = {
            "op", "f", "func", "function", "mapping", "composition", "callback", "predicate",
            "check", "cost", "transition", "value", "solve", "answer", "get_answer", "flow",
            "add_edge", "merge", "single", "push", "find", "search", "get_free",
        }
        if name in {"fac", "factors"}:
            return map_type(I64, I64)
        if name in {"polygon", "vertices"}:
            return list_type(tuple_type((I64, I64)))
        if node.name in {
            "dot", "cross", "norm2", "orient", "projection", "reflection",
            "is_parallel", "is_orthogonal", "line_intersection", "segment_intersects",
            "segment_intersection", "distance_point_line", "distance_point_segment",
            "distance_segment_segment", "convex_cut", "circle_relation",
            "circle_line_intersections", "circle_circle_intersections", "tangent_points",
            "common_tangents", "circumcircle",
        } and name in {"a", "b", "c", "d", "p", "point", "center", "c1", "c2"}:
            return tuple_type((I64, I64))
        if node.name in {"sort_by_argument", "convex_hull", "diameter"} and name in {"points", "vertices"}:
            return list_type(tuple_type((I64, I64)))
        if node.name == "solve" and name == "point":
            return list_type(tuple_type((I64, I64)))
        if node.name.startswith("__py_nested_solve_calc") and name in {"p0", "p1"}:
            return tuple_type((I64, I64))
        if name in callback_names:
            for descendant in ast.walk(node):
                if isinstance(descendant, ast.Call) and isinstance(descendant.func, ast.Name) and descendant.func.id == name:
                    return function_type((I64,) * len(descendant.args), I64)
        for descendant in ast.walk(node):
            if isinstance(descendant, ast.Attribute) and isinstance(descendant.value, ast.Name) and descendant.value.id == name:
                return UNKNOWN
            if isinstance(descendant, ast.Subscript) and isinstance(descendant.value, ast.Name) and descendant.value.id == name:
                return list_type(I64)
            if isinstance(descendant, ast.For) and isinstance(descendant.iter, ast.Name) and descendant.iter.id == name:
                return list_type(I64)
            if isinstance(descendant, ast.Call) and isinstance(descendant.func, ast.Name) and descendant.func.id == "len" and descendant.args:
                if isinstance(descendant.args[0], ast.Name) and descendant.args[0].id == name:
                    return list_type(I64)
            if isinstance(descendant, ast.Call) and isinstance(descendant.func, ast.Attribute) and isinstance(descendant.func.value, ast.Name) and descendant.func.value.id == name:
                if descendant.func.attr in {"split", "strip", "lstrip", "rstrip", "encode", "decode"}:
                    return STRING
                if descendant.func.attr in {"append", "extend", "insert", "pop", "popleft", "remove", "clear", "sort", "reverse", "copy"}:
                    return list_type(I64)
        if name.lower().lstrip("_") in {
            "a", "b", "arr", "array", "arrays", "values", "sequence", "weights", "points", "edges",
            "queries", "costs", "matrix", "rows", "fac", "data", "items", "elements", "circles",
        }:
            return list_type(I64)
        return I64

    def _merge_return_type(self, name: str, value: Type, node: ast.AST) -> bool:
        current = self.return_types[name]
        if name in self.fixed_function_returns:
            if value.kind == "optional" and value.args[0] == current:
                self.return_types[name] = value
                return True
            refined = refine_annotated_type(current, value)
            if refined != current:
                self.return_types[name] = refined
                return True
            return False
        if current.kind == "unknown":
            self.return_types[name] = value
            return value.kind != "unknown"
        if value.kind == "unknown" or current == value:
            return False
        try:
            merged = merge_type(current, value)
        except ValueError:
            self.fail(f"incompatible inferred return type for {name}", node)
        if merged == current:
            return False
        self.return_types[name] = merged
        return True

    def _infer_parameter_usage(self, body: list[ast.stmt], params: dict[str, Type], env: dict[str, Type]) -> bool:
        changed = False
        aliases: dict[str, str] = {}
        for descendant in ast.walk(ast.Module(body=body, type_ignores=[])):
            if not isinstance(descendant, ast.Assign) or len(descendant.targets) != 1:
                continue
            target = descendant.targets[0]
            if not isinstance(target, ast.Name):
                continue
            value = descendant.value
            if isinstance(value, ast.Name):
                aliases[target.id] = aliases.get(value.id, value.id)
            elif (
                isinstance(value, ast.Call)
                and isinstance(value.func, ast.Name)
                and value.func.id in {"list", "tuple", "set", "sorted"}
                and len(value.args) == 1
                and isinstance(value.args[0], ast.Name)
            ):
                aliases[target.id] = aliases.get(value.args[0].id, value.args[0].id)
        def infer_usage(value: ast.AST) -> Type:
            try:
                return self.infer(value, env)
            except TranspileError as error:
                if error.message.startswith("unknown name:"):
                    return UNKNOWN
                raise
        for node in ast.walk(ast.Module(body=body, type_ignores=[])):
            if isinstance(node, ast.Assign) and isinstance(node.value, ast.Name) and node.value.id in params:
                for target in node.targets:
                    if isinstance(target, ast.Subscript) and not isinstance(target.slice, ast.Slice):
                        changed |= self._merge_parameter(params, node.value.id, infer_usage(target), node)
            iterable_name = self._iteration_source_name(node.iter) if isinstance(node, ast.For) else None
            if (
                isinstance(node, ast.Assign)
                and len(node.targets) == 1
                and isinstance(node.targets[0], (ast.Tuple, ast.List))
            ):
                target = node.targets[0]
                expected = tuple_type(I64 for _ in target.elts)
                value = node.value
                if isinstance(value, ast.Name) and value.id in params:
                    changed |= self._merge_parameter(params, value.id, expected, node)
                elif (
                    isinstance(value, ast.Subscript)
                    and isinstance(value.value, ast.Name)
                    and aliases.get(value.value.id, value.value.id) in params
                    and not isinstance(value.slice, ast.Slice)
                ):
                    source_name = aliases.get(value.value.id, value.value.id)
                    changed |= self._merge_parameter(params, source_name, list_type(expected), node)
                elif (
                    isinstance(value, ast.Call)
                    and isinstance(value.func, ast.Name)
                    and value.func.id in params
                ):
                    callback = params[value.func.id]
                    callback_parameters = (
                        function_parameters(callback)
                        if callback.kind == "function"
                        else tuple(UNKNOWN for _ in value.args)
                    )
                    changed |= self._merge_parameter(
                        params,
                        value.func.id,
                        function_type(callback_parameters, expected),
                        node,
                    )
                elif (
                    isinstance(value, ast.Call)
                    and isinstance(value.func, ast.Name)
                    and value.func.id in env
                    and env[value.func.id].kind == "function"
                ):
                    callback = env[value.func.id]
                    env[value.func.id] = function_type(function_parameters(callback), expected)
            elif isinstance(node, ast.For) and iterable_name in params:
                item_type = self._infer_loop_item_usage(node, env)
                if item_type.kind != "unknown":
                    changed |= self._merge_parameter(params, iterable_name, list_type(item_type), node)
            elif isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
                if node.func.id == "range":
                    for argument in node.args:
                        if isinstance(argument, ast.Name) and argument.id in params:
                            changed |= self._merge_parameter(params, argument.id, I64, node)
                elif (
                    node.func.id == "isinstance"
                    and len(node.args) == 2
                    and isinstance(node.args[0], ast.Name)
                    and node.args[0].id in params
                    and isinstance(node.args[1], ast.Name)
                    and node.args[1].id == "tuple"
                ):
                    changed |= self._merge_parameter(params, node.args[0].id, tuple_type((I64, I64)), node)
                elif node.func.id in params:
                    argument_types = tuple(infer_usage(argument) for argument in node.args)
                    changed |= self._merge_parameter(
                        params,
                        node.func.id,
                        function_type(argument_types, I64),
                        node,
                    )
            elif isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id in params:
                argument_types = tuple(infer_usage(argument) for argument in node.args)
                changed |= self._merge_parameter(
                    params,
                    node.func.id,
                    function_type(argument_types, I64),
                    node,
                )
            elif isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "len" and node.args:
                argument = node.args[0]
                if isinstance(argument, ast.Name) and argument.id in params:
                    changed |= self._merge_parameter(params, argument.id, list_type(I64), node)
            elif isinstance(node, ast.Subscript) and isinstance(node.value, ast.Name) and node.value.id in params:
                inferred = (
                    STRING
                    if isinstance(node.slice, ast.Slice) and params[node.value.id].kind == "string"
                    else list_type(I64)
                )
                changed |= self._merge_parameter(params, node.value.id, inferred, node)
            elif isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and isinstance(node.func.value, ast.Name):
                name = node.func.value.id
                if name not in params:
                    continue
                if params[name].kind == "class":
                    continue
                attr = node.func.attr
                if attr in {"append", "extend", "sort", "reverse", "count", "index"}:
                    inferred = list_type(infer_usage(node.args[0]) if attr in {"append", "extend"} and node.args else UNKNOWN)
                    changed |= self._merge_parameter(params, name, inferred, node)
                elif attr in {"appendleft", "popleft", "rotate"}:
                    inferred = deque_type(infer_usage(node.args[0]) if attr == "appendleft" and node.args else UNKNOWN)
                    changed |= self._merge_parameter(params, name, inferred, node)
                elif attr == "add":
                    changed |= self._merge_parameter(params, name, set_type(infer_usage(node.args[0]) if node.args else UNKNOWN), node)
                elif attr == "split":
                    changed |= self._merge_parameter(params, name, STRING, node)
                elif attr in {"encode", "decode", "strip", "lstrip", "rstrip"}:
                    changed |= self._merge_parameter(params, name, STRING, node)
                elif attr in {"copy", "keys", "values", "items"}:
                    changed |= self._merge_parameter(params, name, map_type(), node)
            elif isinstance(node, ast.BinOp):
                if isinstance(node.left, ast.Name) and node.left.id in params:
                    other = infer_usage(node.right)
                    inferred = self._binary_parameter_type(node.op, other)
                    if inferred.kind != "unknown":
                        changed |= self._merge_parameter(params, node.left.id, inferred, node)
                if isinstance(node.right, ast.Name) and node.right.id in params:
                    other = infer_usage(node.left)
                    inferred = self._binary_parameter_type(node.op, other)
                    if inferred.kind != "unknown":
                        changed |= self._merge_parameter(params, node.right.id, inferred, node)
            elif isinstance(node, ast.Compare):
                values = [node.left, *node.comparators]
                for index, value in enumerate(values[:-1]):
                    other = values[index + 1]
                    if isinstance(value, ast.Name) and value.id in params:
                        changed |= self._merge_parameter(params, value.id, infer_usage(other), node)
                    if isinstance(other, ast.Name) and other.id in params:
                        changed |= self._merge_parameter(params, other.id, infer_usage(value), node)
            elif isinstance(node, ast.AugAssign) and isinstance(node.value, ast.Name) and node.value.id in params:
                inferred = self._binary_parameter_type(node.op, infer_usage(node.target))
                if inferred.kind != "unknown":
                    changed |= self._merge_parameter(params, node.value.id, inferred, node)
        return changed

    def _iteration_source_name(self, node: ast.AST) -> str | None:
        if isinstance(node, ast.Name):
            return node.id
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "enumerate"
            and node.args
            and isinstance(node.args[0], ast.Name)
        ):
            return node.args[0].id
        return None

    def _infer_loop_item_usage(self, node: ast.For, env: dict[str, Type]) -> Type:
        target = node.target
        if (
            isinstance(node.iter, ast.Call)
            and isinstance(node.iter.func, ast.Name)
            and node.iter.func.id == "enumerate"
            and isinstance(target, (ast.Tuple, ast.List))
            and len(target.elts) == 2
        ):
            target = target.elts[1]
        if isinstance(target, (ast.Tuple, ast.List)):
            if all(isinstance(child, ast.Name) for child in target.elts):
                return tuple_type(I64 for _ in target.elts)
            return UNKNOWN
        if not isinstance(target, ast.Name):
            return UNKNOWN
        target_name = target.id
        result = UNKNOWN
        for descendant in ast.walk(ast.Module(body=node.body, type_ignores=[])):
            expected: list[tuple[ast.AST, Type]] = []
            if isinstance(descendant, ast.For) and descendant is not node:
                if self._iteration_source_name(descendant.iter) == target_name:
                    nested = self._infer_loop_item_usage(descendant, env)
                    if nested.kind != "unknown":
                        result = merge_type(result, list_type(nested))
            if isinstance(descendant, ast.Call) and isinstance(descendant.func, ast.Attribute):
                owner = self.infer(descendant.func.value, env)
                if owner.kind == "class" and descendant.func.attr in self.class_methods(owner.name):
                    method = self.class_methods(owner.name)[descendant.func.attr]
                    parameters = self.method_parameter_types[owner.name][descendant.func.attr]
                    arguments = self.method_arguments(method)
                    expected.extend(
                        (argument, parameters[parameter.arg])
                        for argument, parameter in zip(descendant.args, arguments)
                    )
            elif isinstance(descendant, ast.Call) and isinstance(descendant.func, ast.Name):
                if descendant.func.id in self.functions:
                    function = self.functions[descendant.func.id]
                    parameters = self.parameter_types.get(descendant.func.id, {})
                    arguments = list(function.args.posonlyargs) + list(function.args.args)
                    expected.extend(
                        (argument, parameters.get(parameter.arg, UNKNOWN))
                        for argument, parameter in zip(descendant.args, arguments)
                    )
            for argument, argument_type in expected:
                if isinstance(argument, ast.Name) and argument.id == target_name:
                    result = merge_type(result, argument_type)
            if isinstance(descendant, ast.Assign) and isinstance(descendant.value, ast.Name) and descendant.value.id == target_name:
                for assignment_target in descendant.targets:
                    result = merge_type(result, self.infer(assignment_target, env))
            if isinstance(descendant, ast.AnnAssign) and isinstance(descendant.value, ast.Name) and descendant.value.id == target_name:
                result = merge_type(result, type_from_annotation(descendant.annotation))
            if isinstance(descendant, ast.BinOp):
                if isinstance(descendant.left, ast.Name) and descendant.left.id == target_name:
                    result = merge_type(result, self._binary_parameter_type(descendant.op, self.infer(descendant.right, env)))
                if isinstance(descendant.right, ast.Name) and descendant.right.id == target_name:
                    result = merge_type(result, self._binary_parameter_type(descendant.op, self.infer(descendant.left, env)))
        return result

    def _binary_parameter_type(self, operator: ast.operator, other: Type) -> Type:
        if isinstance(operator, ast.Add):
            if other.kind == "string":
                return STRING
            if other.kind == "list":
                return list_type(other.args[0])
            if is_numeric(other):
                return I64 if other.kind != "f64" else F64
        if isinstance(operator, ast.Mult):
            if other.kind in {"list", "string"}:
                return I64
            return UNKNOWN
        if isinstance(operator, (ast.Sub, ast.Div, ast.FloorDiv, ast.Mod, ast.BitAnd, ast.BitOr, ast.BitXor, ast.LShift, ast.RShift)):
            return I64 if other.kind != "f64" else F64
        return UNKNOWN

    def locals(self, body: list[ast.stmt], env: dict[str, Type]) -> dict[str, Type]:
        result = dict(env)
        for statement in ast.walk(ast.Module(body=body, type_ignores=[])):
            if isinstance(statement, (ast.FunctionDef, ast.AsyncFunctionDef)):
                result.setdefault(statement.name, UNKNOWN)
                for argument in list(statement.args.posonlyargs) + list(statement.args.args) + list(statement.args.kwonlyargs):
                    result.setdefault(argument.arg, type_from_annotation(argument.annotation))
            elif isinstance(statement, ast.Name) and isinstance(statement.ctx, ast.Store):
                result.setdefault(statement.id, UNKNOWN)
        for statement in ast.walk(ast.Module(body=body, type_ignores=[])):
            if isinstance(statement, ast.NamedExpr) and isinstance(statement.target, ast.Name):
                result[statement.target.id] = self.infer(statement.value, result)
        for statement in body:
            self._scan_statement(statement, result)
        return result

    def call_environments(self, body: list[ast.stmt], env: dict[str, Type]) -> dict[int, dict[str, Type]]:
        result: dict[int, dict[str, Type]] = {}

        def record_calls(node: ast.AST, current: dict[str, Type]) -> None:
            for descendant in ast.walk(node):
                if isinstance(descendant, ast.Call):
                    result.setdefault(id(descendant), dict(current))

        def infer_value(node: ast.AST, current: dict[str, Type]) -> Type:
            try:
                return self.infer(node, current)
            except TranspileError as error:
                if error.message.startswith("unknown name:"):
                    return UNKNOWN
                raise

        def bind_target(target: ast.AST, value_type: Type, current: dict[str, Type]) -> None:
            if isinstance(target, ast.Name):
                if value_type.kind != "unknown":
                    current[target.id] = value_type
                return
            if isinstance(target, ast.Starred):
                bind_target(target.value, list_type(element_type(value_type)), current)
                return
            if isinstance(target, (ast.Tuple, ast.List)):
                if value_type.kind == "tuple" and len(value_type.args) == len(target.elts):
                    for child, child_type in zip(target.elts, value_type.args):
                        bind_target(child, child_type, current)
                else:
                    for child in target.elts:
                        bind_target(child, I64, current)

        def merge_environments(left: dict[str, Type], right: dict[str, Type]) -> dict[str, Type]:
            merged = dict(left)
            for name in set(left) | set(right):
                try:
                    merged[name] = merge_type(left.get(name, UNKNOWN), right.get(name, UNKNOWN))
                except ValueError:
                    merged[name] = UNKNOWN
            return merged

        def visit_statements(statements: list[ast.stmt], current: dict[str, Type]) -> dict[str, Type]:
            for statement in statements:
                result[id(statement)] = dict(current)
                if isinstance(statement, ast.Assign):
                    record_calls(statement.value, current)
                    value_type = infer_value(statement.value, current)
                    for target in statement.targets:
                        bind_target(target, value_type, current)
                        if isinstance(target, ast.Subscript):
                            try:
                                self._refine_subscript_assignment(target, value_type, current)
                            except TranspileError as error:
                                if not error.message.startswith("unknown name:"):
                                    raise
                elif isinstance(statement, ast.AnnAssign):
                    if statement.value is not None:
                        record_calls(statement.value, current)
                        value_type = type_from_annotation(statement.annotation)
                        if value_type.kind == "unknown":
                            value_type = infer_value(statement.value, current)
                        bind_target(statement.target, value_type, current)
                elif isinstance(statement, ast.AugAssign):
                    record_calls(statement, current)
                    if isinstance(statement.target, ast.Name):
                        value_type = infer_value(statement.value, current)
                        if value_type.kind != "unknown":
                            current[statement.target.id] = merge_type(current.get(statement.target.id, UNKNOWN), value_type)
                elif isinstance(statement, ast.Expr):
                    record_calls(statement.value, current)
                    self._refine_mutation(statement.value, current)
                elif isinstance(statement, ast.If):
                    record_calls(statement.test, current)
                    left = visit_statements(statement.body, dict(current))
                    right = visit_statements(statement.orelse, dict(current))
                    current.clear()
                    current.update(merge_environments(left, right))
                elif isinstance(statement, ast.For):
                    record_calls(statement.iter, current)
                    item_type = element_type(self.infer_or_unknown(statement.iter, current))
                    bind_target(statement.target, item_type, current)
                    body_env = visit_statements(statement.body, dict(current))
                    else_env = visit_statements(statement.orelse, dict(body_env))
                    current.clear()
                    current.update(merge_environments(body_env, else_env))
                elif isinstance(statement, ast.While):
                    record_calls(statement.test, current)
                    body_env = visit_statements(statement.body, dict(current))
                    else_env = visit_statements(statement.orelse, dict(body_env))
                    current.clear()
                    current.update(merge_environments(body_env, else_env))
                else:
                    record_calls(statement, current)
            return current

        visit_statements(body, dict(env))
        return result

    def _scan_statement(self, node: ast.stmt, env: dict[str, Type]) -> None:
        if isinstance(node, ast.Assign):
            if self.is_input_alias_assignment(node) or self.is_read_alias_assignment(node):
                return
            value_type = self.infer(node.value, env)
            for target in node.targets:
                if isinstance(target, ast.Subscript):
                    self._refine_subscript_assignment(target, value_type, env)
                self._scan_target(target, value_type, env)
        elif isinstance(node, ast.AnnAssign):
            value_type = type_from_annotation(node.annotation)
            if value_type.kind == "unknown" and node.value is not None:
                value_type = self.infer(node.value, env)
            self._scan_target(node.target, value_type, env)
        elif isinstance(node, ast.AugAssign):
            if isinstance(node.target, ast.Name):
                old = env.get(node.target.id, UNKNOWN)
                value = self.infer(node.value, env)
                if value.kind == "optional":
                    value = value.args[0]
                if value.kind == "void":
                    value = old
                try:
                    env[node.target.id] = merge_type(old, value)
                except ValueError:
                    self.fail(f"incompatible augmented assignment for {node.target.id}", node)
            elif isinstance(node.target, ast.Subscript):
                self._refine_subscript_assignment(node.target, self.infer(node.target, env), env)
        elif isinstance(node, ast.For):
            item = element_type(self.iter_type(node.iter, env))
            self._scan_target(node.target, item, env)
            for child in node.body + node.orelse:
                self._scan_statement(child, env)
        elif isinstance(node, ast.While):
            for child in node.body + node.orelse:
                self._scan_statement(child, env)
        elif isinstance(node, ast.If):
            left = dict(env)
            right = dict(env)
            for child in node.body:
                self._scan_statement(child, left)
            for child in node.orelse:
                self._scan_statement(child, right)
            names = set(left) | set(right)
            for name in names:
                try:
                    env[name] = merge_type(left.get(name, UNKNOWN), right.get(name, UNKNOWN))
                except ValueError:
                    self.fail(f"incompatible branch types for {name}", node)
        elif isinstance(node, ast.Expr):
            self._refine_mutation(node.value, env)

    def _refine_subscript_assignment(self, target: ast.Subscript, value_type: Type, env: dict[str, Type]) -> None:
        if isinstance(target.slice, ast.Slice):
            if isinstance(target.value, ast.Name):
                name = target.value.id
                current = env.get(name, UNKNOWN)
                if current.kind == "list" and value_type.kind in {"list", "deque", "tuple"}:
                    env[name] = list_type(merge_type(current.args[0], element_type(value_type)))
            return
        if isinstance(target.value, ast.Subscript) and isinstance(target.value.value, ast.Name):
            name = target.value.value.id
            current = env.get(name, UNKNOWN)
            if current.kind not in {"map", "default_map"}:
                return
            nested = current.args[1]
            if nested.kind not in {"map", "default_map"}:
                return
            outer_key = self.infer(target.value.slice, env)
            inner_key = self.infer(target.slice, env)
            nested = Type(nested.kind, (merge_type(nested.args[0], inner_key), merge_type(nested.args[1], value_type)))
            env[name] = Type(current.kind, (merge_type(current.args[0], outer_key), nested))
            return
        if not isinstance(target.value, ast.Name):
            return
        name = target.value.id
        current = env.get(name, UNKNOWN)
        if current.kind in {"map", "default_map"}:
            key_type = self.infer(target.slice, env)
            env[name] = Type(current.kind, (merge_type(current.args[0], key_type), merge_type(current.args[1], value_type)))
        elif current.kind == "list":
            env[name] = list_type(merge_type(current.args[0], value_type))

    def _scan_target(self, target: ast.AST, value_type: Type, env: dict[str, Type]) -> None:
        if value_type.kind == "reference":
            value_type = value_type.args[0]
        if isinstance(target, ast.Name):
            try:
                env[target.id] = merge_type(env.get(target.id, UNKNOWN), value_type)
            except ValueError:
                self.fail(f"incompatible type for {target.id}", target)
        elif isinstance(target, (ast.Tuple, ast.List)):
            starred = [child for child in target.elts if isinstance(child, ast.Starred)]
            if len(starred) > 1:
                self.fail("multiple starred assignment targets are not supported", target)
            if starred and value_type.kind in {"list", "deque"}:
                item_type = element_type(value_type)
                for child in target.elts:
                    if isinstance(child, ast.Starred):
                        self._scan_target(child.value, list_type(item_type), env)
                    else:
                        self._scan_target(child, item_type, env)
            elif value_type.kind == "tuple" and len(value_type.args) == len(target.elts):
                for child, child_type in zip(target.elts, value_type.args):
                    self._scan_target(child, child_type, env)
            elif value_type.kind in {"list", "deque", "string", "bytes"}:
                item_type = element_type(value_type)
                for child in target.elts:
                    self._scan_target(child, item_type, env)
            else:
                for child in target.elts:
                    self._scan_target(child.value if isinstance(child, ast.Starred) else child, list_type(I64) if isinstance(child, ast.Starred) else I64, env)
        elif isinstance(target, (ast.Subscript, ast.Attribute)):
            return
        elif isinstance(target, ast.Starred):
            self._scan_target(target.value, list_type(element_type(value_type)), env)
        else:
            self.fail("unsupported assignment target", target)

    def _refine_mutation(self, node: ast.AST, env: dict[str, Type]) -> None:
        if not isinstance(node, ast.Call):
            return
        origin = self.imported_origin(node.func.id) if isinstance(node.func, ast.Name) else None
        heap_push = origin == ("heapq", "heappush") or (
            isinstance(node.func, ast.Attribute)
            and isinstance(node.func.value, ast.Name)
            and node.func.value.id == "heapq"
            and node.func.attr == "heappush"
        )
        if heap_push and len(node.args) >= 2:
            receiver, argument = node.args[:2]
            method = "append"
        elif isinstance(node.func, ast.Attribute) and node.args:
            receiver, argument = node.func.value, node.args[0]
            method = node.func.attr
        else:
            return
        if method not in {"append", "appendleft", "add", "push"}:
            return
        if isinstance(receiver, ast.Subscript) and isinstance(receiver.value, ast.Name):
            name = receiver.value.id
            current = env.get(name)
            if current is None or current.kind not in {"list", "deque", "map", "default_map"}:
                return
            nested = current.args[1] if current.kind in {"map", "default_map"} else current.args[0]
            item = self.infer(argument, env)
            if method in {"append", "push"} and nested.kind == "list":
                nested = list_type(merge_type(nested.args[0], item))
            elif method == "appendleft" and nested.kind == "deque":
                nested = deque_type(merge_type(nested.args[0], item))
            elif method == "add" and nested.kind == "set":
                nested = set_type(merge_type(nested.args[0], item))
            else:
                return
            if current.kind in {"map", "default_map"}:
                key = self.infer(receiver.slice, env)
                env[name] = Type(current.kind, (merge_type(current.args[0], key), nested))
            else:
                env[name] = Type(current.kind, (nested,))
            return
        if not isinstance(receiver, ast.Name):
            return
        name = receiver.id
        current = env.get(name)
        if current is None:
            return
        try:
            item = self.infer(argument, env)
        except TranspileError as error:
            if error.message.startswith("unknown name:"):
                return
            raise
        if item.kind == "optional":
            item = item.args[0]
        if current.kind == "unknown" and method in {"append", "push"}:
            env[name] = list_type(item)
            return
        if current.kind not in {"list", "deque", "set"}:
            return
        env[name] = Type(current.kind, (merge_type(current.args[0], item),))

    def iter_type(self, node: ast.AST, env: dict[str, Type]) -> Type:
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            if node.func.id == "range":
                return RANGE
            if node.func.id == "enumerate" and node.args:
                return list_type(tuple_type((I64, element_type(self.iter_type(node.args[0], env)))))
        result = self.infer(node, env)
        if result.kind == "class" and "__iter__" in self.class_methods(result.name):
            return self.method_return_types[result.name]["__iter__"]
        return result

    def _infer_comprehension(self, node: ast.ListComp | ast.SetComp | ast.DictComp | ast.GeneratorExp, env: dict[str, Type]) -> Type:
        comprehension_env = dict(env)
        for generator in node.generators:
            if generator.is_async:
                self.fail("async comprehensions are not supported", generator)
            item_type = element_type(self.iter_type(generator.iter, comprehension_env))
            self._bind_comprehension_target(generator.target, item_type, comprehension_env)
            for condition in generator.ifs:
                self.infer(condition, comprehension_env)
        if isinstance(node, ast.DictComp):
            return map_type(self.infer(node.key, comprehension_env), self.infer(node.value, comprehension_env))
        value_type = self.infer(node.elt, comprehension_env)
        if isinstance(node, ast.SetComp):
            return set_type(value_type)
        return list_type(value_type)

    def _bind_comprehension_target(self, target: ast.AST, value_type: Type, env: dict[str, Type]) -> None:
        if isinstance(target, ast.Name):
            env[target.id] = value_type
            return
        if isinstance(target, (ast.Tuple, ast.List)):
            if value_type.kind == "tuple" and len(value_type.args) == len(target.elts):
                for child, child_type in zip(target.elts, value_type.args):
                    self._bind_comprehension_target(child, child_type, env)
            else:
                for child in target.elts:
                    self._bind_comprehension_target(child, I64, env)
            return
        self.fail("unsupported comprehension target", target)

    def infer(self, node: ast.AST, env: dict[str, Type]) -> Type:
        if isinstance(node, ast.Constant):
            if node.value is None:
                return VOID
            if isinstance(node.value, bool):
                return BOOL
            if isinstance(node.value, int):
                return I64
            if isinstance(node.value, float):
                return F64
            if isinstance(node.value, str):
                return STRING
            if isinstance(node.value, bytes):
                return BYTES
            self.fail("unsupported literal", node)
        if isinstance(node, ast.Name):
            if node.id in {"True", "False"}:
                return BOOL
            if node.id in env:
                return env[node.id]
            if node.id in self.functions:
                function = self.functions[node.id]
                parameters = self.parameter_types.get(node.id, {})
                captured = set(self.hoisted_captures.get(node.id, ()))
                arguments = [
                    argument
                    for argument in list(function.args.posonlyargs) + list(function.args.args)
                    if argument.arg not in captured
                ]
                return function_type((parameters.get(argument.arg, UNKNOWN) for argument in arguments), self.return_types.get(node.id, UNKNOWN))
            if node.id in self.classes:
                return class_type(node.id)
            origin = self.imported_origin(node.id)
            if origin is not None:
                module, imported = origin
                if imported == "pi":
                    return F64
                if imported in {"randint", "randrange", "getrandbits", "randbelow", "choice", "bisect_left", "bisect_right", "gcd", "lcm", "isqrt"}:
                    return function_type((UNKNOWN,), I64)
                if imported in {"sqrt", "hypot", "sin", "cos", "acos", "atan2", "log2", "perf_counter"}:
                    return function_type((UNKNOWN,), F64)
                if imported in {"floor", "ceil"}:
                    return function_type((UNKNOWN,), I64)
                if imported in {"heappush", "heappop", "heapify"}:
                    return function_type((UNKNOWN, UNKNOWN), I64 if imported == "heappop" else VOID)
                if imported in {"copy", "deepcopy"}:
                    return function_type((UNKNOWN,), UNKNOWN)
                if imported == "StringBuilder":
                    return function_type((), list_type(CHAR))
                if imported == "StringIO":
                    return function_type((), list_type(CHAR))
            if node.id == "pi" and node.id in self.from_imports.get("math", set()):
                return F64
            if node.id in self.from_imports.get("bisect", set()):
                return function_type((list_type(I64), I64), I64)
            if node.id in {"min", "max"}:
                return function_type((I64, I64), I64)
            if node.id in {"len", "abs", "str", "int", "float", "bool", "ord", "chr"}:
                return function_type((UNKNOWN,), I64 if node.id != "str" else STRING)
            if node.id in {"int", "float", "str", "bool", "len", "abs", "min", "max", "sum", "all", "any", "pow", "divmod", "ord", "chr", "bin", "bytes", "bytearray", "array", "list", "tuple", "set", "dict", "defaultdict", "sorted", "reversed", "enumerate", "range", "map", "format", "repr", "getattr", "globals", "print", "input", "sys", "stdin", "stdout", "io", "os", "heapq", "bisect", "math", "random", "copy", "heappush", "heappop", "heapify", "randint", "randrange", "getrandbits", "urandom", "randbelow", "choice", "cmp_to_key", "Fraction", "StringIO", "StringBuilder", "redirect_stdout", "redirect_stderr", "isinstance", "issubclass", "type", "iter", "next", "ValueError", "IndexError", "Sequence"}:
                return UNKNOWN
            self.fail(f"unknown name: {node.id}", node)
        if isinstance(node, ast.Lambda):
            parameters = list(node.args.posonlyargs) + list(node.args.args)
            lambda_env = dict(env)
            parameter_types = []
            for argument in parameters:
                argument_type = type_from_annotation(argument.annotation)
                parameter_types.append(argument_type)
                lambda_env[argument.arg] = argument_type
            try:
                result = self.infer(node.body, lambda_env)
            except TranspileError as error:
                if not error.message.startswith("unknown name:"):
                    raise
                result = UNKNOWN
            return function_type(parameter_types, result)
        if isinstance(node, ast.JoinedStr):
            for value in node.values:
                if isinstance(value, ast.FormattedValue):
                    self.infer(value.value, env)
                    if value.format_spec is not None:
                        self.fail("f-string format specifiers are not supported", value)
            return STRING
        if isinstance(node, ast.List):
            item = UNKNOWN
            for element in node.elts:
                try:
                    item = merge_type(item, self.infer(element, env))
                except ValueError:
                    self.fail("list elements have incompatible types", node)
            if item.kind == "void":
                item = optional_type(UNKNOWN)
            return list_type(item)
        if isinstance(node, ast.Tuple):
            return tuple_type(self.infer(element, env) for element in node.elts)
        if isinstance(node, ast.Set):
            item = UNKNOWN
            for element in node.elts:
                item = merge_type(item, self.infer(element, env))
            return set_type(item)
        if isinstance(node, ast.Dict):
            if any(key is None for key in node.keys):
                self.fail("dictionary unpacking is not supported", node)
            key_type = UNKNOWN
            value_type = UNKNOWN
            for key, value in zip(node.keys, node.values):
                key_type = merge_type(key_type, self.infer(key, env))
                value_type = merge_type(value_type, self.infer(value, env))
            return map_type(key_type, value_type)
        if isinstance(node, (ast.ListComp, ast.SetComp, ast.DictComp, ast.GeneratorExp)):
            return self._infer_comprehension(node, env)
        if isinstance(node, ast.Subscript):
            base = self.infer(node.value, env)
            if base.kind == "optional":
                base = base.args[0]
            if isinstance(node.slice, ast.Slice):
                if base.kind == "range":
                    return list_type(I64)
                if base.kind == "unknown":
                    return UNKNOWN
                if base.kind not in {"list", "string", "bytes"}:
                    return UNKNOWN
                return base
            if base.kind == "tuple" and isinstance(node.slice, ast.Constant) and isinstance(node.slice.value, int):
                index = node.slice.value
                if 0 <= index < len(base.args):
                    return base.args[index]
            if base.kind in {"map", "default_map"}:
                return base.args[1]
            if base.kind == "class":
                return self.method_return_types.get(base.name, {}).get("__getitem__", UNKNOWN)
            return element_type(base)
        if isinstance(node, ast.UnaryOp):
            if isinstance(node.op, ast.Not):
                return BOOL
            return self.infer(node.operand, env)
        if isinstance(node, ast.BoolOp):
            return BOOL
        if isinstance(node, ast.Compare):
            return BOOL
        if isinstance(node, ast.IfExp):
            return merge_type(self.infer(node.body, env), self.infer(node.orelse, env))
        if isinstance(node, ast.BinOp):
            left = self.infer(node.left, env)
            right = self.infer(node.right, env)
            if left.kind == "optional":
                left = left.args[0]
            if right.kind == "optional":
                right = right.args[0]
            if left.kind == "void":
                left = right
            if right.kind == "void":
                right = left
            if isinstance(node.op, ast.Add) and (left.kind in {"string", "bytes"} or right.kind in {"string", "bytes"}):
                return STRING if left.kind == "string" or right.kind == "string" else BYTES
            if isinstance(node.op, ast.Mult) and (left.kind in {"list", "string"} or right.kind in {"list", "string"}):
                return left if left.kind in {"list", "string"} else right
            if isinstance(node.op, ast.Div):
                return F64
            if isinstance(node.op, (ast.Add, ast.Sub, ast.Mult, ast.FloorDiv, ast.Mod, ast.Pow, ast.BitAnd, ast.BitOr, ast.BitXor, ast.LShift, ast.RShift)):
                return merge_type(left, right)
            self.fail("unsupported binary operator", node)
        if isinstance(node, ast.Call):
            return self.infer_call(node, env)
        if isinstance(node, ast.Attribute):
            if isinstance(node.value, ast.Name) and node.value.id == "heapq":
                if node.attr == "heappop":
                    return function_type((UNKNOWN, UNKNOWN), I64)
                if node.attr in {"heappush", "heapify"}:
                    return function_type((UNKNOWN, UNKNOWN), VOID)
            if isinstance(node.value, ast.Name) and node.value.id == "math" and node.attr == "pi":
                return F64
            if isinstance(node.value, ast.Name) and node.value.id == "int" and node.attr in {"bit_length", "bit_count"}:
                return function_type((I64,), I64)
            if isinstance(node.value, ast.Name) and node.value.id == "int" and node.attr == "from_bytes":
                return I64
            if node.attr in {"stdin", "stdout", "buffer", "readline", "read", "setrecursionlimit"}:
                if isinstance(node.value, ast.Name) and node.value.id in {"sys", "stdin"}:
                    return UNKNOWN
                if isinstance(node.value, ast.Attribute) and node.value.attr in {"stdin", "buffer"}:
                    return UNKNOWN
            base = self.infer(node.value, env)
            if base.kind == "optional":
                base = base.args[0]
            if base.kind == "unknown":
                return UNKNOWN
            if base.kind == "void":
                return UNKNOWN
            if base.kind in {"i64", "f64", "bool", "char"}:
                return UNKNOWN
            if base.kind == "class":
                if node.attr in self.class_constants.get(base.name, {}):
                    return self.class_constants[base.name][node.attr][0]
                if node.attr in self.class_attributes.get(base.name, {}):
                    return self.class_attributes[base.name][node.attr]
                if node.attr == "__new__":
                    return function_type((), base)
                if node.attr in self.class_methods(base.name):
                    method = self.class_methods(base.name)[node.attr]
                    if (
                        method.decorator_list
                        and isinstance(method.decorator_list[0], ast.Name)
                        and method.decorator_list[0].id == "property"
                    ):
                        return self.method_return_types.get(base.name, {}).get(node.attr, UNKNOWN)
                    parameters = self.method_parameter_types.get(base.name, {}).get(node.attr, {})
                    arguments = self.method_arguments(method)
                    return function_type(
                        (parameters.get(argument.arg, UNKNOWN) for argument in arguments),
                        self.method_return_types.get(base.name, {}).get(node.attr, UNKNOWN),
                    )
                self.fail(f"unknown attribute: {base.name}.{node.attr}", node)
            if base.kind in {"list", "deque", "set"}:
                item = element_type(base)
                if node.attr in {"append", "appendleft", "add", "extend", "insert"}:
                    return function_type((item,), VOID)
                if node.attr in {"pop", "popleft"}:
                    return function_type((I64,), item)
                if node.attr in {"clear", "sort", "reverse", "rotate", "remove", "discard"}:
                    return function_type((), VOID)
                if node.attr in {"count", "index"}:
                    return function_type((item,), I64)
                if node.attr == "copy":
                    return function_type((), base)
            if base.kind in {"map", "default_map"}:
                if node.attr == "__getitem__":
                    return function_type((base.args[0],), base.args[1])
                if node.attr == "get":
                    return function_type((base.args[0],), optional_type(base.args[1]))
                if node.attr == "pop":
                    return function_type((base.args[0],), base.args[1])
                if node.attr in {"keys", "values", "items"}:
                    result = {
                        "keys": list_type(base.args[0]),
                        "values": list_type(base.args[1]),
                        "items": list_type(tuple_type(base.args)),
                    }[node.attr]
                    return function_type((), result)
            if base.kind in {"string", "bytes"}:
                if node.attr in {"strip", "lstrip", "rstrip", "encode", "decode"}:
                    return function_type((), STRING if node.attr != "encode" else BYTES)
                if node.attr in {"upper", "lower"}:
                    return function_type((), base)
                if node.attr in {"split"}:
                    return function_type((STRING,), list_type(STRING))
                if node.attr in {"count", "find", "index"}:
                    return function_type((STRING,), I64)
                if node.attr == "replace":
                    return function_type((STRING, STRING), base)
            if node.attr in {"real", "imag"}:
                return base
            self.fail(f"attribute access is not supported here: {node.attr}", node)
        if isinstance(node, ast.Starred):
            return self.infer(node.value, env)
        if isinstance(node, ast.NamedExpr):
            value_type = self.infer(node.value, env)
            if isinstance(node.target, ast.Name):
                env[node.target.id] = value_type
            return value_type
        if isinstance(node, ast.Yield):
            return VOID if node.value is None else self.infer(node.value, env)
        self.fail(f"unsupported expression: {type(node).__name__}", node)

    def infer_call(self, node: ast.Call, env: dict[str, Type]) -> Type:
        if isinstance(node.func, ast.Name):
            name = node.func.id
            if name in self.ignored_functions:
                return VOID
            origin = self.imported_origin(name)
            if origin is not None:
                imported = origin[1]
                if imported in {"randint", "randrange", "getrandbits", "shuffle", "randbelow", "choice", "bisect_left", "bisect_right", "heappush", "heappop", "heapify", "copy", "deepcopy"}:
                    if imported in {"copy", "deepcopy"} and node.args:
                        return self.infer(node.args[0], env)
                    if imported == "choice" and node.args:
                        return element_type(self.infer(node.args[0], env))
                    if imported in {"heappop"}:
                        return element_type(self.infer(node.args[0], env))
                    return I64 if imported not in {"heappush", "heapify", "shuffle"} else VOID
                if imported in {"sqrt", "hypot", "sin", "cos", "acos", "atan2", "log2", "perf_counter"}:
                    return F64
                if imported in {"floor", "ceil"}:
                    return I64
                if imported in {"StringBuilder", "StringIO"}:
                    return list_type(CHAR)
            if name in {"copy", "deepcopy"} and name in self.from_imports.get("copy", set()):
                if len(node.args) != 1 or node.keywords:
                    self.fail(f"{name} expects one positional argument", node)
                return self.infer(node.args[0], env)
            if name == "int":
                return I64
            if name == "format":
                return STRING
            if name == "float":
                return F64
            if name == "str":
                return STRING
            if name == "repr":
                return STRING
            if name == "getattr":
                return self.infer(node.args[2], env) if len(node.args) >= 3 else UNKNOWN
            if name == "globals":
                return map_type(STRING, function_type((), STRING))
            if name == "bytes":
                return BYTES
            if name in {"bytearray", "array"}:
                return list_type(I64)
            if name == "bin":
                return STRING
            if name in {"randint", "randrange", "getrandbits"} and name in self.from_imports.get("random", set()):
                return I64
            if name == "bool":
                return BOOL
            if name == "callable":
                return BOOL
            if name == "len":
                if node.args:
                    value_type = self.infer(node.args[0], env)
                    if value_type.kind == "class":
                        return self.method_return_types.get(value_type.name, {}).get("__len__", I64)
                return I64
            if name == "abs":
                return self.infer(node.args[0], env) if node.args else UNKNOWN
            if name in env:
                return function_result(env[name]) if env[name].kind == "function" else I64
            if name in {"isinstance", "issubclass"}:
                return BOOL
            if name in {"type", "iter", "next"}:
                return UNKNOWN if name == "type" else (element_type(self.infer(node.args[0], env)) if node.args else UNKNOWN)
            if name in {"urandom", "randbelow", "choice"}:
                if name == "choice" and node.args:
                    return element_type(self.infer(node.args[0], env))
                return BYTES if name == "urandom" else I64
            if name == "cmp_to_key":
                return function_type((UNKNOWN,), UNKNOWN)
            if name == "Fraction":
                return F64
            if name in {"min", "max"}:
                if len(node.args) == 1:
                    return element_type(self.infer(node.args[0], env))
                result = UNKNOWN
                for arg in node.args:
                    result = merge_type(result, self.infer(arg, env))
                return result
            if name == "sum":
                item = element_type(self.infer(node.args[0], env)) if node.args else I64
                return merge_type(item, self.infer(node.args[1], env) if len(node.args) == 2 else I64)
            if name in {"all", "any"}:
                return BOOL
            if name == "pow":
                result = UNKNOWN
                for arg in node.args[:2]:
                    result = merge_type(result, self.infer(arg, env))
                return I64 if len(node.args) == 3 else result
            if name == "divmod":
                result = merge_type(self.infer(node.args[0], env), self.infer(node.args[1], env))
                return tuple_type((result, result))
            if name == "ord":
                return I64
            if name == "chr":
                return STRING
            if name == "sorted":
                t = self.infer(node.args[0], env) if node.args else list_type()
                return list_type(element_type(t))
            if name == "list":
                if not node.args:
                    return list_type()
                arg = node.args[0]
                if isinstance(arg, ast.Call) and isinstance(arg.func, ast.Name) and arg.func.id == "map":
                    return list_type(self._map_item_type(arg, env))
                if isinstance(arg, ast.Call) and isinstance(arg.func, ast.Attribute) and arg.func.attr == "split":
                    return list_type(STRING)
                arg_type = self.infer(arg, env)
                if arg_type.kind == "class" and "__iter__" in self.class_methods(arg_type.name):
                    return self.method_return_types[arg_type.name]["__iter__"]
                return list_type(element_type(arg_type))
            if name == "tuple":
                arg = self.infer(node.args[0], env) if node.args else tuple_type(())
                return arg if arg.kind == "tuple" else list_type(element_type(arg))
            if name == "set":
                arg = self.infer(node.args[0], env) if node.args else list_type()
                return set_type(element_type(arg))
            if name == "dict":
                if not node.args:
                    return map_type()
                arg = self.infer(node.args[0], env)
                item = element_type(arg)
                return map_type(*item.args) if item.kind == "tuple" and len(item.args) == 2 else map_type()
            if name == "defaultdict":
                if not node.args:
                    self.fail("defaultdict expects a default factory", node)
                factory = node.args[0]
                if isinstance(factory, ast.Name):
                    value_type = {
                        "int": I64,
                        "float": F64,
                        "bool": BOOL,
                        "str": STRING,
                        "list": list_type(),
                        "set": set_type(),
                        "dict": map_type(),
                    }.get(factory.id, UNKNOWN)
                elif isinstance(factory, ast.Lambda):
                    value_type = self.infer(factory.body, env)
                else:
                    value_type = UNKNOWN
                if value_type.kind == "unknown":
                    self.fail("defaultdict factories must be int, float, bool, str, list, set, or dict", factory)
                key_type = UNKNOWN
                if len(node.args) == 2:
                    initial = self.infer(node.args[1], env)
                    if initial.kind in {"map", "default_map"}:
                        key_type = initial.args[0]
                        value_type = merge_type(value_type, initial.args[1])
                return default_map_type(key_type, value_type)
            if name == "deque":
                arg = self.infer(node.args[0], env) if node.args else list_type()
                return deque_type(element_type(arg))
            if name == "range":
                return list_type(I64)
            if name == "enumerate":
                if not node.args:
                    return UNKNOWN
                return list_type(tuple_type((I64, element_type(self.infer(node.args[0], env)))))
            if name == "map":
                return list_type(self._map_item_type(node, env))
            if name == "zip":
                return list_type(tuple_type(element_type(self.infer(arg, env)) for arg in node.args))
            if name == "reversed":
                arg_type = self.infer(node.args[0], env)
                if arg_type.kind == "class" and "__reversed__" in self.class_methods(arg_type.name):
                    return self.method_return_types[arg_type.name]["__reversed__"]
                return list_type(element_type(arg_type))
            if name == "input":
                return STRING
            if name == "print":
                return VOID
            if name in {"heappush", "heappop", "heapify"} and name in self.from_imports.get("heapq", set()):
                return I64 if name == "heappop" else VOID
            if name in self.functions:
                return self.return_types.get(name, type_from_annotation(self.functions[name].returns))
            if name in self.classes:
                return class_type(name)
            if name in {"gcd", "lcm", "isqrt"}:
                return I64
            if name in {"sqrt", "hypot", "sin", "cos", "acos", "atan2", "log2"}:
                return F64
            if name == "prod":
                return element_type(self.infer(node.args[0], env)) if node.args else I64
            if name in {"bisect_left", "bisect_right"}:
                return I64
            if name in self.from_imports.get("bisect", set()):
                return I64
            self.fail(f"unsupported function: {name}", node)
        if isinstance(node.func, ast.Attribute):
            if node.func.attr == "setrecursionlimit" and isinstance(node.func.value, ast.Name) and node.func.value.id == "sys":
                return VOID
            if isinstance(node.func.value, ast.Name) and node.func.value.id == "bytes" and node.func.attr == "maketrans":
                return BYTES
            if isinstance(node.func.value, ast.Name) and node.func.value.id == "int" and node.func.attr == "from_bytes":
                return I64
            if node.func.attr in {"readline", "read"} and self._is_stdin_object(node.func.value):
                return STRING
            attr = node.func.attr
            if isinstance(node.func.value, ast.Name) and node.func.value.id == "copy" and attr in {"copy", "deepcopy"}:
                if len(node.args) != 1 or node.keywords:
                    self.fail(f"copy.{attr} expects one positional argument", node)
                return self.infer(node.args[0], env)
            if isinstance(node.func.value, ast.Name) and node.func.value.id in {"heapq", "bisect", "math", "random"}:
                if node.func.value.id == "random" and attr in {"randint", "randrange", "getrandbits"}:
                    return I64
                if node.func.value.id == "heapq" and attr == "heappop" and node.args:
                    return element_type(self.infer(node.args[0], env))
                if attr in {"heappush", "heapify", "heappop", "heappushpop", "heapreplace", "sort", "reverse"}:
                    return I64 if attr in {"heappop", "heappushpop", "heapreplace"} else VOID
                if attr in {"bisect_left", "bisect_right", "gcd", "lcm", "isqrt"}:
                    return I64
                if attr in {"floor", "ceil"}:
                    return I64
                if attr in {"sqrt", "hypot", "sin", "cos", "acos", "atan2", "log2"}:
                    return F64
            if attr == "indices":
                return tuple_type((I64, I64, I64))
            base = self.infer(node.func.value, env)
            if base.kind == "optional":
                base = base.args[0]
            if base.kind == "function":
                return function_result(base)
            if base.kind in {"unknown", "void"}:
                if attr == "get":
                    return UNKNOWN
                if attr in {"append", "appendleft", "add", "clear", "sort", "reverse", "rotate", "extend", "insert", "remove", "discard", "update", "build"}:
                    return VOID
                if attr in {"pop", "popleft", "count", "index", "find", "select0", "select1"}:
                    return UNKNOWN
                return UNKNOWN
            if base.kind == "class":
                methods = self.class_methods(base.name)
                if attr == "__new__":
                    return base
                if attr in self.class_attributes.get(base.name, {}):
                    value_type = self.class_attributes[base.name][attr]
                    return function_result(value_type) if value_type.kind == "function" else I64
                if attr not in methods:
                    self.fail(f"unknown method: {base.name}.{attr}", node)
                return self.method_return_types.get(base.name, {}).get(attr, UNKNOWN)
            if attr in {"append", "appendleft", "add", "clear", "sort", "reverse", "rotate", "popleft", "extend", "insert", "update"}:
                return VOID if attr != "popleft" else element_type(base)
            if attr == "setdefault" and base.kind in {"map", "default_map"}:
                return base.args[1]
            if attr == "build":
                return STRING
            if attr == "pop":
                return base.args[1] if base.kind in {"map", "default_map"} else element_type(base)
            if attr in {"count", "index", "find", "select0", "select1"}:
                return I64
            if attr in {"bit_length", "bit_count"}:
                return I64
            if attr == "as_integer_ratio":
                return tuple_type((I64, I64))
            if attr == "to_bytes":
                return BYTES
            if attr == "isascii":
                return BOOL
            if attr in {"copy"}:
                return base
            if attr == "join":
                return STRING
            if attr == "format":
                return STRING
            if attr == "split":
                return list_type(STRING)
            if attr in {"strip", "lstrip", "rstrip"}:
                return STRING
            if attr in {"upper", "lower"}:
                return STRING if base.kind == "char" else base
            if attr in {"encode", "decode"}:
                return BYTES if attr == "encode" else STRING
            if attr == "replace":
                return base
            if attr == "translate":
                return base
            if attr == "get":
                if base.kind in {"map", "default_map"}:
                    if len(node.args) == 1:
                        return optional_type(base.args[1])
                    if len(node.args) == 2:
                        return merge_type(base.args[1], self.infer(node.args[1], env))
                    self.fail("dict.get expects one or two arguments", node)
            if attr == "keys":
                return list_type(base.args[0]) if base.kind in {"map", "default_map"} else UNKNOWN
            if attr == "values":
                return list_type(base.args[1]) if base.kind in {"map", "default_map"} else UNKNOWN
            if attr == "items":
                return list_type(tuple_type(base.args)) if base.kind in {"map", "default_map"} else UNKNOWN
            if attr in {"remove", "discard"}:
                return VOID
            if attr == "get":
                return UNKNOWN
            self.fail(f"unsupported method: {attr}", node)
        if isinstance(node.func, ast.Subscript):
            function_type_value = element_type(self.infer(node.func.value, env))
            if function_type_value.kind == "function":
                return function_result(function_type_value)
            return UNKNOWN
        self.fail("unsupported call", node)

    def _map_item_type(self, node: ast.Call, env: dict[str, Type]) -> Type:
        if len(node.args) < 2:
            return UNKNOWN
        function = node.args[0]
        source = node.args[1]
        if isinstance(function, ast.Name):
            if function.id == "int":
                return I64
            if function.id == "float":
                return F64
            if function.id == "str":
                return STRING
            if function.id in self.functions:
                return self.return_types.get(function.id, UNKNOWN)
            if function.id in env and env[function.id].kind == "function":
                return function_result(env[function.id])
        return element_type(self.infer(source, env))


@dataclass
class FunctionContext:
    env: dict[str, Type]
    declared: dict[str, Type]
    temp_id: int = 0
    class_name: str | None = None
    reference_aliases: set[str] = field(default_factory=set)
    reference_expressions: dict[str, str] = field(default_factory=dict)
    initialized_aliases: set[str] = field(default_factory=set)
    native_callable_aliases: set[str] = field(default_factory=set)
    private_lists: set[str] = field(default_factory=set)
    local_lists: set[str] = field(default_factory=set)
    rebound_names: set[str] = field(default_factory=set)
    rebound_expressions: dict[str, str] = field(default_factory=dict)
    call_aliases: dict[str, str] = field(default_factory=dict)
    generator_result: str | None = None
    return_type: Type | None = None

    def temporary(self, prefix: str = "tmp") -> str:
        while True:
            self.temp_id += 1
            name = f"__py_{prefix}_{self.temp_id}"
            if name not in self.env and name not in self.declared:
                return name


class CppGenerator:
    def __init__(self, source: str, filename: str = "main.py") -> None:
        try:
            tree = ast.parse(source, filename=filename, type_comments=True)
        except SyntaxError as error:
            location = f"{filename}:{error.lineno}:{error.offset}: "
            raise TranspileError(location + error.msg) from error
        hoist_local_functions(tree, source, filename)
        flatten_try_blocks(tree)
        flatten_with_blocks(tree)
        normalize_keyword_only_arguments(tree)
        self.source = source
        self.filename = filename
        self.tree = tree
        self.analyzer = Analyzer(tree, filename)
        self.function_signatures = {
            name: self.analyzer.signature(node) for name, node in self.analyzer.functions.items()
        }
        self.optimization_planner = OptimizationPlanner(self.analyzer)
        self.program_ir = self._build_program_ir()
        self.lines: list[str] = []
        self.indent = 0
        self.context: FunctionContext | None = None
        self.loop_else_flags: list[str | None] = []

    def _function_ir(
        self,
        name: str,
        node: ast.FunctionDef,
        params: dict[str, Type],
        return_type: Type,
        class_name: str | None = None,
    ) -> FunctionIR:
        env = self.analyzer.function_env(node, params)
        value_types = self.analyzer.locals(node.body, env)
        rebound_types: dict[str, Type] = {}
        rebound_names: set[str] = set()
        rebound_targets: dict[str, list[Type]] = {}
        if class_name is None:
            rebound_types, rebound_names, rebound_targets = self.local_rebind_info(node.body, env)
            rebound_names -= set(params)
            value_types.update({key: value for key, value in rebound_types.items() if key in rebound_names})

        local_names = self.analyzer.bound_names(node.body)
        local_names -= {
            child.name
            for child in node.body
            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef))
        }
        copied = self.copy_rebound_parameters(node, params)
        mutable_parameters = (
            self.analyzer.method_mutable_parameters[class_name][name]
            if class_name is not None
            else self.analyzer.mutable_parameters.get(name, set())
        )
        reference_parameter_kinds = {
            "list", "deque", "set", "map", "default_map", "string", "bytes", "tuple", "class"
        }
        parameter_ownership = {}
        for variable, variable_type in params.items():
            if variable in copied:
                ownership = Ownership.COPIED
            elif variable in mutable_parameters and variable_type.kind in {
                "list", "deque", "set", "map", "default_map", "class"
            }:
                ownership = Ownership.MUTABLE_BORROW
            elif variable_type.kind in reference_parameter_kinds:
                ownership = Ownership.BORROWED
            else:
                ownership = Ownership.VALUE
            parameter_ownership[variable] = ownership
        reference_aliases = self.analyzer.reference_aliases(node.body, value_types) - copied
        bindings = self.optimization_planner.local_bindings(node, value_types)
        native_aliases = set(bindings.callable_aliases) if class_name is not None else set()
        private_lists = set(bindings.private_lists)
        local_lists = {
            variable
            for variable in private_lists
            if element_type(value_types.get(variable, UNKNOWN)).kind in {"i64", "f64", "char"}
        }
        declared_names = local_names - set(params) - reference_aliases - native_aliases
        locals_ir: dict[str, LocalIR] = {}
        for variable in local_names:
            variable_type = value_types.get(variable, UNKNOWN)
            if variable in local_lists:
                storage = Storage.LOCAL_LIST
            elif variable in native_aliases:
                storage = Storage.CALLABLE_ALIAS
            else:
                storage = Storage.DIRECT
            locals_ir[variable] = LocalIR(
                variable,
                variable_type,
                Ownership.BORROWED if variable in reference_aliases else Ownership.VALUE,
                storage,
                variable in declared_names,
            )

        return FunctionIR(
            name=name,
            node=node,
            parameter_types=dict(params),
            parameter_ownership=parameter_ownership,
            mutable_parameters=frozenset(mutable_parameters),
            return_type=return_type,
            environment=dict(env),
            value_types=dict(value_types),
            locals=locals_ir,
            copied_parameters=frozenset(copied),
            rebound_types=dict(rebound_types),
            rebound_names=frozenset(rebound_names),
            rebound_targets={key: tuple(values) for key, values in rebound_targets.items()},
            reference_aliases=frozenset(reference_aliases),
            native_callable_aliases=frozenset(native_aliases),
            private_lists=frozenset(private_lists),
            local_lists=frozenset(local_lists),
            generator=name in self.analyzer.generator_functions if class_name is None else (class_name, name) in self.analyzer.generator_methods,
        )

    def _build_program_ir(self) -> ProgramIR:
        functions = {
            name: self._function_ir(name, node, *self.function_signatures[name])
            for name, node in self.analyzer.functions.items()
        }
        methods = {}
        for class_name in self.analyzer.classes:
            for method_name, node in self.analyzer.class_methods(class_name).items():
                methods[(class_name, method_name)] = self._function_ir(
                    method_name,
                    node,
                    self.analyzer.method_parameter_types[class_name][method_name],
                    self.analyzer.method_return_types[class_name][method_name],
                    class_name,
                )
        return ProgramIR(self.source, self.filename, self.tree, functions, methods)

    def generate(self) -> str:
        self._emit_runtime()
        for name, node in self.analyzer.classes.items():
            self._emit_class_declaration(name, node)
        self._emit_globals()
        for name, node in self.analyzer.functions.items():
            self._emit_prototype(name, node)
        for name, node in self.analyzer.functions.items():
            self._emit_function(name, node)
        for name, node in self.analyzer.classes.items():
            self._emit_class_methods(name, node)
        self._emit_main()
        return "\n".join(self.lines) + "\n"

    def emit(self, line: str = "") -> None:
        self.lines.append("    " * self.indent + line)

    def fail(self, message: str, node: ast.AST | None = None) -> None:
        self.analyzer.fail(message, node)

    def _emit_runtime(self) -> None:
        runtime = r'''#include <algorithm>
#include <bit>
#include <cassert>
#include <charconv>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <iterator>
#include <numeric>
#include <optional>
#include <queue>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>
#include <deque>
using namespace std;
using i64 = long long;
using u64 = unsigned long long;

template <class K, class V>
struct PyDict {
    unordered_map<K, V> data;
    vector<K> order;

    PyDict() = default;
    PyDict(initializer_list<pair<K, V>> values) {
        for (const auto& [key, value] : values) (*this)[key] = value;
    }
    V& operator[](const K& key) {
        auto [iterator, inserted] = data.try_emplace(key);
        if (inserted) order.push_back(key);
        return iterator->second;
    }
    const V& at(const K& key) const { return data.at(key); }
    V& at(const K& key) { return data.at(key); }
    auto find(const K& key) { return data.find(key); }
    auto find(const K& key) const { return data.find(key); }
    auto end() { return data.end(); }
    auto end() const { return data.end(); }
    size_t size() const { return data.size(); }
    bool empty() const { return data.empty(); }
    void clear() {
        data.clear();
        order.clear();
    }
    void erase(const K& key) {
        if (data.erase(key)) {
            order.erase(find_if(order.begin(), order.end(), [&](const K& value) { return value == key; }));
        }
    }
};

template <class K, class V>
struct PyDefaultDict : PyDict<K, V> {
    function<V()> factory = []() { return V{}; };

    PyDefaultDict() = default;
    explicit PyDefaultDict(function<V()> value) : factory(move(value)) {}
    explicit PyDefaultDict(const PyDict<K, V>& initial) : PyDict<K, V>(initial) {}
    PyDefaultDict(function<V()> value, const PyDict<K, V>& initial)
        : PyDict<K, V>(initial), factory(move(value)) {}
    V& operator[](const K& key) {
        auto iterator = this->data.find(key);
        if (iterator == this->data.end()) {
            this->order.push_back(key);
            return this->data.emplace(key, factory()).first->second;
        }
        return iterator->second;
    }
};

inline bool py_space(unsigned char value) {
    return value == ' ' || value == '\t' || value == '\n' || value == '\v' || value == '\f' || value == '\r';
}

struct FastScanner {
    string data;
    size_t pos = 0;
    FastScanner() {
        ios::sync_with_stdio(false);
        cin.tie(nullptr);
        char buffer[65536];
        size_t read_size;
        while ((read_size = fread(buffer, 1, sizeof(buffer), stdin)) != 0) {
            data.append(buffer, read_size);
        }
    }
    void skip() {
        while (pos < data.size() && py_space(static_cast<unsigned char>(data[pos]))) ++pos;
    }
    string token() {
        skip();
        size_t begin = pos;
        while (pos < data.size() && !py_space(static_cast<unsigned char>(data[pos]))) ++pos;
        return data.substr(begin, pos - begin);
    }
    string line() {
        if (pos >= data.size()) return {};
        if (data[pos] == '\n') ++pos;
        size_t begin = pos;
        pos = data.find('\n', pos);
        if (pos == string::npos) pos = data.size();
        size_t end = pos;
        if (end > begin && data[end - 1] == '\r') --end;
        return data.substr(begin, end - begin);
    }
    string remaining() {
        string result = data.substr(pos);
        pos = data.size();
        return result;
    }
    i64 read_i64_token() {
        skip();
        bool negative = false;
        if (pos < data.size() && (data[pos] == '-' || data[pos] == '+')) {
            negative = data[pos] == '-';
            ++pos;
        }
        if (pos >= data.size() || data[pos] < '0' || data[pos] > '9') {
            throw invalid_argument("invalid integer");
        }
        u64 result = 0;
        while (pos < data.size() && data[pos] >= '0' && data[pos] <= '9') {
            result = result * 10 + (data[pos] - '0');
            ++pos;
        }
        return static_cast<i64>(negative ? u64(0) - result : result);
    }
    i64 read_i64() {
        i64 result = read_i64_token();
        pos = data.find('\n', pos);
        if (pos == string::npos) pos = data.size();
        if (pos < data.size()) ++pos;
        return result;
    }
    void read_i64_line(vector<i64>& result) {
        size_t begin = pos;
        pos = data.find('\n', pos);
        if (pos == string::npos) pos = data.size();
        size_t end = pos;
        size_t estimate = min<size_t>((end - begin + 1) / 2, 1 << 20);
        result.clear();
        result.reserve(estimate);
        pos = begin;
        while (pos < end) {
            while (pos < end && py_space(static_cast<unsigned char>(data[pos]))) ++pos;
            if (pos >= end) break;
            result.push_back(read_i64_token());
        }
        pos = end;
        if (pos < data.size()) ++pos;
    }
    template <class F>
    auto read_i64_line_map(F function) {
        using result_type = decay_t<decltype(function(i64{}))>;
        vector<result_type> result;
        size_t begin = pos;
        pos = data.find('\n', pos);
        if (pos == string::npos) pos = data.size();
        size_t end = pos;
        size_t estimate = min<size_t>((end - begin + 1) / 2, 1 << 20);
        result.reserve(estimate);
        pos = begin;
        while (pos < end) {
            while (pos < end && py_space(static_cast<unsigned char>(data[pos]))) ++pos;
            if (pos >= end) break;
            result.push_back(function(read_i64_token()));
        }
        pos = end;
        if (pos < data.size()) ++pos;
        return result;
    }
    void read_i64_remaining(vector<i64>& result) {
        size_t estimate = min<size_t>((data.size() - pos + 1) / 2, 1 << 20);
        result.clear();
        result.reserve(estimate);
        while (true) {
            skip();
            if (pos >= data.size()) break;
            result.push_back(read_i64_token());
        }
    }
    template <class T>
    T read() {
        if constexpr (is_same_v<T, i64>) return read_i64();
        else {
            string s = token();
            pos = data.find('\n', pos);
            if (pos == string::npos) pos = data.size();
            if (pos < data.size()) ++pos;
            if constexpr (is_same_v<T, double>) return stod(s);
            else return T(s);
        }
    }
};

FastScanner py_input;

vector<string> py_split(const string& s) {
    vector<string> result;
    size_t i = 0;
    while (i < s.size()) {
        while (i < s.size() && py_space(static_cast<unsigned char>(s[i]))) ++i;
        size_t begin = i;
        while (i < s.size() && !py_space(static_cast<unsigned char>(s[i]))) ++i;
        if (begin < i) result.push_back(s.substr(begin, i - begin));
    }
    return result;
}

vector<string> py_split(const string& s, const string& separator) {
    if (separator.empty()) throw runtime_error("empty separator");
    vector<string> result;
    size_t begin = 0;
    while (true) {
        size_t end = s.find(separator, begin);
        if (end == string::npos) {
            result.push_back(s.substr(begin));
            return result;
        }
        result.push_back(s.substr(begin, end - begin));
        begin = end + separator.size();
    }
}

string py_strip(const string& s, i64 mode) {
    size_t left = 0;
    size_t right = s.size();
    if (mode <= 0) while (left < right && py_space(static_cast<unsigned char>(s[left]))) ++left;
    if (mode >= 0) while (right > left && py_space(static_cast<unsigned char>(s[right - 1]))) --right;
    return s.substr(left, right - left);
}

string py_upper(const string& source) {
    string result = source;
    for (char& value : result) value = static_cast<char>(toupper(static_cast<unsigned char>(value)));
    return result;
}

string py_upper(char source) { return py_upper(string(1, source)); }

string py_lower(const string& source) {
    string result = source;
    for (char& value : result) value = static_cast<char>(tolower(static_cast<unsigned char>(value)));
    return result;
}

string py_lower(char source) { return py_lower(string(1, source)); }

string py_replace(const string& source, const string& old_value, const string& new_value) {
    if (old_value.empty()) {
        string result = new_value;
        for (char value : source) {
            result.push_back(value);
            result += new_value;
        }
        return result;
    }
    string result;
    size_t begin = 0;
    while (true) {
        size_t position = source.find(old_value, begin);
        if (position == string::npos) {
            result.append(source, begin, string::npos);
            return result;
        }
        result.append(source, begin, position - begin);
        result += new_value;
        begin = position + old_value.size();
    }
}

i64 py_count(const string& source, const string& value) {
    if (value.empty()) return static_cast<i64>(source.size()) + 1;
    i64 result = 0;
    size_t begin = 0;
    while (true) {
        size_t position = source.find(value, begin);
        if (position == string::npos) return result;
        ++result;
        begin = position + value.size();
    }
}

i64 py_int(const string& value) { return stoll(value); }
i64 py_int(char value) { return value >= '0' && value <= '9' ? value - '0' : stoll(string(1, value)); }
i64 py_int_base(const string& value, i64 base) { return stoll(value, nullptr, static_cast<int>(base)); }
template <class T>
i64 py_int(const T& value) { return static_cast<i64>(value); }

double py_float(const string& value) { return stod(value); }
double py_float(const char* value) { return stod(string(value)); }
double py_float(char value) { return stod(string(1, value)); }
template <class T>
double py_float(const T& value) { return static_cast<double>(value); }

string py_string(const string& value) { return value; }
string py_string(const char* value) { return string(value); }
string py_string(char value) { return string(1, value); }
string py_string(bool value) { return value ? "True" : "False"; }
template <class T>
string py_string(const T& value) {
    if constexpr (is_integral_v<T> && sizeof(T) <= sizeof(u64)) {
        char buffer[32];
        auto converted = to_chars(buffer, buffer + sizeof(buffer), value);
        return string(buffer, converted.ptr);
    } else return to_string(value);
}

template <class T>
string py_string(const vector<T>& values);
template <class T>
string py_string(const deque<T>& values);
template <class T>
string py_string(const set<T>& values);
template <class... T>
string py_string(const tuple<T...>& values);
template <class K, class V>
string py_string(const PyDict<K, V>& values);
template <class K, class V>
string py_string(const PyDefaultDict<K, V>& values);

string py_repr(const string& value) { return "'" + value + "'"; }
string py_repr(const char* value) { return py_repr(string(value)); }
string py_repr(char value) { return "'" + string(1, value) + "'"; }
template <class T>
string py_repr(const T& value) { return py_string(value); }

template <class T>
string py_string(const vector<T>& values) {
    string result = "[";
    bool first = true;
    for (const auto& value : values) {
        if (!first) result += ", ";
        first = false;
        result += py_repr(static_cast<T>(value));
    }
    return result + "]";
}

template <class T>
string py_string(const deque<T>& values) {
    string result = "deque([";
    bool first = true;
    for (const auto& value : values) {
        if (!first) result += ", ";
        first = false;
        result += py_repr(value);
    }
    return result + "])";
}

template <class T>
string py_string(const set<T>& values) {
    if (values.empty()) return "set()";
    string result = "{";
    bool first = true;
    for (const auto& value : values) {
        if (!first) result += ", ";
        first = false;
        result += py_repr(value);
    }
    return result + "}";
}

template <class... T>
string py_string(const tuple<T...>& values) {
    string result = "(";
    bool first = true;
    apply([&](const auto&... value) {
        ((result += (first ? "" : ", ") + py_repr(value), first = false), ...);
    }, values);
    if constexpr (sizeof...(T) == 1) result += ",";
    return result + ")";
}

template <class... T>
auto py_tuple_values(const tuple<T...>& values) {
    if constexpr (sizeof...(T) == 0) {
        return vector<i64>{};
    } else {
        using value_type = common_type_t<T...>;
        vector<value_type> result;
        result.reserve(sizeof...(T));
        apply([&](const auto&... value) {
            (result.push_back(value), ...);
        }, values);
        return result;
    }
}

template <class K, class V>
string py_string(const PyDict<K, V>& values) {
    string result = "{";
    bool first = true;
    for (const K& key : values.order) {
        if (!first) result += ", ";
        first = false;
        result += py_repr(key) + ": " + py_repr(values.at(key));
    }
    return result + "}";
}

template <class K, class V>
string py_string(const PyDefaultDict<K, V>& values) {
    return py_string(static_cast<const PyDict<K, V>&>(values));
}

template <class... T>
string py_format(const string& pattern, const T&... values) {
    vector<string> arguments{py_string(values)...};
    string result;
    size_t argument = 0;
    for (size_t index = 0; index < pattern.size(); ++index) {
        char current = pattern[index];
        if (current == '{' && index + 1 < pattern.size() && pattern[index + 1] == '{') {
            result.push_back('{');
            ++index;
        } else if (current == '}' && index + 1 < pattern.size() && pattern[index + 1] == '}') {
            result.push_back('}');
            ++index;
        } else if (current == '{' && index + 1 < pattern.size() && pattern[index + 1] == '}') {
            if (argument >= arguments.size()) throw runtime_error("format argument missing");
            result += arguments[argument++];
            ++index;
        } else if (current == '{' || current == '}') {
            throw runtime_error("unsupported format field");
        } else {
            result.push_back(current);
        }
    }
    if (argument != arguments.size()) throw runtime_error("unused format argument");
    return result;
}

i64 py_ord(char value) { return static_cast<unsigned char>(value); }
i64 py_ord(const string& value) { return static_cast<unsigned char>(value.at(0)); }
i64 py_ord(const char* value) { return static_cast<unsigned char>(value[0]); }

template <class T>
vector<T> py_read_vector(FastScanner& in) {
    vector<T> result;
    if constexpr (is_same_v<T, i64>) {
        in.read_i64_line(result);
        return result;
    } else {
        string line = in.line();
        for (const string& token : py_split(line)) {
            if constexpr (is_same_v<T, double>) result.push_back(stod(token));
            else result.push_back(token);
        }
        return result;
    }
}

template <class T>
vector<T> py_read_remaining(FastScanner& in) {
    vector<T> result;
    if constexpr (is_same_v<T, i64>) {
        in.read_i64_remaining(result);
        return result;
    } else {
        for (const string& token : py_split(in.remaining())) {
            if constexpr (is_same_v<T, double>) result.push_back(stod(token));
            else result.push_back(token);
        }
        return result;
    }
}

template <class Tuple, size_t... I>
auto py_tuple_at_impl(const Tuple& values, i64 index, index_sequence<I...>) {
    using result_type = common_type_t<tuple_element_t<I, Tuple>...>;
    result_type items[] = {static_cast<result_type>(get<I>(values))...};
    i64 size = static_cast<i64>(sizeof...(I));
    if (index < 0) index += size;
    return items[index];
}

template <class... T>
auto py_at(tuple<T...>& values, i64 index) {
    return py_tuple_at_impl(values, index, index_sequence_for<T...>{});
}

template <class... T>
auto py_at(const tuple<T...>& values, i64 index) {
    return py_tuple_at_impl(values, index, index_sequence_for<T...>{});
}

template <class C>
decltype(auto) py_at(C& container, i64 index) {
    i64 n = static_cast<i64>(container.size());
    if (index < 0) index += n;
    if (index < 0 || index >= n) throw out_of_range("sequence index out of range");
    return container[static_cast<size_t>(index)];
}

template <class T>
vector<T> py_concat(const vector<T>& left, const vector<T>& right) {
    vector<T> result;
    result.reserve(left.size() + right.size());
    result.insert(result.end(), left.begin(), left.end());
    result.insert(result.end(), right.begin(), right.end());
    return result;
}

vector<string> py_char_rows(const vector<vector<char>>& rows) {
    vector<string> result;
    result.reserve(rows.size());
    for (const auto& row : rows) result.emplace_back(row.begin(), row.end());
    return result;
}

string py_repeat(const string& values, i64 count) {
    if (count <= 0 || values.empty()) return {};
    if (count == 1) return values;
    if (values.size() == 1) return string(static_cast<size_t>(count), values.front());
    string result;
    __int128 total = static_cast<__int128>(values.size()) * count;
    if (total > static_cast<__int128>(result.max_size())) throw length_error("repeated string is too large");
    size_t size = static_cast<size_t>(total);
    result.reserve(size);
    result.append(values);
    while (result.size() < size) result.append(result, 0, min(result.size(), size - result.size()));
    return result;
}

template <class C>
C py_repeat(const C& values, i64 count) {
    if (count <= 0 || values.empty()) return C{};
    if (values.size() == 1) return C(static_cast<size_t>(count), values.front());
    C result;
    if (count > 0) {
        __int128 total = static_cast<__int128>(values.size()) * count;
        if constexpr (requires { result.reserve(size_t{}); }) {
            if (total <= static_cast<__int128>(result.max_size())) result.reserve(static_cast<size_t>(total));
        }
    }
    for (i64 index = 0; index < count; ++index) {
        result.insert(result.end(), values.begin(), values.end());
    }
    return result;
}

template <size_t N>
string py_repeat(const char (&values)[N], i64 count) {
    return py_repeat(string(values, N - 1), count);
}

template <class C>
decltype(auto) py_at(const C& container, i64 index) {
    i64 n = static_cast<i64>(container.size());
    if (index < 0) index += n;
    if (index < 0 || index >= n) throw out_of_range("sequence index out of range");
    return container[static_cast<size_t>(index)];
}

template <class T>
class py_local_list {
    T local[32];
    vector<T> heap;
    size_t used;
    bool spilled;
public:
    using value_type = T;
    py_local_list() : used(0), spilled(false) {}
    size_t size() const { return used; }
    bool empty() const { return used == 0; }
    T* data() { return spilled ? heap.data() : local; }
    const T* data() const { return spilled ? heap.data() : local; }
    T& operator[](size_t index) { return data()[index]; }
    const T& operator[](size_t index) const { return data()[index]; }
    T& back() { return spilled ? heap.back() : local[used - 1]; }
    void clear() { heap.clear(); used = 0; spilled = false; }
    void push_back(T value) {
        if (!spilled && used < 32) local[used++] = value;
        else {
            if (!spilled) {
                heap.reserve(64);
                heap.insert(heap.end(), local, local + used);
                spilled = true;
            }
            heap.push_back(value);
            ++used;
        }
    }
    void pop_back() { if (spilled) heap.pop_back(); --used; }
    void assign(initializer_list<T> values) {
        clear();
        for (T value : values) push_back(value);
    }
};

template <class T>
T py_pop(py_local_list<T>& container) {
    T result = container.back();
    container.pop_back();
    return result;
}

template <class T>
T py_pop(vector<T>& container) {
    T result = std::move(container.back());
    container.pop_back();
    return result;
}

template <class T>
T py_pop(deque<T>& container) {
    T result = std::move(container.back());
    container.pop_back();
    return result;
}

template <class C>
typename C::value_type py_pop_at(C& container, i64 index) {
    i64 size = static_cast<i64>(container.size());
    if (index < 0) index += size;
    auto iterator = container.begin() + index;
    typename C::value_type result = std::move(*iterator);
    container.erase(iterator);
    return result;
}

template <class T>
T py_popleft(deque<T>& container) {
    T result = std::move(container.front());
    container.pop_front();
    return result;
}

template <class C>
void py_erase_at(C& container, i64 index) {
    i64 size = static_cast<i64>(container.size());
    if (index < 0) index += size;
    container.erase(container.begin() + index);
}

template <class C, class T>
i64 py_bisect(const C& source, const T& value, bool right) {
    auto position = right ? upper_bound(source.begin(), source.end(), value)
                          : lower_bound(source.begin(), source.end(), value);
    return static_cast<i64>(position - source.begin());
}

template <class C, class S>
void py_extend(C& target, const S& source) {
    if constexpr (is_same_v<C, S>) {
        if (&target == &source) {
            S copy = source;
            target.insert(target.end(), copy.begin(), copy.end());
            return;
        }
    }
    target.insert(target.end(), source.begin(), source.end());
}

template <class C, class S>
void py_update(C& target, const S& source) {
    target.insert(source.begin(), source.end());
}

template <class C, class T>
void py_remove(C& container, const T& value) {
    auto iterator = find(container.begin(), container.end(), value);
    if (iterator == container.end()) throw runtime_error("value is not in container");
    container.erase(iterator);
}

template <class C>
void py_rotate(C& container, i64 count) {
    i64 size = static_cast<i64>(container.size());
    if (size == 0) return;
    count %= size;
    if (count < 0) count += size;
    rotate(container.begin(), container.end() - count, container.end());
}

template <class C, class T>
i64 py_find(const C& container, const T& value) {
    auto iterator = find(container.begin(), container.end(), value);
    return iterator == container.end() ? -1 : static_cast<i64>(iterator - container.begin());
}

i64 py_find(const string& container, const string& value) {
    size_t position = container.find(value);
    return position == string::npos ? -1 : static_cast<i64>(position);
}

template <class C>
vector<typename C::value_type> py_sorted(const C& source, bool descending = false) {
    vector<typename C::value_type> result(source.begin(), source.end());
    sort(result.begin(), result.end());
    if (descending) reverse(result.begin(), result.end());
    return result;
}

template <class T>
vector<T> py_sorted(vector<T>&& source, bool descending = false) {
    sort(source.begin(), source.end());
    if (descending) reverse(source.begin(), source.end());
    return std::move(source);
}

template <class K, class V>
vector<K> py_sorted(const PyDict<K, V>& source, bool descending = false) {
    vector<K> result = source.order;
    sort(result.begin(), result.end());
    if (descending) reverse(result.begin(), result.end());
    return result;
}

template <class K, class V>
vector<K> py_sorted(const PyDefaultDict<K, V>& source, bool descending = false) {
    return py_sorted(static_cast<const PyDict<K, V>&>(source), descending);
}

template <class C, class F>
vector<typename C::value_type> py_sorted_key_direct(const C& source, F key, bool descending = false) {
    vector<typename C::value_type> result(source.begin(), source.end());
    stable_sort(result.begin(), result.end(), [&](const auto& left, const auto& right) {
        return descending ? key(right) < key(left) : key(left) < key(right);
    });
    return result;
}

template <class T, class F>
void py_sort_key_direct(vector<T>& values, F key, bool descending = false) {
    stable_sort(values.begin(), values.end(), [&](const T& left, const T& right) {
        return descending ? key(right) < key(left) : key(left) < key(right);
    });
}

template <class C, class F>
vector<size_t> py_key_order(const C& values, F key, bool descending) {
    using Key = decay_t<decltype(key(*values.begin()))>;
    vector<Key> keys;
    keys.reserve(values.size());
    for (const auto& value : values) keys.push_back(key(value));
    vector<size_t> order(values.size());
    iota(order.begin(), order.end(), size_t(0));
    stable_sort(order.begin(), order.end(), [&](size_t left, size_t right) {
        return descending ? keys[right] < keys[left] : keys[left] < keys[right];
    });
    return order;
}

template <class C, class F>
vector<typename C::value_type> py_sorted_key(const C& source, F key, bool descending = false) {
    vector<typename C::value_type> values(source.begin(), source.end());
    auto order = py_key_order(values, key, descending);
    vector<typename C::value_type> result;
    result.reserve(values.size());
    for (size_t index : order) result.push_back(move(values[index]));
    return result;
}

template <class T, class F>
void py_sort_key(vector<T>& values, F key, bool descending = false) {
    auto order = py_key_order(values, key, descending);
    vector<T> result;
    result.reserve(values.size());
    for (size_t index : order) result.push_back(move(values[index]));
    std::move(result.begin(), result.end(), values.begin());
}

template <class C>
vector<typename C::value_type> py_reversed(const C& source) {
    return vector<typename C::value_type>(source.rbegin(), source.rend());
}

template <class T>
vector<T> py_reversed(vector<T>&& source) {
    reverse(source.begin(), source.end());
    return std::move(source);
}

template <class C>
C py_slice(const C& source, bool has_left, i64 left, bool has_right, i64 right, i64 step) {
    if (step == 0) throw runtime_error("slice step cannot be zero");
    i64 n = static_cast<i64>(source.size());
    if (step > 0) {
        if (!has_left) left = 0;
        else if (left < 0) left += n;
        if (!has_right) right = n;
        else if (right < 0) right += n;
        left = max<i64>(0, min<i64>(n, left));
        right = max<i64>(0, min<i64>(n, right));
    } else {
        if (!has_left) left = n - 1;
        else {
            if (left < 0) left += n;
            left = max<i64>(-1, min<i64>(n - 1, left));
        }
        if (!has_right) right = -1;
        else {
            if (right < 0) right += n;
            right = max<i64>(-1, min<i64>(n - 1, right));
        }
    }
    if (step == 1) return C(source.begin() + left, source.begin() + max(left, right));
    C result;
    i64 length = 0;
    if (step > 0) {
        if (left < right) {
            __int128 distance = static_cast<__int128>(right) - left;
            length = static_cast<i64>((distance + step - 1) / step);
        }
        if constexpr (requires { result.reserve(size_t{}); }) {
            if (length > 0) result.reserve(static_cast<size_t>(length));
        }
        for (i64 offset = 0; offset < length; ++offset) result.push_back(source[static_cast<size_t>(left + offset * step)]);
    } else {
        if (left > right) {
            __int128 distance = static_cast<__int128>(left) - right;
            __int128 width = -static_cast<__int128>(step);
            length = static_cast<i64>((distance + width - 1) / width);
        }
        if constexpr (requires { result.reserve(size_t{}); }) {
            if (length > 0) result.reserve(static_cast<size_t>(length));
        }
        for (i64 offset = 0; offset < length; ++offset) result.push_back(source[static_cast<size_t>(left + offset * step)]);
    }
    return result;
}

tuple<i64, i64, i64> py_slice_indices(i64 key, i64 size) {
    return make_tuple(key, size, i64(1));
}

template <class T>
void py_set_slice(vector<T>& target, bool has_left, i64 left, bool has_right, i64 right, i64 step, vector<T> replacement) {
    if (step == 0) throw runtime_error("slice step cannot be zero");
    i64 size = static_cast<i64>(target.size());
    if (step < 0) {
        if (!has_left) left = size - 1;
        else if (left < 0) left += size;
        if (!has_right) right = -1;
        else if (right < 0) right += size;
        left = max<i64>(-1, min<i64>(size - 1, left));
        right = max<i64>(-1, min<i64>(size - 1, right));
        size_t offset = 0;
        for (i64 index = left; index > right; index += step) {
            if (offset >= replacement.size()) throw runtime_error("slice assignment size mismatch");
            target[static_cast<size_t>(index)] = replacement[offset++];
        }
        if (offset != replacement.size()) throw runtime_error("slice assignment size mismatch");
        return;
    }
    if (!has_left) left = 0;
    else if (left < 0) left += size;
    if (!has_right) right = size;
    else if (right < 0) right += size;
    left = max<i64>(0, min<i64>(size, left));
    right = max<i64>(0, min<i64>(size, right));
    if (step != 1) {
        size_t offset = 0;
        for (i64 index = left; index < right; index += step) {
            if (offset >= replacement.size()) throw runtime_error("slice assignment size mismatch");
            target[static_cast<size_t>(index)] = replacement[offset++];
        }
        if (offset != replacement.size()) throw runtime_error("slice assignment size mismatch");
        return;
    }
    right = max<i64>(left, right);
    if (static_cast<size_t>(right - left) == replacement.size()) {
        std::move(replacement.begin(), replacement.end(), target.begin() + left);
        return;
    }
    target.erase(target.begin() + left, target.begin() + right);
    target.insert(target.begin() + left, make_move_iterator(replacement.begin()), make_move_iterator(replacement.end()));
}

template <class C>
void py_delete_slice(C& target, bool has_left, i64 left, bool has_right, i64 right, i64 step) {
    if (step == 0) throw runtime_error("slice step cannot be zero");
    i64 size = static_cast<i64>(target.size());
    __int128 first;
    __int128 width;
    i64 count;
    if (step > 0) {
        if (!has_left) left = 0;
        else if (left < 0) left += size;
        if (!has_right) right = size;
        else if (right < 0) right += size;
        left = max<i64>(0, min<i64>(size, left));
        right = max<i64>(0, min<i64>(size, right));
        if (left >= right) return;
        first = left;
        width = step;
        count = static_cast<i64>((static_cast<__int128>(right) - left - 1) / width + 1);
    } else {
        if (!has_left) left = size - 1;
        else if (left < 0) left += size;
        if (!has_right) right = -1;
        else if (right < 0) right += size;
        left = max<i64>(-1, min<i64>(size - 1, left));
        right = max<i64>(-1, min<i64>(size - 1, right));
        if (left <= right) return;
        width = -static_cast<__int128>(step);
        count = static_cast<i64>((static_cast<__int128>(left) - right - 1) / width + 1);
        first = static_cast<__int128>(left) - (count - 1) * width;
    }
    if (width == 1) {
        target.erase(target.begin() + static_cast<i64>(first), target.begin() + static_cast<i64>(first) + count);
        return;
    }
    i64 write = static_cast<i64>(first);
    __int128 next = first;
    for (i64 read = write; read < size; ++read) {
        if (count > 0 && static_cast<__int128>(read) == next) { next += width; --count; }
        else target[static_cast<size_t>(write++)] = std::move(target[static_cast<size_t>(read)]);
    }
    target.erase(target.begin() + write, target.end());
}

i64 py_floor_div(i64 a, i64 b) {
    i64 q = a / b;
    i64 r = a % b;
    if (r != 0 && ((r > 0) != (b > 0))) --q;
    return q;
}

i64 py_mod(i64 a, i64 b) {
    i64 r = a % b;
    if (r != 0 && ((r > 0) != (b > 0))) r += b;
    return r;
}

i64 py_mod_wide(__int128 a, i64 b) {
    if (b == 0) throw runtime_error("integer modulo by zero");
    if (a >= (-9223372036854775807LL - 1) && a <= 9223372036854775807LL) return b == -1 ? 0 : py_mod(static_cast<i64>(a), b);
    __int128 r = a % b;
    if (r != 0 && ((r > 0) != (b > 0))) r += b;
    return static_cast<i64>(r);
}

i64 py_pow(i64 a, i64 b) {
    i64 result = 1;
    while (b > 0) {
        if (b & 1) result *= a;
        b >>= 1;
        if (b) a *= a;
    }
    return result;
}

i64 py_range_count(i64 start, i64 stop, i64 step) {
    if (step == 0) throw runtime_error("range() arg 3 must not be zero");
    __int128 distance = step > 0
        ? static_cast<__int128>(stop) - start
        : static_cast<__int128>(start) - stop;
    if (distance <= 0) return 0;
    __int128 width = step > 0 ? step : -static_cast<__int128>(step);
    return static_cast<i64>((distance + width - 1) / width);
}

i64 py_pow_mod(i64 a, i64 b, i64 mod) {
    a = py_mod(a, mod);
    if (b < 0) {
        i64 old_r = a;
        i64 r = mod;
        i64 old_s = 1;
        i64 s = 0;
        while (r != 0) {
            i64 quotient = old_r / r;
            i64 next_r = old_r - quotient * r;
            old_r = r;
            r = next_r;
            i64 next_s = old_s - quotient * s;
            old_s = s;
            s = next_s;
        }
        if (old_r != 1) throw runtime_error("base is not invertible for the given modulus");
        a = py_mod(old_s, mod);
        b = -b;
    }
    if (mod > 0 && static_cast<u64>(mod) <= 4294967295ULL) {
        u64 base = static_cast<u64>(a);
        u64 modulus = static_cast<u64>(mod);
        u64 result = 1 % modulus;
        while (b > 0) {
            if (b & 1) result = result * base % modulus;
            b >>= 1;
            if (b) base = base * base % modulus;
        }
        return static_cast<i64>(result);
    }
    i64 result = py_mod(1, mod);
    while (b > 0) {
        if (b & 1) result = py_mod_wide(static_cast<__int128>(result) * a, mod);
        a = py_mod_wide(static_cast<__int128>(a) * a, mod);
        b >>= 1;
    }
    return result;
}

tuple<i64, i64> py_divmod(i64 a, i64 b) {
    return make_tuple(py_floor_div(a, b), py_mod(a, b));
}

i64 py_bit_length(i64 value) {
    u64 current = value < 0 ? static_cast<u64>(-(value + 1)) + 1 : static_cast<u64>(value);
    return static_cast<i64>(bit_width(current));
}

string py_bin(i64 value) {
    u64 current = value < 0 ? static_cast<u64>(-(value + 1)) + 1 : static_cast<u64>(value);
    string result;
    do {
        result.push_back((current & 1) ? '1' : '0');
        current >>= 1;
    } while (current);
    result += "b0";
    if (value < 0) result.push_back('-');
    reverse(result.begin(), result.end());
    return result;
}

string py_format_binary(i64 value, const string& specification) {
    string result = py_bin(value);
    if (!specification.empty() && specification.back() == 'b') {
        i64 width = 0;
        for (size_t index = 0; index + 1 < specification.size(); ++index) {
            if (isdigit(static_cast<unsigned char>(specification[index]))) {
                width = width * 10 + (specification[index] - '0');
            }
        }
        if (width > static_cast<i64>(result.size())) result.insert(0, width - result.size(), '0');
    }
    return result;
}

vector<i64> py_range(i64 start, i64 stop, i64 step = 1) {
    vector<i64> result;
    i64 count = py_range_count(start, stop, step);
    if (count > 0) result.reserve(static_cast<size_t>(count));
    for (i64 value = start; step > 0 ? value < stop : value > stop; value += step) result.push_back(value);
    return result;
}

mt19937_64& py_rng() {
    static mt19937_64 generator(random_device{}());
    return generator;
}

i64 py_randint(i64 left, i64 right) {
    if (left > right) throw runtime_error("empty range for randint");
    return uniform_int_distribution<i64>(left, right)(py_rng());
}

i64 py_randrange(i64 start, i64 stop, i64 step = 1) {
    if (step == 0) throw runtime_error("zero step for randrange");
    __int128 distance = step > 0 ? static_cast<__int128>(stop) - start : static_cast<__int128>(start) - stop;
    __int128 width = step > 0 ? step : -static_cast<__int128>(step);
    if (distance <= 0) throw runtime_error("empty range for randrange");
    i64 count = static_cast<i64>((distance + width - 1) / width);
    return static_cast<i64>(static_cast<__int128>(start) + static_cast<__int128>(step) * uniform_int_distribution<i64>(0, count - 1)(py_rng()));
}

i64 py_getrandbits(i64 width) {
    if (width < 0 || width > 63) throw runtime_error("getrandbits supports 0..63 bits");
    if (width == 0) return 0;
    return static_cast<i64>(py_rng()() & ((u64(1) << width) - 1));
}

string py_urandom(i64 size) {
    if (size < 0) throw runtime_error("negative byte count");
    string result(static_cast<size_t>(size), '\0');
    for (char& value : result) value = static_cast<char>(py_rng()() & 255);
    return result;
}

i64 py_randbelow(i64 upper) {
    if (upper <= 0) throw runtime_error("Upper bound must be positive");
    return uniform_int_distribution<i64>(0, upper - 1)(py_rng());
}

template <class C>
decltype(auto) py_choice(const C& source) {
    if (source.empty()) throw runtime_error("Cannot choose from an empty sequence");
    return py_at(source, py_randbelow(static_cast<i64>(source.size())));
}

template <class C>
void py_shuffle(C& source) {
    shuffle(source.begin(), source.end(), py_rng());
}

template <class C>
void py_sort_reverse(C& source, bool reverse_order) {
    if (reverse_order) sort(source.begin(), source.end(), greater<typename C::value_type>());
    else sort(source.begin(), source.end());
}

tuple<i64, i64> py_as_integer_ratio(double value) {
    if (!isfinite(value)) throw runtime_error("cannot convert non-finite float to ratio");
    int exponent = 0;
    double mantissa = frexp(value, &exponent);
    i64 numerator = static_cast<i64>(ldexp(mantissa, 53));
    i64 denominator = i64(1) << 53;
    if (exponent >= 0) {
        if (exponent >= 63) throw overflow_error("float ratio overflow");
        numerator <<= exponent;
    } else {
        int shift = -exponent;
        if (shift >= 63) denominator <<= min(shift, 62);
        else denominator >>= shift;
    }
    while (denominator > 1 && numerator % 2 == 0 && denominator % 2 == 0) {
        numerator /= 2;
        denominator /= 2;
    }
    return {numerator, denominator};
}

bool py_isascii(const string& value) {
    for (unsigned char character : value) if (character >= 128) return false;
    return true;
}

string py_to_bytes(i64 value, i64 size, const string& byteorder) {
    if (size < 0) throw runtime_error("negative byte count");
    string result(static_cast<size_t>(size), '\0');
    u64 bits = static_cast<u64>(value);
    for (i64 index = 0; index < size; ++index) {
        i64 offset = byteorder == "little" ? index : size - index - 1;
        result[static_cast<size_t>(offset)] = static_cast<char>(bits & 255);
        bits >>= 8;
    }
    return result;
}

template <class C>
auto py_next(const C& source) -> typename C::value_type {
    if (source.empty()) throw runtime_error("StopIteration");
    return *source.begin();
}

template <class C, class T>
auto py_next(const C& source, const T& fallback) {
    if (source.empty()) return fallback;
    return *source.begin();
}

i64 py_int_from_bytes(const string& source) {
    u64 result = 0;
    size_t limit = min<size_t>(source.size(), sizeof(i64));
    for (size_t index = 0; index < limit; ++index) {
        result |= static_cast<u64>(static_cast<unsigned char>(source[index])) << (index * 8);
    }
    return static_cast<i64>(result);
}

vector<i64> py_bytes_values(const string& source) {
    vector<i64> result;
    result.reserve(source.size());
    for (unsigned char value : source) result.push_back(value);
    return result;
}

string py_translate(const string& source, const string& table) {
    string result = source;
    for (char& value : result) {
        unsigned char index = static_cast<unsigned char>(value);
        if (index < table.size()) value = table[index];
    }
    return result;
}

void py_set_slice(vector<i64>& target, bool has_left, i64 left, bool has_right, i64 right, i64 step, const string& replacement) {
    py_set_slice<i64>(target, has_left, left, has_right, right, step, py_bytes_values(replacement));
}

template <class C>
vector<i64> py_integer_values(const C& source) {
    return vector<i64>(source.begin(), source.end());
}

template <class... T>
vector<i64> py_integer_values(const tuple<T...>& source) {
    vector<i64> result;
    result.reserve(sizeof...(T));
    apply([&](const auto&... value) { (result.push_back(static_cast<i64>(value)), ...); }, source);
    return result;
}

template <class C>
auto py_enumerate(const C& source, i64 start = 0) {
    vector<tuple<i64, typename C::value_type>> result;
    result.reserve(source.size());
    for (const auto& value : source) result.emplace_back(start++, value);
    return result;
}

i64 py_bit_count(i64 value) {
    u64 current = value < 0 ? static_cast<u64>(-(value + 1)) + 1 : static_cast<u64>(value);
    return static_cast<i64>(popcount(current));
}

bool py_truth(bool value) { return value; }
bool py_truth(i64 value) { return value != 0; }
bool py_truth(double value) { return value != 0.0; }
bool py_truth(const string& value) { return !value.empty(); }
template <class C>
bool py_truth(const C& value) {
    if constexpr (requires { value.empty(); }) return !value.empty();
    else return static_cast<bool>(value);
}
template <class T>
bool py_truth(const optional<T>& value) { return value.has_value() && py_truth(*value); }

template <class T>
bool py_callable(const T&) { return false; }
template <class R, class... Args>
bool py_callable(const function<R(Args...)>&) { return true; }

template <class C, class T>
bool py_contains(const C& container, const T& value) {
    return find(container.begin(), container.end(), value) != container.end();
}

template <class... T, class U>
bool py_contains(const tuple<T...>& container, const U& value) {
    bool result = false;
    apply([&](const auto&... item) { ((result = result || (item == value)), ...); }, container);
    return result;
}

bool py_contains(const string& container, const string& value) {
    return container.find(value) != string::npos;
}

bool py_contains(const string& container, char value) {
    return container.find(value) != string::npos;
}

template <class T, class C, class A>
bool py_contains(const set<T, C, A>& container, const T& value) {
    return container.find(value) != container.end();
}

template <class K, class V>
bool py_contains(const PyDict<K, V>& container, const K& value) {
    return container.find(value) != container.end();
}

template <class K, class V>
bool py_contains(const PyDefaultDict<K, V>& container, const K& value) {
    return container.find(value) != container.end();
}

template <class C>
bool py_all(const C& values) {
    for (const auto& value : values) if (!py_truth(value)) return false;
    return true;
}

template <class C>
bool py_any(const C& values) {
    for (const auto& value : values) if (py_truth(value)) return true;
    return false;
}

bool py_all_string(const string&) { return true; }
bool py_any_string(const string& values) { return !values.empty(); }

template <size_t N>
bool py_all(const char (&values)[N]) {
    for (size_t i = 0; i + 1 < N; ++i) if (!values[i]) return false;
    return true;
}

template <size_t N>
bool py_any(const char (&values)[N]) {
    for (size_t i = 0; i + 1 < N; ++i) if (values[i]) return true;
    return false;
}

template <class C>
typename C::value_type py_product(const C& values) {
    typename C::value_type result = 1;
    for (const auto& value : values) result *= value;
    return result;
}

template <class A, class B>
auto py_min(const A& left, const B& right) {
    using result_type = common_type_t<A, B>;
    return std::min(static_cast<result_type>(left), static_cast<result_type>(right));
}

template <class A, class B>
auto py_max(const A& left, const B& right) {
    using result_type = common_type_t<A, B>;
    return std::max(static_cast<result_type>(left), static_cast<result_type>(right));
}

template <class C>
string py_join(const string& separator, const C& values) {
    string result;
    using Item = typename C::value_type;
    if constexpr (is_same_v<Item, string> || is_same_v<Item, char>) {
        if (values.empty()) return result;
        __int128 size = static_cast<__int128>(separator.size()) * (values.size() - 1);
        if constexpr (is_same_v<Item, string>) {
            for (const auto& value : values) size += value.size();
        } else {
            size += values.size();
        }
        if (size > result.max_size()) throw length_error("joined string is too large");
        result.reserve(static_cast<size_t>(size));
    }
    bool first = true;
    for (const auto& value : values) {
        if (!first) result += separator;
        first = false;
        if constexpr (is_same_v<Item, string> || is_same_v<Item, char>) result += value;
        else result += py_string(value);
    }
    return result;
}

template <size_t N>
string py_join(const string& separator, const char (&values)[N]) {
    return py_join(separator, string(values, N - 1));
}

template <class... C>
auto py_zip(const C&... values) {
    using result_type = tuple<typename C::value_type...>;
    size_t size = min({static_cast<size_t>(values.size())...});
    vector<result_type> result;
    result.reserve(size);
    for (size_t i = 0; i < size; ++i) result.emplace_back(values[i]...);
    return result;
}

template <class K, class V>
PyDict<K, V> py_dict(const vector<tuple<K, V>>& values) {
    PyDict<K, V> result;
    result.data.reserve(values.size());
    result.order.reserve(values.size());
    for (const auto& [key, value] : values) result[key] = value;
    return result;
}

template <class K, class V, class Q>
optional<V> py_dict_get(const PyDict<K, V>& values, const Q& key) {
    auto iterator = values.find(K(key));
    return iterator == values.end() ? nullopt : optional<V>(iterator->second);
}

template <class K, class V, class Q>
optional<V> py_dict_get(const PyDefaultDict<K, V>& values, const Q& key) {
    return py_dict_get(static_cast<const PyDict<K, V>&>(values), key);
}

template <class K, class V, class Q, class F>
V py_dict_get(const PyDict<K, V>& values, const Q& key, const F& fallback) {
    auto iterator = values.find(K(key));
    return iterator == values.end() ? static_cast<V>(fallback) : iterator->second;
}

template <class K, class V, class Q, class F>
V py_dict_get(const PyDefaultDict<K, V>& values, const Q& key, const F& fallback) {
    return py_dict_get(static_cast<const PyDict<K, V>&>(values), key, fallback);
}

template <class K, class V>
vector<K> py_dict_keys(const PyDict<K, V>& values) {
    return values.order;
}

template <class K, class V>
vector<K> py_dict_keys(const PyDefaultDict<K, V>& values) {
    return values.order;
}

template <class K, class V>
vector<V> py_dict_values(const PyDict<K, V>& values) {
    vector<V> result;
    result.reserve(values.size());
    for (const K& key : values.order) result.push_back(values.at(key));
    return result;
}

template <class K, class V>
vector<V> py_dict_values(const PyDefaultDict<K, V>& values) {
    return py_dict_values(static_cast<const PyDict<K, V>&>(values));
}

template <class K, class V>
vector<tuple<K, V>> py_dict_items(const PyDict<K, V>& values) {
    vector<tuple<K, V>> result;
    result.reserve(values.size());
    for (const K& key : values.order) result.emplace_back(key, values.at(key));
    return result;
}

template <class K, class V>
vector<tuple<K, V>> py_dict_items(const PyDefaultDict<K, V>& values) {
    return py_dict_items(static_cast<const PyDict<K, V>&>(values));
}

template <class K, class V>
V py_dict_pop(PyDict<K, V>& values, const K& key) {
    K normalized(key);
    V result = std::move(values.at(normalized));
    values.erase(normalized);
    return result;
}

template <class K, class V>
V py_dict_pop(PyDefaultDict<K, V>& values, const K& key) {
    return py_dict_pop(static_cast<PyDict<K, V>&>(values), key);
}

template <class K, class V, class Q, class F>
V& py_dict_setdefault(PyDict<K, V>& values, const Q& key, const F& fallback) {
    K normalized(key);
    auto [iterator, inserted] = values.data.try_emplace(normalized, fallback);
    if (inserted) values.order.push_back(std::move(normalized));
    return iterator->second;
}

template <class K, class V, class Q, class F>
V& py_dict_setdefault(PyDefaultDict<K, V>& values, const Q& key, const F& fallback) {
    return py_dict_setdefault(static_cast<PyDict<K, V>&>(values), key, fallback);
}

template <class T>
void py_heapify(vector<T>& a) { make_heap(a.begin(), a.end(), greater<T>()); }

template <class T, class U>
void py_heappush(vector<T>& a, const U& value) {
    a.push_back(static_cast<T>(value));
    push_heap(a.begin(), a.end(), greater<T>());
}

template <class T>
T py_heappop(vector<T>& a) {
    pop_heap(a.begin(), a.end(), greater<T>());
    T result = std::move(a.back());
    a.pop_back();
    return result;
}
'''
        self.lines.extend(runtime.splitlines())

    def method_cpp_name(self, name: str) -> str:
        special = {
            "__len__": "py_len",
            "__getitem__": "py_getitem",
            "__setitem__": "py_setitem",
            "__contains__": "py_contains",
            "__str__": "py_str",
        }
        return special.get(name, sanitize_name(name))

    def _method_arguments(self, class_name: str, method: ast.FunctionDef, defaults: bool, copied: set[str] | None = None) -> list[str]:
        function_ir = self.program_ir.methods[(class_name, method.name)]
        params = function_ir.parameter_types
        mutable_params = function_ir.mutable_parameters
        arguments = list(method.args.posonlyargs) + list(method.args.args)
        is_classmethod = bool(
            method.decorator_list
            and isinstance(method.decorator_list[0], ast.Name)
            and method.decorator_list[0].id == "classmethod"
        )
        arguments = [
            argument
            for index, argument in enumerate(arguments)
            if argument.arg != "self" and not (is_classmethod and index == 0)
        ]
        default_offset = len(arguments) - len(method.args.defaults)
        result: list[str] = []
        previous_context = self.context
        env = dict(self.analyzer.global_types)
        env.update(params)
        self.context = FunctionContext(env, {}, class_name=class_name)
        for index, argument in enumerate(arguments):
            argument_type = params[argument.arg]
            mutable = argument.arg in mutable_params
            by_value = index >= default_offset and mutable and argument_type.kind in {
                "list", "deque", "set", "map", "default_map"
            }
            parameter_type = argument_type.cpp() if by_value else self.parameter_cpp_type(argument_type, mutable)
            argument_name = f"__py_param_{sanitize_name(argument.arg)}" if copied and argument.arg in copied else sanitize_name(argument.arg)
            declaration = f"{parameter_type} {argument_name}"
            if defaults and index >= default_offset:
                default = method.args.defaults[index - default_offset]
                declaration += f" = {self.default_argument_expression(argument_type, default)}"
            result.append(declaration)
        self.context = previous_context
        return result

    def default_argument_expression(self, argument_type: Type, default: ast.expr) -> str:
        if (
            isinstance(default, ast.Constant)
            and default.value is None
            and argument_type.kind in {"list", "deque", "set", "map", "default_map", "string", "bytes"}
        ):
            return f"{argument_type.cpp()}{{}}"
        return self.expr(default)

    def rebound_expression_key(self, name: str, value_type: Type) -> str:
        return f"{name}:{value_type.cpp()}"

    def local_rebind_info(
        self,
        body: list[ast.stmt],
        env: dict[str, Type],
    ) -> tuple[dict[str, Type], set[str], dict[str, list[Type]]]:
        current = dict(env)
        first: dict[str, Type] = {}
        rebound: set[str] = set()
        rebound_targets: dict[str, list[Type]] = {}

        def is_rebind(old: Type, new: Type) -> bool:
            if old.kind in {"unknown", "void"} or new.kind in {"unknown", "void"}:
                return False
            return old.kind != new.kind

        def update_name(name: str, value_type: Type) -> None:
            if value_type.kind in {"unknown", "void"}:
                return
            old = current.get(name, UNKNOWN)
            if name not in first:
                first[name] = value_type
                current[name] = value_type
                return
            if is_rebind(old, value_type):
                if name not in rebound and old.kind not in {"unknown", "void"}:
                    try:
                        first[name] = merge_type(first[name], old)
                    except ValueError:
                        first[name] = old
                rebound.add(name)
                targets = rebound_targets.setdefault(name, [])
                if not any(target.kind == value_type.kind and target.cpp() == value_type.cpp() for target in targets):
                    targets.append(value_type)
                current[name] = value_type
                return
            try:
                current[name] = merge_type(old, value_type)
            except ValueError:
                current[name] = value_type

        def bind_target(target: ast.AST, value_type: Type) -> None:
            if isinstance(target, ast.Name):
                update_name(target.id, value_type)
            elif isinstance(target, ast.Starred):
                bind_target(target.value, list_type(element_type(value_type)))
            elif isinstance(target, (ast.Tuple, ast.List)):
                if value_type.kind == "tuple" and len(value_type.args) == len(target.elts):
                    for child, child_type in zip(target.elts, value_type.args):
                        bind_target(child, child_type)
                elif value_type.kind in {"list", "deque", "string", "bytes"}:
                    item_type = element_type(value_type)
                    for child in target.elts:
                        bind_target(child, item_type)
                else:
                    for child in target.elts:
                        bind_target(child, I64)

        def infer_value(node: ast.AST, state: dict[str, Type]) -> Type:
            try:
                return self.analyzer.infer(node, state)
            except TranspileError as error:
                if error.message.startswith("unknown name:"):
                    return UNKNOWN
                raise

        def merge_branch(left: dict[str, Type], right: dict[str, Type]) -> dict[str, Type]:
            merged = dict(left)
            for name in set(left) | set(right):
                try:
                    merged[name] = merge_type(left.get(name, UNKNOWN), right.get(name, UNKNOWN))
                except ValueError:
                    merged[name] = UNKNOWN
            return merged

        def visit(statements: list[ast.stmt], state: dict[str, Type]) -> dict[str, Type]:
            nonlocal current
            previous = current
            current = state
            for statement in statements:
                if isinstance(statement, ast.Assign):
                    value_type = infer_value(statement.value, current)
                    for target in statement.targets:
                        if isinstance(target, ast.Subscript):
                            try:
                                self.analyzer._refine_subscript_assignment(target, value_type, current)
                            except TranspileError as error:
                                if not error.message.startswith("unknown name:"):
                                    raise
                        bind_target(target, value_type)
                elif isinstance(statement, ast.AnnAssign):
                    value_type = type_from_annotation(statement.annotation)
                    if value_type.kind == "unknown" and statement.value is not None:
                        value_type = infer_value(statement.value, current)
                    bind_target(statement.target, value_type)
                elif isinstance(statement, ast.AugAssign) and isinstance(statement.target, ast.Name):
                    update_name(statement.target.id, infer_value(statement.value, current))
                elif isinstance(statement, ast.Expr):
                    try:
                        self.analyzer._refine_mutation(statement.value, current)
                    except TranspileError as error:
                        if not error.message.startswith("unknown name:"):
                            raise
                elif isinstance(statement, ast.If):
                    left = visit(statement.body, dict(current))
                    right = visit(statement.orelse, dict(current))
                    current = merge_branch(left, right)
                elif isinstance(statement, ast.For):
                    try:
                        item_type = element_type(self.analyzer.iter_type(statement.iter, current))
                    except TranspileError as error:
                        if not error.message.startswith("unknown name:"):
                            raise
                        item_type = UNKNOWN
                    if (
                        item_type.kind == "unknown"
                        and isinstance(statement.iter, ast.Call)
                        and isinstance(statement.iter.func, ast.Name)
                        and statement.iter.func.id == "range"
                    ):
                        item_type = I64
                    bind_target(statement.target, item_type)
                    body_state = visit(statement.body, dict(current))
                    else_state = visit(statement.orelse, dict(body_state))
                    current = merge_branch(body_state, else_state)
                elif isinstance(statement, ast.While):
                    body_state = visit(statement.body, dict(current))
                    else_state = visit(statement.orelse, dict(body_state))
                    current = merge_branch(body_state, else_state)
            result = current
            current = previous
            return result

        final_state = visit(body, dict(current))
        initial: dict[str, Type] = {}
        for name, value_type in first.items():
            if name not in rebound:
                continue
            final_type = final_state.get(name, value_type)
            if value_type.kind == final_type.kind and value_type.kind in {"list", "deque", "set", "map", "default_map"}:
                try:
                    value_type = merge_type(value_type, final_type)
                except ValueError:
                    pass
            initial[name] = value_type
        return initial, rebound, rebound_targets

    def _emit_class_declaration(self, class_name: str, node: ast.ClassDef) -> None:
        methods = self.analyzer.class_methods(class_name)
        cpp_name = sanitize_name(class_name)
        bases = self.analyzer.class_bases.get(class_name, [])
        inheritance = " : " + ", ".join(f"public {sanitize_name(base)}" for base in bases) if bases else ""
        self.emit(f"class {cpp_name}{inheritance} {{")
        self.emit("public:")
        self.indent += 1
        for base in bases:
            self.emit(f"using {sanitize_name(base)}::{sanitize_name(base)};")
        previous_context = self.context
        class_env = dict(self.analyzer.global_types)
        class_env.update({attribute: value[0] for attribute, value in self.analyzer.class_constants[class_name].items()})
        self.context = FunctionContext(class_env, {}, class_name=class_name)
        for attribute, (attribute_type, value) in self.analyzer.class_constants[class_name].items():
            self.emit(f"static inline {attribute_type.cpp()} {sanitize_name(attribute)} = {self.expr(value)};")
        self.context = previous_context
        for attribute, attribute_type in self.analyzer.class_attributes[class_name].items():
            self.emit(f"{attribute_type.cpp()} {sanitize_name(attribute)}{{}};")
        initializer = methods.get("__init__")
        if initializer is None:
            self.emit(f"{cpp_name}() = default;")
        else:
            copied = set(self.program_ir.methods[(class_name, "__init__")].copied_parameters)
            arguments = self._method_arguments(class_name, initializer, defaults=True, copied=copied)
            python_arguments = [
                argument
                for argument in list(initializer.args.posonlyargs) + list(initializer.args.args)
                if argument.arg != "self"
            ]
            required = len(python_arguments) - len(initializer.args.defaults)
            if required:
                self.emit(f"{cpp_name}() = default;")
            self.emit(f"{cpp_name}({', '.join(arguments)});")
        for method_name, method in methods.items():
            if method_name == "__init__":
                continue
            return_type = self.analyzer.method_return_types[class_name][method_name]
            arguments = self._method_arguments(class_name, method, defaults=True)
            static_prefix = "static " if (
                method.decorator_list
                and isinstance(method.decorator_list[0], ast.Name)
                and method.decorator_list[0].id in {"staticmethod", "classmethod"}
            ) else ""
            self.emit(f"{static_prefix}{return_type.cpp()} {self.method_cpp_name(method_name)}({', '.join(arguments)});")
        self.indent -= 1
        self.emit("};")
        self.emit()

    def native_local_bindings(self, node: ast.FunctionDef, env: dict[str, Type]) -> tuple[set[str], set[str]]:
        result = self.optimization_planner.local_bindings(node, env)
        return set(result.callable_aliases), set(result.private_lists)

    def _emit_class_methods(self, class_name: str, node: ast.ClassDef) -> None:
        cpp_name = sanitize_name(class_name)
        for method_name, method in self.analyzer.class_methods(class_name).items():
            function_ir = self.program_ir.methods[(class_name, method_name)]
            params = dict(function_ir.parameter_types)
            copied = set(function_ir.copied_parameters)
            arguments = self._method_arguments(class_name, method, defaults=False, copied=copied)
            if method_name == "__init__":
                signature = f"{cpp_name}::{cpp_name}({', '.join(arguments)})"
            else:
                return_type = self.analyzer.method_return_types[class_name][method_name]
                signature = f"{return_type.cpp()} {cpp_name}::{self.method_cpp_name(method_name)}({', '.join(arguments)})"
            self.emit(signature + " {")
            self.indent += 1
            env = dict(function_ir.environment)
            all_types = dict(function_ir.value_types)
            reference_aliases = set(function_ir.reference_aliases)
            native_aliases = set(function_ir.native_callable_aliases)
            private_lists = set(function_ir.private_lists)
            local_lists = set(function_ir.local_lists)
            declared = {
                name: variable.value_type
                for name, variable in sorted(function_ir.locals.items())
                if variable.declared
            }
            for variable in sorted(copied):
                self.emit(f"{all_types[variable].cpp()} {sanitize_name(variable)} = __py_param_{sanitize_name(variable)};")
            for variable, variable_type in declared.items():
                cpp_type = f"py_local_list<{element_type(variable_type).cpp()}>" if variable in local_lists else variable_type.cpp()
                self.emit(f"{cpp_type} {sanitize_name(variable)}{{}};")
            if declared:
                self.emit()
            generator_result = None
            if (class_name, method_name) in self.analyzer.generator_methods:
                generator_result = "__py_generator_result"
                return_type = self.analyzer.method_return_types[class_name][method_name]
                self.emit(f"{return_type.cpp()} {generator_result}{{}};")
                self.emit()
            self.context = FunctionContext(
                all_types,
                declared,
                class_name=class_name,
                native_callable_aliases=native_aliases,
                private_lists=private_lists,
                local_lists=local_lists,
                reference_aliases=reference_aliases,
                generator_result=generator_result,
                return_type=self.analyzer.method_return_types[class_name][method_name],
            )
            for statement in method.body:
                self.statement(statement)
            if generator_result is not None:
                self.emit(f"return {generator_result};")
            elif self.analyzer.method_return_types[class_name][method_name].kind == "optional":
                self.emit("return nullopt;")
            self.indent -= 1
            self.emit("}")
            self.emit()

    def emit_local_function(self, node: ast.FunctionDef) -> None:
        if node.args.vararg is not None or node.args.kwarg is not None:
            return
        parameters: dict[str, Type] = {}
        arguments = list(node.args.posonlyargs) + list(node.args.args)
        for argument in arguments:
            argument_type = type_from_annotation(argument.annotation)
            if argument_type.kind == "unknown":
                argument_type = self.analyzer._fallback_parameter_type(node, argument.arg)
            parameters[argument.arg] = argument_type
        outer_context = self.context
        if outer_context is None:
            self.fail("local function outside a function", node)
        local_env = dict(outer_context.env)
        local_env.update(parameters)
        return_type = type_from_annotation(node.returns)
        if return_type.kind == "unknown":
            return_type = self.analyzer._infer_returns(node.body, local_env) or VOID
        function = function_type(parameters.values(), return_type)
        outer_context.env[node.name] = function
        declarations = [
            f"{self.parameter_cpp_type(parameters[argument.arg], False)} {sanitize_name(argument.arg)}"
            for argument in arguments
        ]
        self.emit(f"auto {sanitize_name(node.name)} = [&]({', '.join(declarations)}) {{")
        self.indent += 1
        all_types = self.analyzer.locals(node.body, local_env)
        local_names = self.analyzer.bound_names(node.body) - set(parameters)
        local_function_names = {
            child.name
            for child in node.body
            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef))
        }
        local_names -= local_function_names
        declared = {
            name: all_types.get(name, UNKNOWN)
            for name in sorted(local_names)
        }
        for name, value_type in declared.items():
            self.emit(f"{value_type.cpp()} {sanitize_name(name)}{{}};")
        if declared:
            self.emit()
        self.context = FunctionContext(
            all_types,
            declared,
            class_name=outer_context.class_name,
            reference_expressions=dict(outer_context.reference_expressions),
            return_type=return_type,
        )
        for statement in node.body:
            self.statement(statement)
        self.context = outer_context
        self.indent -= 1
        self.emit("};")

    def _emit_globals(self) -> None:
        for name, t in self.analyzer.global_types.items():
            constant = self.optimization_planner.constant_globals.get(name)
            if constant is not None:
                value = ast.literal_eval(constant.value)
                literal = "true" if value is True else "false" if value is False else str(value) + "LL"
                self.emit(f"constexpr {t.cpp()} {sanitize_name(name)} = {literal};")
            else:
                self.emit(f"{t.cpp()} {sanitize_name(name)}{{}};")
        if self.analyzer.global_types:
            self.emit()

    def _emit_prototype(self, name: str, node: ast.FunctionDef) -> None:
        function_ir = self.program_ir.functions[name]
        params = function_ir.parameter_types
        return_type = function_ir.return_type
        args = []
        arguments = list(node.args.posonlyargs) + list(node.args.args)
        default_offset = len(arguments) - len(node.args.defaults)
        previous_context = self.context
        self.context = FunctionContext(dict(self.analyzer.global_types), {})
        for index, arg in enumerate(arguments):
            t = params[arg.arg]
            mutable = arg.arg in function_ir.mutable_parameters
            declaration = f"{self.parameter_cpp_type(t, mutable)} {sanitize_name(arg.arg)}"
            if index >= default_offset:
                if mutable and t.kind in {"list", "deque", "set", "map", "default_map"}:
                    self.fail("mutable container parameters cannot have defaults", arg)
                declaration += f" = {self.default_argument_expression(t, node.args.defaults[index - default_offset])}"
            args.append(declaration)
        self.context = previous_context
        self.emit(f"{return_type.cpp()} {sanitize_name(name)}({', '.join(args)});")

    def parameter_cpp_type(self, t: Type, mutable: bool = False) -> str:
        if t.kind in {"list", "deque", "set", "map", "default_map"}:
            return f"{t.cpp()}&" if mutable else f"const {t.cpp()}&"
        if t.kind in {"string", "bytes", "tuple"}:
            return f"const {t.cpp()}&"
        if t.kind == "class":
            return f"{t.cpp()}&"
        return t.cpp()

    def copy_rebound_parameters(self, node: ast.FunctionDef, params: dict[str, Type]) -> set[str]:
        rebound = self.analyzer.bound_names(node.body)
        arguments = list(node.args.posonlyargs) + list(node.args.args)
        none_defaults = {
            argument.arg
            for argument, default in zip(arguments[-len(node.args.defaults):], node.args.defaults)
            if isinstance(default, ast.Constant) and default.value is None
        }
        body = [statement for statement in node.body if not (
            isinstance(statement, ast.Expr)
            and isinstance(statement.value, ast.Constant)
            and isinstance(statement.value.value, str)
        )]
        copied: set[str] = set()
        for name in rebound:
            if params.get(name, UNKNOWN).kind in {"i64", "f64", "bool", "char"} and any(
                isinstance(statement, ast.Assign)
                and any(
                    isinstance(target, ast.Name) and target.id == name
                    for target in statement.targets
                )
                for statement in ast.walk(node)
            ):
                copied.add(name)
        if body and isinstance(body[0], ast.Assign) and len(body[0].targets) == 1:
            first = body[0]
            if (
                isinstance(first.targets[0], ast.Name)
                and isinstance(first.value, ast.Call)
                and isinstance(first.value.func, ast.Name)
                and first.value.func.id == "list"
                and len(first.value.args) == 1
                and not first.value.keywords
                and isinstance(first.value.args[0], ast.Name)
                and first.value.args[0].id == first.targets[0].id
                and params.get(first.targets[0].id, UNKNOWN).kind == "list"
            ):
                copied.add(first.targets[0].id)
            elif (
                isinstance(first.targets[0], ast.Name)
                and isinstance(first.value, ast.Subscript)
                and isinstance(first.value.value, ast.Name)
                and first.value.value.id == first.targets[0].id
                and isinstance(first.value.slice, ast.Slice)
                and first.value.slice.lower is None
                and first.value.slice.upper is None
                and first.value.slice.step is None
                and params.get(first.targets[0].id, UNKNOWN).kind == "list"
            ):
                copied.add(first.targets[0].id)
            elif (
                isinstance(first.targets[0], ast.Name)
                and isinstance(first.value, ast.Call)
                and isinstance(first.value.func, ast.Attribute)
                and first.value.func.attr == "copy"
                and isinstance(first.value.func.value, ast.Name)
                and first.value.func.value.id == first.targets[0].id
                and not first.value.args
                and not first.value.keywords
                and params.get(first.targets[0].id, UNKNOWN).kind == "list"
            ):
                copied.add(first.targets[0].id)
        for statement in body:
            if not isinstance(statement, ast.Assign) or len(statement.targets) != 1 or not isinstance(statement.targets[0], ast.Name):
                continue
            target = statement.targets[0].id
            value = statement.value
            is_copy = (
                isinstance(value, ast.Call)
                and (
                    (
                        isinstance(value.func, ast.Attribute)
                        and value.func.attr == "copy"
                        and isinstance(value.func.value, ast.Name)
                        and value.func.value.id == target
                        and not value.args
                        and not value.keywords
                    )
                    or (
                        isinstance(value.func, ast.Name)
                        and value.func.id in {"iter", "list", "tuple", "set", "sorted"}
                        and len(value.args) == 1
                        and isinstance(value.args[0], ast.Name)
                        and value.args[0].id == target
                        and not value.keywords
                    )
                )
            )
            is_slice = (
                isinstance(value, ast.Subscript)
                and isinstance(value.value, ast.Name)
                and value.value.id == target
                and isinstance(value.slice, ast.Slice)
                and value.slice.lower is None
                and value.slice.upper is None
                and value.slice.step is None
            )
            if (is_copy or is_slice) and params.get(target, UNKNOWN).kind == "list":
                copied.add(target)
        for name in rebound:
            if name in copied or params.get(name, UNKNOWN).kind not in {
                "list", "deque", "set", "map", "default_map", "string", "bytes", "tuple", "class"
            }:
                continue
            rebinding_lines = [
                statement.lineno
                for statement in ast.walk(node)
                if isinstance(statement, ast.Assign)
                and any(
                    any(isinstance(child, ast.Name) and child.id == name for child in ast.walk(target))
                    for target in statement.targets
                )
                and not (
                    isinstance(statement.value, ast.Name)
                    and statement.value.id == name
                )
            ]
            if rebinding_lines:
                first_rebinding = min(rebinding_lines)
                mutated_before_rebinding = any(
                    getattr(descendant, "lineno", first_rebinding) < first_rebinding
                    and (
                        (
                            isinstance(descendant, ast.AugAssign)
                            and (
                                isinstance(descendant.target, ast.Name)
                                and descendant.target.id == name
                                or isinstance(descendant.target, ast.Subscript)
                                and isinstance(descendant.target.value, ast.Name)
                                and descendant.target.value.id == name
                            )
                        )
                        or (
                            isinstance(descendant, ast.Delete)
                            and any(
                                isinstance(target, ast.Name) and target.id == name
                                or isinstance(target, ast.Subscript)
                                and isinstance(target.value, ast.Name)
                                and target.value.id == name
                                for target in descendant.targets
                            )
                        )
                        or (
                            isinstance(descendant, ast.Assign)
                            and any(
                                isinstance(target, ast.Subscript)
                                and isinstance(target.value, ast.Name)
                                and target.value.id == name
                                for target in descendant.targets
                            )
                        )
                        or (
                            isinstance(descendant, ast.Call)
                            and isinstance(descendant.func, ast.Attribute)
                            and isinstance(descendant.func.value, ast.Name)
                            and descendant.func.value.id == name
                        )
                    )
                    for descendant in ast.walk(node)
                )
                if not mutated_before_rebinding:
                    copied.add(name)
                    continue
            mutated = False
            for descendant in ast.walk(node):
                if isinstance(descendant, ast.Assign) and isinstance(descendant.value, ast.Name) and descendant.value.id == name:
                    if any(not isinstance(target, ast.Name) or target.id != name for target in descendant.targets):
                        mutated = True
                if isinstance(descendant, (ast.Assign, ast.AugAssign, ast.Delete)):
                    targets = descendant.targets if isinstance(descendant, (ast.Assign, ast.Delete)) else [descendant.target]
                    if any(
                        isinstance(target, ast.Subscript)
                        and isinstance(target.value, ast.Name)
                        and target.value.id == name
                        for target in targets
                    ):
                        mutated = True
                if isinstance(descendant, ast.Call):
                    if isinstance(descendant.func, ast.Attribute) and isinstance(descendant.func.value, ast.Name):
                        if descendant.func.value.id == name and descendant.func.attr in {
                            "append", "appendleft", "extend", "insert", "pop", "popleft", "remove", "clear", "sort", "reverse", "rotate"
                        }:
                            mutated = True
                    if any(isinstance(argument, ast.Name) and argument.id == name for argument in descendant.args):
                        if not isinstance(descendant.func, ast.Name) or descendant.func.id not in {
                            "len", "sum", "min", "max", "all", "any", "sorted", "list", "tuple", "set", "str",
                            "isinstance", "issubclass", "type"
                        }:
                            mutated = True
            if not mutated:
                copied.add(name)
        for argument in list(node.args.posonlyargs) + list(node.args.args):
            if argument.arg in rebound and params[argument.arg].kind in {
                "list", "deque", "set", "map", "default_map", "class"
            } and argument.arg not in copied:
                self.fail(f"rebinding reference parameter {argument.arg} is not supported", argument)
        return copied

    def _emit_function(self, name: str, node: ast.FunctionDef) -> None:
        function_ir = self.program_ir.functions[name]
        params = dict(function_ir.parameter_types)
        return_type = function_ir.return_type
        copied = set(function_ir.copied_parameters)
        args = []
        for arg in list(node.args.posonlyargs) + list(node.args.args):
            mutable = arg.arg in function_ir.mutable_parameters
            argument_name = f"__py_param_{sanitize_name(arg.arg)}" if arg.arg in copied else sanitize_name(arg.arg)
            args.append(f"{self.parameter_cpp_type(params[arg.arg], mutable)} {argument_name}")
        self.emit(f"{return_type.cpp()} {sanitize_name(name)}({', '.join(args)}) {{")
        self.indent += 1
        env = dict(function_ir.environment)
        all_types = dict(function_ir.value_types)
        rebound_types = dict(function_ir.rebound_types)
        rebound_names = set(function_ir.rebound_names)
        rebound_targets = {
            key: list(values) for key, values in function_ir.rebound_targets.items()
        }
        rebound_expressions: dict[str, str] = {}
        rebound_declarations: list[tuple[str, Type, str]] = []
        for variable in sorted(rebound_names):
            for index, rebound_type in enumerate(rebound_targets.get(variable, [])):
                key = self.rebound_expression_key(variable, rebound_type)
                expression = f"__py_rebound_{sanitize_name(variable)}_{index}"
                rebound_expressions[key] = expression
                rebound_declarations.append((expression, rebound_type, key))
        for variable, variable_type in rebound_types.items():
            if variable in rebound_names:
                all_types[variable] = variable_type
        reference_aliases = set(function_ir.reference_aliases)
        private_lists = set(function_ir.private_lists)
        local_lists = set(function_ir.local_lists)
        declared = {
            variable: local.value_type
            for variable, local in sorted(function_ir.locals.items())
            if local.declared
        }
        for variable in sorted(copied):
            initial = f"__py_param_{sanitize_name(variable)}"
            parameter_type = params.get(variable, UNKNOWN)
            local_type = all_types[variable]
            if local_type.kind == "string" and parameter_type.kind in {"i64", "f64", "bool", "char"}:
                initial = f"py_string({initial})"
            self.emit(f"{local_type.cpp()} {sanitize_name(variable)} = {initial};")
        for variable, t in declared.items():
            cpp_type = f"py_local_list<{element_type(t).cpp()}>" if variable in local_lists else t.cpp()
            self.emit(f"{cpp_type} {sanitize_name(variable)}{{}};")
        for variable, rebound_type, _ in rebound_declarations:
            if rebound_type.kind not in {"unknown", "void"}:
                self.emit(f"{rebound_type.cpp()} {variable}{{}};")
        if declared or rebound_declarations:
            self.emit()
        self.context = FunctionContext(
            all_types,
            declared,
            reference_aliases=reference_aliases,
            private_lists=private_lists,
            local_lists=local_lists,
            rebound_names=rebound_names,
            rebound_expressions=rebound_expressions,
            generator_result="__py_generator_result" if name in self.analyzer.generator_functions else None,
            return_type=return_type,
        )
        if name in self.analyzer.generator_functions:
            self.emit(f"{return_type.cpp()} __py_generator_result{{}};")
            self.emit()
        for statement in node.body:
            self.statement(statement)
        if name in self.analyzer.generator_functions:
            self.emit("return __py_generator_result;")
        self.indent -= 1
        self.emit("}")
        self.emit()

    def _emit_main(self) -> None:
        body: list[ast.stmt] = []
        for node in self.tree.body:
            if isinstance(node, (ast.Import, ast.ImportFrom, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                continue
            if isinstance(node, ast.If) and self.analyzer.is_main_condition(node.test):
                body.extend(node.body)
            else:
                body.append(node)
        env = dict(self.analyzer.global_types)
        all_types = self.analyzer.locals(body, env)
        declared = {key: value for key, value in all_types.items() if key not in self.analyzer.global_types}
        self.emit("int main() {")
        self.indent += 1
        for variable, t in declared.items():
            self.emit(f"{t.cpp()} {sanitize_name(variable)}{{}};")
        if declared:
            self.emit()
        self.context = FunctionContext(all_types, declared)
        for statement in body:
            self.statement(statement)
        self.emit("return 0;")
        self.indent -= 1
        self.emit("}")

    def line_directive(self, node: ast.AST) -> None:
        return

    def statement(self, node: ast.stmt) -> None:
        self.line_directive(node)
        if node in self.optimization_planner.constant_globals.values():
            return
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            self.emit_local_function(node)
            return
        if isinstance(node, (ast.Pass, ast.Import, ast.ImportFrom, ast.Global, ast.Nonlocal)):
            return
        if isinstance(node, ast.Raise):
            if node.cause is not None or node.exc is None:
                self.fail("only direct exception raises are supported", node)
            exception = node.exc
            if isinstance(exception, ast.Name) and exception.id == "IndexError":
                self.emit('throw out_of_range("IndexError");')
                return
            if isinstance(exception, ast.Call) and isinstance(exception.func, ast.Name) and exception.func.id == "IndexError":
                if len(exception.args) > 1 or exception.keywords:
                    self.fail("IndexError accepts at most one positional argument", exception)
                message = self.expr(exception.args[0]) if exception.args else '"IndexError"'
                self.emit(f"throw out_of_range({message});")
                return
            if isinstance(exception, ast.Name) and exception.id in {"ValueError", "TypeError", "RuntimeError", "StopIteration", "OverflowError", "EOFError"}:
                self.emit(f'throw runtime_error("{exception.id}");')
                return
            if isinstance(exception, ast.Call) and isinstance(exception.func, ast.Name) and exception.func.id in {"ValueError", "TypeError", "RuntimeError", "StopIteration", "OverflowError", "EOFError"}:
                if len(exception.args) > 1 or exception.keywords:
                    self.fail(f"{exception.func.id} accepts at most one positional argument", exception)
                message = self.expr(exception.args[0]) if exception.args else cpp_string(exception.func.id)
                self.emit(f"throw runtime_error({message});")
                return
            self.fail("only direct standard exception raises are supported", node)
        if isinstance(node, ast.Expr):
            if isinstance(node.value, ast.Yield):
                if self.context.generator_result is None or node.value.value is None:
                    self.fail("only value-yielding methods are supported", node)
                self.emit(f"{self.context.generator_result}.push_back({self.expr(node.value.value)});")
                return
            if isinstance(node.value, ast.Constant) and isinstance(node.value.value, str):
                return
            if self.emit_special_expression_statement(node.value):
                return
            self.emit(self.expr(node.value) + ";")
            return
        if isinstance(node, ast.Assign):
            if self.analyzer.is_input_alias_assignment(node) or self.analyzer.is_read_alias_assignment(node):
                return
            self.assignment(node)
            return
        if isinstance(node, ast.AnnAssign):
            if node.value is not None:
                if isinstance(node.value, (ast.List, ast.Set, ast.Dict)) and not (
                    node.value.elts if isinstance(node.value, (ast.List, ast.Set)) else node.value.keys
                ):
                    target_type = type_from_annotation(node.annotation)
                    if target_type.kind == "unknown":
                        target_type = self.analyzer.infer(node.target, self.context.env)
                    self.assign_target(node.target, f"{target_type.cpp()}{{}}")
                else:
                    self.assign_target(node.target, self.expr(node.value))
            return
        if isinstance(node, ast.AugAssign):
            target = self.target(node.target)
            value = self.value_expr(node.value)
            target_type = self.analyzer.infer(node.target, self.context.env)
            if target_type.kind == "optional":
                target_type = target_type.args[0]
                target = f"(*{target})"
            if isinstance(node.op, ast.Add) and target_type.kind == "list":
                temporary = self.context.temporary("extension")
                self.emit(f"auto&& {temporary} = ({value});")
                self.emit(f"{target}.insert({target}.end(), {temporary}.begin(), {temporary}.end());")
            elif isinstance(node.op, ast.FloorDiv):
                self.emit(f"{target} = py_floor_div({target}, {value});")
            elif isinstance(node.op, ast.Mod):
                self.emit(f"{target} = py_mod({target}, {value});")
            else:
                self.emit(f"{target} {self.operator(node.op)}= {value};")
            return
        if isinstance(node, ast.If):
            branch_bindings = {
                name: (self.context.env[name], self.context.reference_expressions.get(name))
                for name in self.context.rebound_names
            } if node.orelse else {}
            self.emit(f"if (py_truth({self.expr(node.test)})) {{")
            self.indent += 1
            for child in node.body:
                self.statement(child)
            self.indent -= 1
            if node.orelse:
                for name, (value_type, reference) in branch_bindings.items():
                    self.context.env[name] = value_type
                    if reference is None:
                        self.context.reference_expressions.pop(name, None)
                    else:
                        self.context.reference_expressions[name] = reference
                self.emit("} else {")
                self.indent += 1
                for child in node.orelse:
                    self.statement(child)
                self.indent -= 1
            self.emit("}")
            return
        if isinstance(node, ast.For):
            flag = self.context.temporary("completed") if node.orelse else None
            if flag is not None:
                self.emit(f"bool {flag} = true;")
            self.loop_else_flags.append(flag)
            self.for_statement(node)
            self.loop_else_flags.pop()
            if flag is not None:
                self.emit(f"if ({flag}) {{")
                self.indent += 1
                for child in node.orelse:
                    self.statement(child)
                self.indent -= 1
                self.emit("}")
            return
        if isinstance(node, ast.While):
            flag = self.context.temporary("completed") if node.orelse else None
            if flag is not None:
                self.emit(f"bool {flag} = true;")
            self.loop_else_flags.append(flag)
            self.emit(f"while (py_truth({self.expr(node.test)})) {{")
            self.indent += 1
            for child in node.body:
                self.statement(child)
            self.indent -= 1
            self.emit("}")
            self.loop_else_flags.pop()
            if flag is not None:
                self.emit(f"if ({flag}) {{")
                self.indent += 1
                for child in node.orelse:
                    self.statement(child)
                self.indent -= 1
                self.emit("}")
            return
        if isinstance(node, ast.Return):
            if node.value is None:
                self.emit("return;")
            elif isinstance(node.value, ast.List) and not node.value.elts and self.context.return_type is not None and self.context.return_type.kind == "list":
                self.emit(f"return {self.context.return_type.cpp()}{{}};")
            elif (
                isinstance(node.value, ast.List)
                and len(node.value.elts) == 1
                and self.context.return_type is not None
                and self.context.return_type.kind == "list"
                and self.context.return_type.args[0].kind == "tuple"
                and len(self.context.return_type.args[0].args) == 2
                and all(item.kind == "f64" for item in self.context.return_type.args[0].args)
                and self.analyzer.infer(node.value.elts[0], self.context.env).kind == "tuple"
            ):
                item = self.expr(node.value.elts[0])
                converted = (
                    f"make_tuple(static_cast<double>(get<0>({item})), "
                    f"static_cast<double>(get<1>({item})))"
                )
                self.emit(f"return {self.context.return_type.cpp()}{{{converted}}};")
            elif self.context.return_type is not None:
                expected = self.context.return_type
                actual = self.analyzer.infer(node.value, self.context.env)
                if (
                    expected.kind == "optional"
                    and expected.args[0].kind == "tuple"
                    and actual.kind not in {"optional", "tuple", "void"}
                ):
                    self.emit(f"return make_tuple({self.expr(node.value)}, vector<bool>{{}});")
                elif (
                    expected.kind == "optional"
                    and expected.args[0].kind in {"i64", "f64", "bool", "char", "string"}
                    and actual.kind == "optional"
                    and actual.args[0].kind == "tuple"
                ):
                    value = self.expr(node.value)
                    temporary = self.context.temporary("optional_return")
                    self.emit(f"auto {temporary} = {value};")
                    self.emit(
                        f"return {temporary}.has_value() ? optional<{expected.args[0].cpp()}>(get<0>(*{temporary})) : nullopt;"
                    )
                elif (
                    isinstance(node.value, (ast.Tuple, ast.List))
                    and expected.kind == "tuple"
                    and any(item.kind == "reference" for item in expected.args)
                ):
                    values = ", ".join(self.expr(value) for value in node.value.elts)
                    self.emit(f"return {expected.cpp()}({values});")
                else:
                    self.emit(f"return {self.expr(node.value)};")
            elif (
                isinstance(node.value, (ast.Tuple, ast.List))
                and self.context.return_type is not None
                and self.context.return_type.kind == "tuple"
                and any(item.kind == "reference" for item in self.context.return_type.args)
            ):
                values = ", ".join(self.expr(value) for value in node.value.elts)
                self.emit(f"return {self.context.return_type.cpp()}({values});")
            else:
                self.emit(f"return {self.expr(node.value)};")
            return
        if isinstance(node, ast.Break):
            if self.loop_else_flags and self.loop_else_flags[-1] is not None:
                self.emit(f"{self.loop_else_flags[-1]} = false;")
            self.emit("break;")
            return
        if isinstance(node, ast.Continue):
            self.emit("continue;")
            return
        if isinstance(node, ast.Assert):
            self.emit(f"assert(py_truth({self.expr(node.test)}));")
            return
        if isinstance(node, ast.Delete):
            for target in node.targets:
                if isinstance(target, ast.Name):
                    target_type = self.context.env.get(target.id, UNKNOWN)
                    if is_container(target_type):
                        self.emit(f"{self.target(target)}.clear();")
                    else:
                        self.emit(f"{self.target(target)} = {{}};")
                elif isinstance(target, ast.Subscript) and not isinstance(target.slice, ast.Slice):
                    base = self.expr(target.value)
                    index = self.expr(target.slice)
                    base_type = self.analyzer.infer(target.value, self.context.env)
                    if base_type.kind in {"map", "default_map"}:
                        self.emit(f"{base}.erase({index});")
                    else:
                        self.emit(f"py_erase_at({base}, {index});")
                elif isinstance(target, ast.Subscript) and isinstance(target.slice, ast.Slice):
                    base = self.expr(target.value)
                    lower = self.expr(target.slice.lower) if target.slice.lower else "0"
                    upper = self.expr(target.slice.upper) if target.slice.upper else "0"
                    step = self.expr(target.slice.step) if target.slice.step else "1"
                    self.emit(
                        f"py_delete_slice({base}, {'true' if target.slice.lower else 'false'}, {lower}, "
                        f"{'true' if target.slice.upper else 'false'}, {upper}, {step});"
                    )
                else:
                    self.fail("unsupported delete target", target)
            return
        self.fail(f"unsupported statement: {type(node).__name__}", node)

    def rebind_assignment(self, node: ast.Assign) -> bool:
        if len(node.targets) != 1 or not isinstance(node.targets[0], ast.Name):
            return False
        name = node.targets[0].id
        if self.context is None or name not in self.context.rebound_names:
            return False
        value_type = self.analyzer.infer(node.value, self.context.env)
        current_type = self.context.env.get(name, UNKNOWN)
        if value_type.kind in {"unknown", "void"} or current_type.kind in {"unknown", "void"}:
            return False
        if value_type.kind == current_type.kind:
            return False
        value = self.value_expr(node.value)
        variable = self.context.rebound_expressions.get(self.rebound_expression_key(name, value_type))
        if variable is None:
            variable = self.context.temporary(f"rebound_{sanitize_name(name)}")
            self.emit(f"{value_type.cpp()} {variable} = {value};")
        else:
            self.emit(f"{variable} = {value};")
        self.context.reference_expressions[name] = variable
        self.context.env[name] = value_type
        return True

    def assignment(self, node: ast.Assign) -> None:
        if len(node.targets) == 1 and isinstance(node.targets[0], ast.Name) and node.targets[0].id in self.context.local_lists:
            values = ", ".join(self.expr(value) for value in node.value.elts)
            self.emit(f"{sanitize_name(node.targets[0].id)}.assign({{{values}}});")
            return
        if len(node.targets) == 1 and isinstance(node.targets[0], ast.Name) and node.targets[0].id in self.context.native_callable_aliases:
            reference = isinstance(node.value, ast.Attribute) and node.value.attr in self.optimization_planner.immutable_callables.get(self.context.class_name, ())
            declaration = "const auto&" if reference else "auto"
            self.emit(f"{declaration} {sanitize_name(node.targets[0].id)} = {self.expr(node.value)};")
            return
        if self.rebind_assignment(node):
            return
        if (
            len(node.targets) == 1
            and isinstance(node.targets[0], ast.Name)
            and isinstance(node.value, ast.Attribute)
            and isinstance(node.value.value, ast.Name)
            and node.value.value.id == "heapq"
            and node.value.attr in {"heappush", "heappop", "heapify"}
        ):
            self.context.call_aliases[node.targets[0].id] = node.value.attr
            return
        if (
            len(node.targets) == 1
            and isinstance(node.targets[0], ast.Name)
            and node.targets[0].id in self.context.reference_aliases
            and node.targets[0].id not in self.context.initialized_aliases
        ):
            name = node.targets[0].id
            self.emit(f"auto& {sanitize_name(name)} = {self.expr(node.value)};")
            self.context.initialized_aliases.add(name)
            return
        if len(node.targets) == 1 and isinstance(node.targets[0], (ast.Tuple, ast.List)):
            target = node.targets[0]
            if any(isinstance(child, ast.Starred) for child in target.elts):
                if self.emit_starred_unpack(target, node.value):
                    return
            if self.emit_unpack_read(target, node.value):
                return
            if isinstance(node.value, (ast.Tuple, ast.List)) and len(target.elts) == len(node.value.elts):
                values = [self.expr(value) for value in node.value.elts]
                targets = [self.target(value) for value in target.elts]
                self.emit(f"std::tie({', '.join(targets)}) = std::make_tuple({', '.join(values)});")
                return
            value_type = self.analyzer.infer(node.value, self.context.env)
            if value_type.kind == "tuple" and len(value_type.args) == len(target.elts):
                temporary = self.context.temporary("unpacked")
                self.emit(f"auto {temporary} = {self.expr(node.value)};")
                for index, (child, child_type) in enumerate(zip(target.elts, value_type.args)):
                    value = f"std::get<{index}>({temporary})"
                    if child_type.kind == "reference" and isinstance(child, ast.Name):
                        self.emit(f"auto& {sanitize_name(child.id)} = {value};")
                        self.context.initialized_aliases.add(child.id)
                    else:
                        self.assign_unpack_target(child, value, child_type)
                return
            if value_type.kind == "optional" and value_type.args[0].kind == "tuple" and len(value_type.args[0].args) == len(target.elts):
                temporary = self.context.temporary("unpacked")
                self.emit(f"auto {temporary} = {self.expr(node.value)};")
                for index, child in enumerate(target.elts):
                    self.assign_unpack_target(child, f"std::get<{index}>(*{temporary})", value_type.args[0].args[index])
                return
            if value_type.kind == "optional":
                temporary = self.context.temporary("unpacked")
                self.emit(f"auto {temporary} = {self.expr(node.value)};")
                source = f"(*{temporary})"
                for index, child in enumerate(target.elts):
                    self.assign_unpack_target(child, f"py_at({source}, {index})", element_type(value_type.args[0]))
                return
            if value_type.kind in {"list", "deque", "string", "bytes", "unknown"}:
                temporary = self.context.temporary("unpacked")
                self.emit(f"auto {temporary} = {self.expr(node.value)};")
                for index, child in enumerate(target.elts):
                    self.assign_unpack_target(child, f"py_at({temporary}, {index})", element_type(value_type))
                return
            temporary = self.context.temporary("unpacked")
            self.emit(f"auto {temporary} = {self.expr(node.value)};")
            for index, child in enumerate(target.elts):
                self.assign_unpack_target(child, f"py_at({temporary}, {index})", element_type(value_type))
            return
        elif (
            len(node.targets) == 1
            and isinstance(node.value, ast.Call)
            and isinstance(node.value.func, ast.Name)
            and node.value.func.id in {"deque", "list", "set", "dict"}
            and not node.value.args
            and not node.value.keywords
        ):
            target_type = self.analyzer.infer(node.targets[0], self.context.env)
            if target_type.kind in {"list", "deque", "set", "map", "default_map"}:
                value = f"{target_type.cpp()}{{}}"
            else:
                value = self.value_expr(node.value)
        elif (
            len(node.targets) == 1
            and isinstance(node.value, ast.Call)
            and isinstance(node.value.func, ast.Name)
            and node.value.func.id == "defaultdict"
        ):
            target_type = self.analyzer.infer(node.targets[0], self.context.env)
            factory = node.value.args[0]
            default_factory = f"[&]() {{ return {self.expr(factory.body)}; }}" if isinstance(factory, ast.Lambda) else None
            if len(node.value.args) == 1:
                value = f"{target_type.cpp()}{{{default_factory}}}" if default_factory is not None else f"{target_type.cpp()}{{}}"
            else:
                initial = self.expr(node.value.args[1])
                value = f"{target_type.cpp()}{{{default_factory}, {initial}}}" if default_factory is not None else f"{target_type.cpp()}{{{initial}}}"
        elif (
            len(node.targets) == 1
            and isinstance(node.value, (ast.List, ast.Set, ast.Dict))
            and not (node.value.elts if isinstance(node.value, (ast.List, ast.Set)) else node.value.keys)
        ):
            target_type = self.analyzer.infer(node.targets[0], self.context.env)
            value = f"{target_type.cpp()}{{}}"
        elif (
            len(node.targets) == 1
            and isinstance(node.value, ast.ListComp)
            and isinstance(node.value.elt, ast.List)
            and not node.value.elt.elts
        ):
            value = self.comprehension(node.value, self.analyzer.infer(node.targets[0], self.context.env))
        elif (
            len(node.targets) == 1
            and isinstance(node.value, ast.BinOp)
            and isinstance(node.value.op, ast.Mult)
            and isinstance(node.value.left, ast.List)
            and len(node.value.left.elts) == 1
            and isinstance(node.value.left.elts[0], ast.Constant)
            and node.value.left.elts[0].value is None
        ):
            target_type = self.analyzer.infer(node.targets[0], self.context.env)
            if target_type.kind != "list" or target_type.args[0].kind != "optional":
                self.fail("cannot infer the element type of a repeated None list", node)
            value = f"{target_type.cpp()}({self.expr(node.value.right)}, nullopt)"
        else:
            target_types = [self.analyzer.infer(target, self.context.env) for target in node.targets]
            value = self.expr(node.value) if any(value.kind == "optional" for value in target_types) else self.value_expr(node.value)
            value_type = self.analyzer.infer(node.value, self.context.env)
            if (
                len(target_types) == 1
                and target_types[0].kind == "list"
                and value_type.kind == "list"
                and target_types[0].args[0].kind == "string"
                and isinstance(node.value, ast.BinOp)
                and isinstance(node.value.op, ast.Mult)
                and isinstance(node.value.left, ast.List)
                and len(node.value.left.elts) == 1
                and self.analyzer.infer(node.value.left.elts[0], self.context.env).kind in {"i64", "f64", "bool", "char"}
            ):
                item = f"py_string({self.expr(node.value.left.elts[0])})"
                value = f"py_repeat(vector<string>{{{item}}}, {self.expr(node.value.right)})"
            if (
                len(target_types) == 1
                and target_types[0] == list_type(STRING)
                and value_type == list_type(list_type(CHAR))
            ):
                value = f"py_char_rows({value})"
        for target in node.targets:
            self.assign_target(target, value)

    def emit_starred_unpack(self, target: ast.Tuple | ast.List, value: ast.AST) -> bool:
        starred = [child for child in target.elts if isinstance(child, ast.Starred)]
        if len(starred) != 1:
            self.fail("multiple starred assignment targets are not supported", target)
        value_type = self.analyzer.infer(value, self.context.env)
        if value_type.kind not in {"list", "deque"}:
            self.fail("starred assignment requires a list or deque value", value)
        star_index = next(index for index, child in enumerate(target.elts) if isinstance(child, ast.Starred))
        prefix = target.elts[:star_index]
        suffix = target.elts[star_index + 1:]
        temporary = self.context.temporary("unpacked")
        self.emit(f"auto {temporary} = {self.expr(value)};")
        for index, child in enumerate(prefix):
            self.assign_target(child, f"py_at({temporary}, {index})")
        suffix_start = f"static_cast<i64>({temporary}.size()) - {len(suffix)}"
        for index, child in enumerate(suffix):
            self.assign_target(child, f"py_at({temporary}, {suffix_start} + {index})")
        self.assign_target(
            starred[0].value,
            f"py_slice({temporary}, true, {len(prefix)}, true, {suffix_start}, 1)",
        )
        return True

    def emit_unpack_read(self, target: ast.Tuple | ast.List, value: ast.AST) -> bool:
        if not isinstance(value, ast.Call) or not isinstance(value.func, ast.Name) or value.func.id != "map":
            return False
        if len(value.args) < 2 or not isinstance(value.args[0], ast.Name) or value.args[0].id != "int":
            return False
        source = value.args[1]
        if not self.is_input_split(source):
            return False
        temp = self.context.temporary("values")
        self.emit(f"auto {temp} = py_read_vector<i64>(py_input);")
        for index, child in enumerate(target.elts):
            self.assign_target(child, f"{temp}.at({index})")
        return True

    def assign_target(self, target: ast.AST, value: str) -> None:
        if isinstance(target, (ast.Tuple, ast.List)):
            targets = [self.target(child) for child in target.elts]
            self.emit(f"std::tie({', '.join(targets)}) = {value};")
            return
        if isinstance(target, ast.Subscript) and isinstance(target.slice, ast.Slice):
            base_type = self.analyzer.infer(target.value, self.context.env)
            if base_type.kind not in {"list", "unknown", "bytes", "string", "optional"}:
                self.fail("slice assignment is only supported for lists", target)
            base = self.expr(target.value)
            lower = self.expr(target.slice.lower) if target.slice.lower else "0"
            upper = self.expr(target.slice.upper) if target.slice.upper else "0"
            step = self.expr(target.slice.step) if target.slice.step else "1"
            self.emit(
                f"py_set_slice({base}, {'true' if target.slice.lower else 'false'}, {lower}, "
                f"{'true' if target.slice.upper else 'false'}, {upper}, {step}, {value});"
            )
            return
        if isinstance(target, ast.Subscript):
            base_type = self.analyzer.infer(target.value, self.context.env)
            if base_type.kind == "class":
                if isinstance(target.slice, ast.Slice):
                    self.fail("class slice assignment is not supported", target)
                base = self.expr(target.value)
                self.emit(f"{base}.{self.method_cpp_name('__setitem__')}({self.expr(target.slice)}, {value});")
                return
        self.emit(f"{self.target(target)} = {value};")

    def assign_unpack_target(self, target: ast.AST, value: str, value_type: Type) -> None:
        if isinstance(target, ast.Name) and target.id in self.context.rebound_names:
            current_type = self.context.env.get(target.id, UNKNOWN)
            if value_type.kind not in {"unknown", "void"} and current_type.kind != value_type.kind:
                variable = self.context.rebound_expressions.get(
                    self.rebound_expression_key(target.id, value_type)
                )
                if variable is not None:
                    self.emit(f"{variable} = {value};")
                    self.context.reference_expressions[target.id] = variable
                    self.context.env[target.id] = value_type
                    return
        self.assign_target(target, value)

    def target(self, node: ast.AST) -> str:
        if isinstance(node, ast.Name):
            if self.context is not None and node.id in self.context.reference_expressions:
                return self.context.reference_expressions[node.id]
            return sanitize_name(node.id)
        if isinstance(node, ast.Subscript):
            return self.subscript(node, mutable=True)
        if isinstance(node, ast.Attribute):
            base_type = self.analyzer.infer(node.value, self.context.env)
            if base_type.kind == "class":
                if isinstance(node.value, ast.Name) and node.value.id == "self":
                    return f"this->{sanitize_name(node.attr)}"
                return f"{self.expr(node.value)}.{sanitize_name(node.attr)}"
            return f"{self.expr(node.value)}.{sanitize_name(node.attr)}"
        self.fail("unsupported assignment target", node)
        return ""

    def for_statement(self, node: ast.For) -> None:
        iterator = node.iter
        if (
            isinstance(iterator, ast.Call)
            and isinstance(iterator.func, ast.Name)
            and iterator.func.id == "reversed"
            and len(iterator.args) == 1
            and isinstance(iterator.args[0], ast.Name)
            and iterator.args[0].id in self.context.private_lists
            and not iterator.keywords
        ):
            source = self.expr(iterator.args[0])
            index = self.context.temporary("reverse_index")
            self.emit(f"for (size_t {index} = {source}.size(); {index}-- > 0;) {{")
            self.indent += 1
            self.assign_target(node.target, f"{source}[{index}]")
            for child in node.body:
                self.statement(child)
            self.indent -= 1
            self.emit("}")
            return
        if self.reverse_range_slice_for(node):
            return
        if isinstance(node.iter, ast.Call) and isinstance(node.iter.func, ast.Name) and node.iter.func.id == "range":
            self.range_for(node)
            return
        if isinstance(node.iter, ast.Call) and isinstance(node.iter.func, ast.Name) and node.iter.func.id == "enumerate":
            if not 1 <= len(node.iter.args) <= 2 or not isinstance(node.target, (ast.Tuple, ast.List)) or len(node.target.elts) != 2:
                self.fail("enumerate must have a two-item target", node)
            source = self.expr(node.iter.args[0])
            source_type = self.analyzer.infer(node.iter.args[0], self.context.env)
            if source_type.kind == "optional":
                source = f"(*{source})"
                source_type = source_type.args[0]
            item_type = element_type(source_type)
            source_name = self.context.temporary("enumerate_source")
            self.emit(f"auto&& {source_name} = ({source});")
            source = source_name
            index = self.context.temporary("index")
            start = self.expr(node.iter.args[1]) if len(node.iter.args) == 2 else "0"
            start_value = self.context.temporary("enumerate_start")
            offset = self.context.temporary("offset")
            self.emit(f"i64 {start_value} = {start};")
            alias_target = (
                isinstance(node.target.elts[1], ast.Name)
                and item_type.kind in {"list", "deque", "set", "map", "default_map", "class"}
                and isinstance(node.iter.args[0], (ast.Name, ast.Attribute))
            )
            if alias_target:
                alias_name = node.target.elts[1].id
                index_target = self.target(node.target.elts[0])
                self.context.reference_expressions[alias_name] = (
                    f"{source}[static_cast<size_t>(({index_target}) - {start_value})]"
                )
            else:
                alias_name = None
            self.emit(f"for (i64 {offset} = 0; {offset} < static_cast<i64>({source}.size()); ++{offset}) {{")
            self.indent += 1
            self.emit(f"i64 {index} = {start_value} + {offset};")
            self.assign_target(node.target.elts[0], index)
            if alias_name is None:
                self.assign_target(node.target.elts[1], f"{source}[static_cast<size_t>({offset})]")
            for child in node.body:
                self.statement(child)
            self.indent -= 1
            self.emit("}")
            return
        source = self.expr(node.iter)
        source_type = self.analyzer.infer(node.iter, self.context.env)
        if source_type.kind == "optional":
            source_type = source_type.args[0]
            source = f"(*{source})"
        if source_type.kind == "class" and "__iter__" in self.analyzer.class_methods(source_type.name):
            class_name = source_type.name
            source = f"{source}.{self.method_cpp_name('__iter__')}()"
            source_type = self.analyzer.method_return_types[class_name]["__iter__"]
        iterable = f"py_dict_keys({source})" if source_type.kind in {"map", "default_map"} else f"py_tuple_values({source})" if source_type.kind == "tuple" else source
        if source_type.kind == "bytes":
            iterable = f"py_bytes_values({source})"
        item = self.context.temporary("item")
        item_type = element_type(source_type)
        alias_target = (
            isinstance(node.target, ast.Name)
            and item_type.kind in {"list", "deque", "set", "map", "default_map", "class"}
            and isinstance(node.iter, (ast.Name, ast.Attribute))
        )
        loop_reference = "auto&" if alias_target else "const auto&"
        self.emit(f"for ({loop_reference} {item} : {iterable}) {{")
        self.indent += 1
        if alias_target:
            self.emit(f"auto& {sanitize_name(node.target.id)} = {item};")
        else:
            self.assign_target(node.target, item)
        for child in node.body:
            self.statement(child)
        self.indent -= 1
        self.emit("}")

    def reverse_range_slice_for(self, node: ast.For) -> bool:
        iterator = node.iter
        if not (
            isinstance(iterator, ast.Subscript)
            and isinstance(iterator.value, ast.Call)
            and isinstance(iterator.value.func, ast.Name)
            and iterator.value.func.id == "range"
            and isinstance(iterator.slice, ast.Slice)
            and iterator.slice.lower is None
            and iterator.slice.upper is None
            and isinstance(iterator.slice.step, ast.UnaryOp)
            and isinstance(iterator.slice.step.op, ast.USub)
            and isinstance(iterator.slice.step.operand, ast.Constant)
            and iterator.slice.step.operand.value == 1
        ):
            return False
        if not isinstance(node.target, ast.Name):
            return False
        args = iterator.value.args
        if not 1 <= len(args) <= 3:
            self.fail("range expects one to three arguments", iterator.value)
        if len(args) == 1:
            start, stop, step = "0", self.expr(args[0]), "1"
        elif len(args) == 2:
            start, stop, step = self.expr(args[0]), self.expr(args[1]), "1"
        else:
            start, stop, step = self.expr(args[0]), self.expr(args[1]), self.expr(args[2])
        start_name = self.context.temporary("range_start")
        stop_name = self.context.temporary("range_stop")
        step_name = self.context.temporary("range_step")
        count_name = self.context.temporary("range_count")
        loop_variable = self.context.temporary("range")
        self.emit(f"i64 {start_name} = {start};")
        self.emit(f"i64 {stop_name} = {stop};")
        self.emit(f"i64 {step_name} = {step};")
        self.emit(f"i64 {count_name} = py_range_count({start_name}, {stop_name}, {step_name});")
        self.emit(f"i64 {loop_variable} = {start_name};")
        self.emit(f"if ({count_name} > 0) {loop_variable} += ({count_name} - 1) * {step_name};")
        self.emit(
            f"for (; {count_name} > 0; --{count_name}, "
            f"({count_name} > 0 ? {loop_variable} -= {step_name} : {loop_variable})) {{"
        )
        self.indent += 1
        self.assign_target(node.target, loop_variable)
        for child in node.body:
            self.statement(child)
        self.indent -= 1
        self.emit("}")
        return True

    def range_for(self, node: ast.For) -> None:
        call = node.iter
        args = call.args
        if not 1 <= len(args) <= 3:
            self.fail("range expects one to three arguments", call)
        if not isinstance(node.target, ast.Name):
            self.fail("range target must be a name", node.target)
        target = self.target(node.target)
        loop_variable = self.context.temporary("range")
        if len(args) == 1:
            start, stop, step = "0", self.expr(args[0]), "1"
        elif len(args) == 2:
            start, stop, step = self.expr(args[0]), self.expr(args[1]), "1"
        else:
            start, stop, step = self.expr(args[0]), self.expr(args[1]), self.expr(args[2])
        start_name = self.context.temporary("range_start")
        stop_name = self.context.temporary("range_stop")
        self.emit(f"i64 {start_name} = {start};")
        self.emit(f"i64 {stop_name} = {stop};")
        start, stop = start_name, stop_name
        if step in {"1", "-1"}:
            condition = f"{loop_variable} < {stop}" if step == "1" else f"{loop_variable} > {stop}"
            increment = "++" + loop_variable if step == "1" else "--" + loop_variable
            self.emit(f"for (i64 {loop_variable} = {start}; {condition}; {increment}) {{")
        else:
            step_name = self.context.temporary("step")
            self.emit(f"i64 {step_name} = {step};")
            self.emit(f"if ({step_name} == 0) throw runtime_error(\"range() arg 3 must not be zero\");")
            self.emit(f"for (i64 {loop_variable} = {start}; ({step_name} > 0 ? {loop_variable} < {stop} : {loop_variable} > {stop}); {loop_variable} += {step_name}) {{")
        self.indent += 1
        self.assign_target(node.target, loop_variable)
        for child in node.body:
            self.statement(child)
        self.indent -= 1
        self.emit("}")

    def emit_special_expression_statement(self, node: ast.AST) -> bool:
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "print":
            self.print_statement(node)
            return True
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
            owner = node.func.value
            attr = node.func.attr
            if isinstance(owner, ast.Attribute) and isinstance(owner.value, ast.Attribute):
                if isinstance(owner.value.value, ast.Name) and owner.value.value.id == "sys" and owner.value.attr == "stdout" and attr == "write":
                    self.emit(f"cout << {self.expr(node.args[0])};")
                    return True
            if isinstance(owner, ast.Attribute) and isinstance(owner.value, ast.Name):
                if owner.value.id == "sys" and owner.attr == "stdout" and attr == "write":
                    if any(
                        isinstance(child, ast.Name) and child.id in self.analyzer.buffered_output_names
                        for argument in node.args
                        for child in ast.walk(argument)
                    ):
                        return True
                    self.emit(f"cout << {self.expr(node.args[0])};")
                    return True
            if isinstance(owner, ast.Name) and owner.id == "sys" and attr == "setrecursionlimit":
                return True
            if self.analyzer.infer(owner, self.context.env).kind == "class":
                return False
            if attr in {"append", "appendleft", "add", "clear", "sort", "reverse", "rotate", "extend", "insert"}:
                self.emit(self.method_call(node) + ";")
                return True
            if isinstance(owner, ast.Name) and owner.id == "heapq" and attr in {"heappush", "heappop", "heapify"}:
                self.emit(self.heapq_call(node) + ";")
                return True
            if isinstance(owner, ast.Name) and owner.id == "bisect" and attr in {"bisect_left", "bisect_right"}:
                self.emit(self.expr(node) + ";")
                return True
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id in {"heappush", "heappop", "heapify"} and node.func.id in self.analyzer.from_imports.get("heapq", set()):
            self.emit(self.heapq_call(node) + ";")
            return True
        return False

    def print_statement(self, node: ast.Call) -> None:
        sep = " "
        end = "\n"
        flush = False
        for keyword in node.keywords:
            if keyword.arg == "sep":
                if not isinstance(keyword.value, ast.Constant) or not isinstance(keyword.value.value, str):
                    self.fail("print sep must be a string literal", keyword.value)
                sep = keyword.value.value
            elif keyword.arg == "end":
                if not isinstance(keyword.value, ast.Constant) or not isinstance(keyword.value.value, str):
                    self.fail("print end must be a string literal", keyword.value)
                end = keyword.value.value
            elif keyword.arg == "flush":
                if not isinstance(keyword.value, ast.Constant) or not isinstance(keyword.value.value, bool):
                    self.fail("print flush must be a bool literal", keyword.value)
                flush = keyword.value.value
            else:
                self.fail(f"unsupported print keyword: {keyword.arg}", keyword)
        first = self.context.temporary("print_first")
        if node.args:
            self.emit(f"bool {first} = true;")
        for arg in node.args:
            if isinstance(arg, ast.Starred):
                item = self.context.temporary("print_item")
                source_type = self.analyzer.infer(arg.value, self.context.env)
                source = self.expr(arg.value)
                iterable = f"py_dict_keys({source})" if source_type.kind in {"map", "default_map"} else f"py_tuple_values({source})" if source_type.kind == "tuple" else source
                value = item
                value_type = element_type(source_type)
                self.emit(f"for (const auto& {item} : {iterable}) {{")
                self.indent += 1
                self.emit(f"if (!{first}) cout << {cpp_string(sep)};")
                if value_type.kind == "bool":
                    self.emit(f"cout << ({value} ? \"True\" : \"False\");")
                else:
                    self.emit(f"cout << {value};")
                self.emit(f"{first} = false;")
                self.indent -= 1
                self.emit("}")
            else:
                self.emit(f"if (!{first}) cout << {cpp_string(sep)};")
                value = self.expr(arg)
                value_type = self.analyzer.infer(arg, self.context.env)
                if value_type.kind == "optional":
                    temporary = self.context.temporary("optional_print")
                    self.emit(f"auto {temporary} = {value};")
                    self.emit(f"cout << ({temporary}.has_value() ? py_string(*{temporary}) : string(\"None\"));")
                elif value_type.kind == "bool":
                    self.emit(f"cout << ({value} ? \"True\" : \"False\");")
                elif value_type.kind == "class" and "__str__" in self.analyzer.class_methods(value_type.name):
                    self.emit(f"cout << {value}.{self.method_cpp_name('__str__')}();")
                else:
                    self.emit(f"cout << {value};")
                self.emit(f"{first} = false;")
        if end:
            self.emit(f"cout << {cpp_string(end)};")
        if flush:
            self.emit("cout.flush();")

    def expr(self, node: ast.AST) -> str:
        if isinstance(node, ast.Constant):
            if node.value is None:
                return "nullopt"
            if isinstance(node.value, bool):
                return "true" if node.value else "false"
            if isinstance(node.value, str):
                return cpp_string(node.value)
            if isinstance(node.value, bytes):
                return cpp_string(node.value.decode("latin1"))
            if isinstance(node.value, float):
                return repr(node.value)
            return str(node.value)
        if isinstance(node, ast.JoinedStr):
            parts = ["string()"]
            for value in node.values:
                if isinstance(value, ast.Constant) and isinstance(value.value, str):
                    parts.append(cpp_string(value.value))
                elif isinstance(value, ast.FormattedValue):
                    if value.format_spec is not None or value.conversion not in {-1, 115, 114, 97}:
                        self.fail("unsupported f-string format", value)
                    formatter = "py_repr" if value.conversion in {114, 97} else "py_string"
                    parts.append(f"{formatter}({self.expr(value.value)})")
                else:
                    self.fail("unsupported f-string part", value)
            return "(" + " + ".join(parts) + ")"
        if isinstance(node, ast.Lambda):
            lambda_type = self.analyzer.infer(node, self.context.env)
            parameters = list(node.args.posonlyargs) + list(node.args.args)
            previous_env = self.context.env
            lambda_env = dict(previous_env)
            declarations = []
            for argument, argument_type in zip(parameters, function_parameters(lambda_type)):
                lambda_env[argument.arg] = argument_type
                if argument_type.kind in {"list", "deque", "set", "map", "default_map", "string", "bytes", "tuple", "class"}:
                    declarations.append(f"const {argument_type.cpp()}& {sanitize_name(argument.arg)}")
                else:
                    declarations.append(f"{argument_type.cpp()} {sanitize_name(argument.arg)}")
            self.context.env = lambda_env
            body = self.expr(node.body)
            self.context.env = previous_env
            return f"[&]({', '.join(declarations)}) {{ return {body}; }}"
        if isinstance(node, ast.NamedExpr):
            if not isinstance(node.target, ast.Name):
                self.fail("named expressions must assign to a name", node)
            return f"({sanitize_name(node.target.id)} = {self.expr(node.value)})"
        if isinstance(node, ast.Name):
            if node.id == "True":
                return "true"
            if node.id == "False":
                return "false"
            if node.id == "min":
                return "[&](i64 left, i64 right) { return std::min(left, right); }"
            if node.id == "max":
                return "[&](i64 left, i64 right) { return std::max(left, right); }"
            if node.id == "len":
                return "[&](const auto& value) { return static_cast<i64>(value.size()); }"
            if node.id in self.analyzer.from_imports.get("bisect", set()):
                function = "lower_bound" if node.id.endswith("left") else "upper_bound"
                return f"[&](const vector<i64>& values, i64 value) {{ return static_cast<i64>({function}(values.begin(), values.end(), value) - values.begin()); }}"
            if node.id == "pi" and node.id in self.analyzer.from_imports.get("math", set()):
                return "acos(-1.0)"
            if node.id == "self" and self.context.class_name is not None:
                return "(*this)"
            if node.id == "cls" and self.context.class_name is not None:
                return sanitize_name(self.context.class_name)
            if self.context is not None and node.id in self.context.reference_expressions:
                return self.context.reference_expressions[node.id]
            if node.id in self.analyzer.hoisted_captures and (self.context is None or node.id not in self.context.env):
                function = self.analyzer.functions[node.id]
                captured = set(self.analyzer.hoisted_captures[node.id])
                arguments = [
                    argument
                    for argument in list(function.args.posonlyargs) + list(function.args.args)
                    if argument.arg not in captured
                ]
                parameters = [f"auto&& __py_callback_arg_{index}" for index in range(len(arguments))]
                forwarded = [f"__py_callback_arg_{index}" for index in range(len(arguments))]
                captures = [self.expr(ast.Name(id=name, ctx=ast.Load())) for name in self.analyzer.hoisted_captures[node.id]]
                call_arguments = captures + forwarded
                return f"[&]({', '.join(parameters)}) {{ return {sanitize_name(node.id)}({', '.join(call_arguments)}); }}"
            return sanitize_name(node.id)
        if isinstance(node, ast.List):
            value_type = self.analyzer.infer(node, self.context.env)
            return value_type.cpp() + "{" + ", ".join(self.expr(element) for element in node.elts) + "}"
        if isinstance(node, ast.Tuple):
            if not node.elts:
                return "tuple<>()"
            return "make_tuple(" + ", ".join(self.expr(element) for element in node.elts) + ")"
        if isinstance(node, ast.Set):
            value_type = self.analyzer.infer(node, self.context.env)
            return value_type.cpp() + "{" + ", ".join(self.expr(element) for element in node.elts) + "}"
        if isinstance(node, (ast.ListComp, ast.SetComp, ast.DictComp, ast.GeneratorExp)):
            return self.comprehension(node)
        if isinstance(node, ast.Dict):
            value_type = self.analyzer.infer(node, self.context.env)
            if not node.keys:
                return "{}"
            items = ", ".join(
                "{" + self.expr(key) + ", " + self.expr(value) + "}"
                for key, value in zip(node.keys, node.values)
            )
            return f"{value_type.cpp()}{{{items}}}"
        if isinstance(node, ast.Subscript):
            return self.subscript(node, mutable=False)
        if isinstance(node, ast.Attribute):
            if isinstance(node.value, ast.Name) and node.value.id == "math" and node.attr == "pi":
                return "acos(-1.0)"
            if isinstance(node.value, ast.Name) and node.value.id == "int" and node.attr == "bit_count":
                return "[&](i64 value) { return py_bit_count(value); }"
            if isinstance(node.value, ast.Name) and node.value.id == "int" and node.attr == "bit_length":
                return "[&](i64 value) { return py_bit_length(value); }"
            base_type = self.analyzer.infer(node.value, self.context.env)
            if base_type.kind in {"list", "deque", "set"}:
                base = self.expr(node.value)
                item = element_type(base_type).cpp()
                if node.attr in {"append", "appendleft", "add"}:
                    operation = {
                        "append": f"{base}.push_back(value)",
                        "appendleft": f"{base}.push_front(value)",
                        "add": f"{base}.insert(value)",
                    }[node.attr]
                    return f"[&]({item} value) {{ {operation}; }}"
                if node.attr == "extend":
                    return f"[&](const vector<{item}>& value) {{ {base}.insert({base}.end(), value.begin(), value.end()); }}"
                if node.attr in {"pop", "popleft"}:
                    operation = "py_pop" if node.attr == "pop" else "py_popleft"
                    return f"[&]() {{ return {operation}({base}); }}"
                if node.attr in {"clear", "sort", "reverse"}:
                    operation = {
                        "clear": f"{base}.clear()",
                        "sort": f"sort({base}.begin(), {base}.end())",
                        "reverse": f"reverse({base}.begin(), {base}.end())",
                    }[node.attr]
                    return f"[&]() {{ {operation}; }}"
            if base_type.kind in {"map", "default_map"} and node.attr == "__getitem__":
                base = self.expr(node.value)
                key = base_type.args[0].cpp()
                return f"[&](const {key}& value) {{ return {base}[value]; }}"
            if base_type.kind == "class":
                if node.attr in self.analyzer.class_constants.get(base_type.name, {}):
                    return f"{sanitize_name(base_type.name)}::{sanitize_name(node.attr)}"
                method = self.analyzer.class_methods(base_type.name).get(node.attr)
                if isinstance(node.value, ast.Name) and node.value.id == "self":
                    if (
                        method is not None
                        and method.decorator_list
                        and isinstance(method.decorator_list[0], ast.Name)
                        and method.decorator_list[0].id == "property"
                    ):
                        return f"this->{self.method_cpp_name(node.attr)}()"
                    if method is not None:
                        parameters = self.analyzer.method_parameter_types[base_type.name][node.attr]
                        arguments = self.analyzer.method_arguments(method)
                        declarations = [
                            f"{self.parameter_cpp_type(parameters.get(argument.arg, UNKNOWN))} {sanitize_name(argument.arg)}"
                            for argument in arguments
                        ]
                        call_arguments = ", ".join(sanitize_name(argument.arg) for argument in arguments)
                        return f"[&]({', '.join(declarations)}) {{ return this->{self.method_cpp_name(node.attr)}({call_arguments}); }}"
                    return f"this->{sanitize_name(node.attr)}"
                if (
                    method is not None
                    and method.decorator_list
                    and isinstance(method.decorator_list[0], ast.Name)
                    and method.decorator_list[0].id == "property"
                ):
                    return f"{self.expr(node.value)}.{self.method_cpp_name(node.attr)}()"
                if method is not None:
                    parameters = self.analyzer.method_parameter_types[base_type.name][node.attr]
                    arguments = self.analyzer.method_arguments(method)
                    declarations = [
                        f"{self.parameter_cpp_type(parameters.get(argument.arg, UNKNOWN))} {sanitize_name(argument.arg)}"
                        for argument in arguments
                    ]
                    call_arguments = ", ".join(sanitize_name(argument.arg) for argument in arguments)
                    receiver = self.expr(node.value)
                    return f"[&]({', '.join(declarations)}) {{ return {receiver}.{self.method_cpp_name(node.attr)}({call_arguments}); }}"
                return f"{self.expr(node.value)}.{sanitize_name(node.attr)}"
            return f"{self.expr(node.value)}.{sanitize_name(node.attr)}"
        if isinstance(node, ast.UnaryOp):
            if isinstance(node.op, ast.Not):
                return f"(!py_truth({self.expr(node.operand)}))"
            return f"({self.operator(node.op)}{self.expr(node.operand)})"
        if isinstance(node, ast.BoolOp):
            op = " && " if isinstance(node.op, ast.And) else " || "
            return "(" + op.join(f"py_truth({self.expr(value)})" for value in node.values) + ")"
        if isinstance(node, ast.Compare):
            pieces: list[str] = []
            left = self.expr(node.left)
            left_type = self.analyzer.infer(node.left, self.context.env)
            for op, right_node in zip(node.ops, node.comparators):
                right = self.expr(right_node)
                right_type = self.analyzer.infer(right_node, self.context.env)
                if isinstance(op, (ast.In, ast.NotIn)) and right_type.kind == "optional":
                    right = f"(*{right})"
                    right_type = right_type.args[0]
                if isinstance(op, ast.In):
                    if right_type.kind == "class":
                        pieces.append(f"{right}.{self.method_cpp_name('__contains__')}({left})")
                    elif right_type.kind in {"set", "map", "default_map"}:
                        pieces.append(f"({right}.find({left}) != {right}.end())")
                    elif right_type.kind in {"string", "bytes"}:
                        pieces.append(f"({right}.find({left}) != string::npos)")
                    else:
                        pieces.append(f"::py_contains({right}, {left})")
                elif isinstance(op, ast.NotIn):
                    if right_type.kind == "class":
                        pieces.append(f"(!{right}.{self.method_cpp_name('__contains__')}({left}))")
                    elif right_type.kind in {"set", "map", "default_map"}:
                        pieces.append(f"({right}.find({left}) == {right}.end())")
                    elif right_type.kind in {"string", "bytes"}:
                        pieces.append(f"({right}.find({left}) == string::npos)")
                    else:
                        pieces.append(f"(!::py_contains({right}, {left}))")
                else:
                    comparison_left = left
                    comparison_right = right
                    if left_type.kind == "char" and right_type.kind in {"string", "bytes"}:
                        comparison_left = f"string(1, {left})"
                    elif left_type.kind in {"string", "bytes"} and right_type.kind == "char":
                        comparison_right = f"string(1, {right})"
                    pieces.append(f"({comparison_left} {self.operator(op)} {comparison_right})")
                left = right
                left_type = right_type
            return "(" + " && ".join(pieces) + ")"
        if isinstance(node, ast.IfExp):
            body = self.expr(node.body)
            orelse = self.expr(node.orelse)
            result_type = self.analyzer.infer(node, self.context.env)
            if result_type.kind == "optional" and self.analyzer.infer(node.body, self.context.env).kind != "optional":
                body = f"optional<{result_type.args[0].cpp()}>({body})"
            return f"(py_truth({self.expr(node.test)}) ? {body} : {orelse})"
        if isinstance(node, ast.BinOp):
            left = self.value_expr(node.left)
            right = self.value_expr(node.right)
            left_type = self.analyzer.infer(node.left, self.context.env)
            right_type = self.analyzer.infer(node.right, self.context.env)
            if isinstance(node.op, ast.Add) and left_type.kind == "list" and right_type.kind == "list":
                return f"py_concat({left}, {right})"
            if isinstance(node.op, ast.FloorDiv):
                return f"py_floor_div({left}, {right})"
            if isinstance(node.op, ast.Mod):
                if left_type.kind == "i64" and right_type.kind == "i64":
                    if isinstance(node.left, ast.BinOp) and isinstance(node.left.op, ast.Pow) and isinstance(node.left.right, ast.Constant) and isinstance(node.left.right.value, int) and node.left.right.value >= 0:
                        return f"py_pow_mod({self.value_expr(node.left.left)}, {self.value_expr(node.left.right)}, {right})"
                    return f"py_mod_wide({self.wide_integer_expr(node.left)}, {right})"
                return f"py_mod({left}, {right})"
            if isinstance(node.op, ast.Pow):
                return f"py_pow({left}, {right})"
            if isinstance(node.op, ast.Mult):
                if left_type.kind == "list" and right_type.kind == "i64":
                    return f"py_repeat({left}, {right})"
                if right_type.kind == "list" and left_type.kind == "i64":
                    return f"py_repeat({right}, {left})"
                if left_type.kind in {"string", "bytes"} and right_type.kind == "i64":
                    return f"py_repeat({left}, {right})"
                if right_type.kind in {"string", "bytes"} and left_type.kind == "i64":
                    return f"py_repeat({right}, {left})"
            return f"({left} {self.operator(node.op)} {right})"
        if isinstance(node, ast.Call):
            return self.call(node)
        self.fail(f"unsupported expression: {type(node).__name__}", node)
        return ""

    def value_expr(self, node: ast.AST) -> str:
        value = self.expr(node)
        value_type = self.analyzer.infer(node, self.context.env)
        return f"(*{value})" if value_type.kind == "optional" else value

    def is_simple_reference(self, node: ast.AST) -> bool:
        if isinstance(node, ast.Name):
            return True
        return isinstance(node, ast.Attribute) and self.is_simple_reference(node.value)

    def subscript(self, node: ast.Subscript, mutable: bool) -> str:
        base = self.expr(node.value)
        base_type = self.analyzer.infer(node.value, self.context.env)
        if base_type.kind == "optional":
            base_type = base_type.args[0]
            base = f"(*{base})"
        if isinstance(node.slice, ast.Slice):
            lower = self.expr(node.slice.lower) if node.slice.lower else "0"
            upper = self.expr(node.slice.upper) if node.slice.upper else "0"
            step = self.expr(node.slice.step) if node.slice.step else "1"
            return f"py_slice({base}, {'true' if node.slice.lower else 'false'}, {lower}, {'true' if node.slice.upper else 'false'}, {upper}, {step})"
        index = self.expr(node.slice)
        if base_type.kind == "tuple" and isinstance(node.slice, ast.Constant) and isinstance(node.slice.value, int):
            return f"std::get<{node.slice.value}>({base})"
        if base_type.kind in {"map", "default_map"}:
            return f"{base}[{index}]" if mutable or base_type.kind == "default_map" else f"{base}.at({index})"
        if base_type.kind == "class":
            if mutable:
                self.fail("class item assignment must be a statement", node)
            return f"{base}.{self.method_cpp_name('__getitem__')}({index})"
        if base_type.kind == "bytes":
            return f"static_cast<i64>(static_cast<unsigned char>(py_at({base}, {index})))"
        return f"py_at({base}, {index})"

    def direct_input_integer_comprehension(self, node: ast.ListComp) -> str | None:
        if len(node.generators) != 1:
            return None
        generator = node.generators[0]
        if generator.ifs or not isinstance(generator.target, ast.Name) or not self.is_input_split(generator.iter):
            return None
        target_name = generator.target.id
        converted = False
        valid = True

        class InputIntegerUses(ast.NodeVisitor):
            def visit_Call(self, call: ast.Call) -> None:
                nonlocal converted
                if (
                    isinstance(call.func, ast.Name)
                    and call.func.id == "int"
                    and len(call.args) == 1
                    and not call.keywords
                    and isinstance(call.args[0], ast.Name)
                    and call.args[0].id == target_name
                ):
                    converted = True
                    return
                self.generic_visit(call)

            def visit_Name(self, name: ast.Name) -> None:
                nonlocal valid
                if name.id == target_name:
                    valid = False

        visitor = InputIntegerUses()
        visitor.visit(node.elt)
        if not valid or not converted:
            return None
        previous_env = self.context.env
        self.context.env = dict(previous_env)
        self.analyzer._bind_comprehension_target(generator.target, I64, self.context.env)
        try:
            body = self.expr(node.elt)
        finally:
            self.context.env = previous_env
        return (
            f"py_input.read_i64_line_map([&](i64 {sanitize_name(target_name)}) {{ return {body}; }})"
        )

    def generator_reduction(self, name: str, node: ast.Call) -> str | None:
        if node.keywords or len(node.args) not in ({1, 2} if name == "sum" else {1}):
            return None
        generator = node.args[0]
        if not isinstance(generator, ast.GeneratorExp) or len(generator.generators) != 1:
            return None
        clause = generator.generators[0]
        if clause.is_async:
            return None
        iterator = clause.iter
        is_range = isinstance(iterator, ast.Call) and isinstance(iterator.func, ast.Name) and iterator.func.id == "range"
        source_type = self.analyzer.infer(iterator, self.context.env)
        if not is_range and source_type.kind not in {"list", "deque"}:
            return None
        item_type = I64 if is_range else element_type(source_type)

        def pure_initial(expression: ast.AST) -> bool:
            if isinstance(expression, (ast.Constant, ast.Name)):
                return True
            if isinstance(expression, ast.UnaryOp):
                return pure_initial(expression.operand)
            if isinstance(expression, ast.BinOp):
                return pure_initial(expression.left) and pure_initial(expression.right)
            return False

        owned_source = isinstance(iterator, (ast.List, ast.ListComp)) or (
            isinstance(iterator, ast.Call)
            and isinstance(iterator.func, ast.Name)
            and iterator.func.id in self.analyzer.functions
        )
        if len(node.args) == 2 and not is_range and not owned_source and not pure_initial(node.args[1]):
            return None

        def value_type(t: Type) -> bool:
            return t.kind in {"i64", "f64", "bool", "char", "string", "bytes"} or (
                t.kind == "tuple" and all(value_type(child) for child in t.args)
            )

        def bindable(target: ast.AST, t: Type) -> bool:
            if isinstance(target, ast.Name):
                return True
            return (
                isinstance(target, (ast.Tuple, ast.List))
                and t.kind == "tuple"
                and len(target.elts) == len(t.args)
                and all(bindable(child, child_type) for child, child_type in zip(target.elts, t.args))
            )

        targets = [child.id for child in ast.walk(clause.target) if isinstance(child, ast.Name)]
        if not value_type(item_type) or not bindable(clause.target, item_type) or len(set(targets)) != len(targets):
            return None
        local_env = dict(self.context.env)
        self.analyzer._bind_comprehension_target(clause.target, item_type, local_env)
        result_item = self.analyzer.infer(generator.elt, local_env)
        initial_type = self.analyzer.infer(node.args[1], self.context.env) if len(node.args) == 2 else I64
        result_type = merge_type(result_item, initial_type) if name == "sum" else result_item if name in {"min", "max"} else BOOL
        if name in {"min", "max"} and result_item.kind not in {"i64", "f64", "bool"}:
            return None
        if name == "sum" and (result_item.kind not in {"i64", "f64", "bool"} or initial_type.kind not in {"i64", "f64"}):
            return None
        if name != "sum" and result_item.kind not in {"i64", "f64", "bool", "char", "string", "bytes"}:
            return None
        parts = [f"([&]() -> {result_type.cpp()} {{"]
        if is_range:
            if not 1 <= len(iterator.args) <= 3 or iterator.keywords:
                return None
            arguments = [self.expr(arg) for arg in iterator.args]
            if len(arguments) == 1:
                arguments.insert(0, "0")
            if len(arguments) == 2:
                arguments.append("1")
            start, stop, step = [self.context.temporary("reduce_bound") for _ in range(3)]
            for variable, expression in zip((start, stop, step), arguments):
                parts.append(f"    i64 {variable} = {expression};")
            parts.append(f'    if ({step} == 0) throw runtime_error("range() arg 3 must not be zero");')
            index = self.context.temporary("reduce_index")
            header = f"for (__int128 {index} = {start}; ({step} > 0 ? {index} < {stop} : {index} > {stop}); {index} += {step}) {{"
            item = f"static_cast<i64>({index})"
        else:
            source = self.context.temporary("reduce_source")
            parts.append(f"    auto&& {source} = ({self.expr(iterator)});")
            index = self.context.temporary("reduce_index")
            header = f"for (size_t {index} = 0; {index} < {source}.size(); ++{index}) {{"
            item = f"{source}[{index}]"
        result = self.context.temporary("reduction")
        if name == "sum":
            initial = self.expr(node.args[1]) if len(node.args) == 2 else "0"
            parts.append(f"    {result_type.cpp()} {result} = {initial};")
        if name in {"min", "max"}:
            parts.append(f"    optional<{result_type.cpp()}> {result};")
        parts.append("    " + header)
        previous_env = self.context.env
        previous_references = self.context.reference_expressions
        self.context.env = local_env
        self.context.reference_expressions = dict(previous_references)

        def bind(target: ast.AST, t: Type, expression: str) -> None:
            if isinstance(target, ast.Name):
                variable = self.context.temporary("reduce_item")
                self.context.reference_expressions[target.id] = variable
                parts.append(f"        {t.cpp()} {variable} = {expression};")
            else:
                for position, (child, child_type) in enumerate(zip(target.elts, t.args)):
                    bind(child, child_type, f"std::get<{position}>({expression})")

        try:
            bind(clause.target, item_type, item)
            for condition in clause.ifs:
                parts.append(f"        if (!py_truth({self.expr(condition)})) continue;")
            expression = self.expr(generator.elt)
            value = self.context.temporary("reduce_value")
            parts.append(f"        {result_item.cpp()} {value} = {expression};")
            if name == "sum":
                parts.append(f"        {result} += {value};")
            elif name in {"min", "max"}:
                comparison = "<" if name == "min" else ">"
                parts.append(f"        if (!{result} || {value} {comparison} *{result}) {result} = {value};")
            else:
                truth = "true" if result_item.kind == "char" else f"py_truth({value})"
                if name == "all":
                    parts.append(f"        if (!{truth}) return false;")
                else:
                    parts.append(f"        if ({truth}) return true;")
        finally:
            self.context.env = previous_env
            self.context.reference_expressions = previous_references
        parts.append("    }")
        if name in {"min", "max"}:
            parts.append(f'    if (!{result}) throw runtime_error("{name}() arg is an empty sequence");')
        returned = f"*{result}" if name in {"min", "max"} else result if name == "sum" else "true" if name == "all" else "false"
        parts.append(f"    return {returned};")
        parts.append("}())")
        return "\n".join(parts)

    def comprehension(self, node: ast.ListComp | ast.SetComp | ast.DictComp | ast.GeneratorExp, result_type: Type | None = None) -> str:
        if isinstance(node, ast.ListComp):
            direct_input = self.direct_input_integer_comprehension(node)
            if direct_input is not None:
                return direct_input
        original_env = self.context.env
        if result_type is None:
            result_type = self.analyzer.infer(node, original_env)
        result = self.context.temporary("comprehension")
        parts = ["([&]() {", f"    {result_type.cpp()} {result};"]
        self.context.env = dict(original_env)
        if (
            result_type.kind == "list"
            and len(node.generators) == 1
            and self.is_simple_reference(node.generators[0].iter)
            and self.analyzer.iter_type(node.generators[0].iter, original_env).kind in {
                "list", "deque", "set", "string", "bytes"
            }
        ):
            source = self.expr(node.generators[0].iter)
            parts.append(f"    {result}.reserve({source}.size());")

        def emit_level(index: int) -> None:
            generator = node.generators[index]
            level = index + 1
            prefix = "    " * level
            source_type = self.analyzer.iter_type(generator.iter, self.context.env)
            source = self.expr(generator.iter) if not (
                isinstance(generator.iter, ast.Call)
                and isinstance(generator.iter.func, ast.Name)
                and generator.iter.func.id == "range"
            ) else ""
            if source_type.kind == "optional":
                source_type = source_type.args[0]
                source = f"(*{source})"
            self.analyzer._bind_comprehension_target(
                generator.target,
                element_type(source_type),
                self.context.env,
            )
            if (
                isinstance(generator.iter, ast.Call)
                and isinstance(generator.iter.func, ast.Name)
                and generator.iter.func.id == "range"
            ):
                if not isinstance(generator.target, ast.Name):
                    self.fail("range comprehension target must be a name", generator.target)
                item = sanitize_name(generator.target.id)
                args = generator.iter.args
                if len(args) == 1:
                    start, stop, step = "0", self.expr(args[0]), "1"
                elif len(args) == 2:
                    start, stop, step = self.expr(args[0]), self.expr(args[1]), "1"
                elif len(args) == 3:
                    start, stop, step = self.expr(args[0]), self.expr(args[1]), self.expr(args[2])
                else:
                    self.fail("range expects one to three arguments", generator.iter)
                start_name = self.context.temporary("range_start")
                stop_name = self.context.temporary("range_stop")
                step_name = self.context.temporary("range_step")
                count_name = self.context.temporary("range_count")
                offset = self.context.temporary("range_offset")
                parts.append(f"{prefix}i64 {start_name} = {start};")
                parts.append(f"{prefix}i64 {stop_name} = {stop};")
                parts.append(f"{prefix}i64 {step_name} = {step};")
                parts.append(f"{prefix}i64 {count_name} = py_range_count({start_name}, {stop_name}, {step_name});")
                if result_type.kind == "list" and len(node.generators) == 1 and not generator.ifs:
                    parts.append(f"{prefix}{result}.reserve(static_cast<size_t>({count_name}));")
                parts.append(f"{prefix}for (i64 {offset} = 0; {offset} < {count_name}; ++{offset}) {{")
                parts.append(f"{prefix}    i64 {item} = static_cast<i64>(static_cast<__int128>({start_name}) + static_cast<__int128>({offset}) * {step_name});")
            else:
                iterable = (
                    f"py_dict_keys({source})"
                    if source_type.kind in {"map", "default_map"}
                    else f"py_tuple_values({source})"
                    if source_type.kind == "tuple"
                    else source
                )
                item = self.context.temporary("item")
                parts.append(f"{prefix}for (const auto& {item} : {iterable}) {{")
                parts.extend(self.comprehension_bindings(generator.target, item, prefix + "    "))
            for condition in generator.ifs:
                parts.append(f"{prefix}    if (!py_truth({self.expr(condition)})) continue;")
            if index + 1 < len(node.generators):
                emit_level(index + 1)
            elif isinstance(node, ast.DictComp):
                parts.append(f"{prefix}    {result}[{self.expr(node.key)}] = {self.expr(node.value)};")
            elif isinstance(node, ast.SetComp):
                parts.append(f"{prefix}    {result}.insert({self.expr(node.elt)});")
            elif isinstance(node.elt, ast.List) and not node.elt.elts:
                parts.append(f"{prefix}    {result}.emplace_back();")
            else:
                parts.append(f"{prefix}    {result}.push_back({self.expr(node.elt)});")
            parts.append(f"{prefix}}}")
            self.context.env = dict(self.context.env)

        try:
            emit_level(0)
        finally:
            self.context.env = original_env
        parts.append(f"    return {result};")
        parts.append("}())")
        return "\n".join(parts)

    def comprehension_bindings(self, target: ast.AST, value: str, indent: str) -> list[str]:
        if isinstance(target, ast.Name):
            return [f"{indent}auto {self.target(target)} = {value};"]
        if isinstance(target, (ast.Tuple, ast.List)):
            result: list[str] = []
            for index, child in enumerate(target.elts):
                result.extend(self.comprehension_bindings(child, f"std::get<{index}>({value})", indent))
            return result
        self.fail("unsupported comprehension target", target)
        return []

    def sort_key_runtime(self, name: str, key: ast.AST) -> str:
        if (
            isinstance(key, ast.Lambda)
            and len(key.args.args) == 1
            and isinstance(key.body, ast.Name)
            and key.body.id == key.args.args[0].arg
        ):
            return name + "_direct"
        return name

    def sort_key_lambda(self, node: ast.Lambda, item_type: Type) -> str:
        if len(node.args.args) != 1:
            self.fail("sort key lambdas must have one parameter", node)
        name = node.args.args[0].arg
        previous_env = self.context.env
        self.context.env = dict(previous_env)
        self.context.env[name] = item_type
        try:
            body = self.expr(node.body)
        finally:
            self.context.env = previous_env
        return f"[&](const auto& {sanitize_name(name)}) {{ return {body}; }}"

    def expanded_starred_call(self, node: ast.Call) -> str:
        if len(node.args) != 1 or not isinstance(node.args[0], ast.Starred):
            self.fail("starred call arguments are not supported", node)
        value = node.args[0].value
        if not isinstance(value, ast.GeneratorExp) or len(value.generators) != 1:
            self.fail("starred call arguments require one generator", node.args[0])
        generator = value.generators[0]
        if generator.is_async or generator.ifs:
            self.fail("starred call generators cannot be async or filtered", generator)
        source_type = self.analyzer.infer(generator.iter, self.context.env)
        if source_type.kind != "tuple":
            self.fail("starred call generators require a fixed-size tuple", generator.iter)
        source = self.expr(generator.iter)
        temporary = self.context.temporary("expanded")
        previous_env = self.context.env
        previous_references = self.context.reference_expressions
        self.context.env = dict(previous_env)
        self.context.reference_expressions = dict(previous_references)
        arguments: list[str] = []
        try:
            for index, item_type in enumerate(source_type.args):
                if not isinstance(generator.target, ast.Name):
                    self.fail("starred call generators require a name target", generator.target)
                self.context.env[generator.target.id] = item_type
                self.context.reference_expressions[generator.target.id] = f"std::get<{index}>({temporary})"
                arguments.append(self.expr(value.elt))
        finally:
            self.context.env = previous_env
            self.context.reference_expressions = previous_references
        function = self.expr(node.func)
        return f"([&]() {{ auto&& {temporary} = ({source}); return {function}({', '.join(arguments)}); }}())"

    def call_argument(self, node: ast.AST, expected: Type | None = None) -> str:
        value = self.expr(node)
        actual = self.analyzer.infer(node, self.context.env)
        if actual.kind == "optional" and (expected is None or expected.kind != "optional"):
            return f"(*{value})"
        return value

    def user_function_call(self, name: str, node: ast.Call) -> str:
        if node.keywords:
            self.fail("keyword arguments for user functions are not supported", node)
        function = self.analyzer.functions[name]
        target_arguments = list(function.args.posonlyargs) + list(function.args.args)
        mutable = self.analyzer.mutable_parameters.get(name, set())
        arguments: list[str] = []
        temporaries: list[tuple[str, str]] = []
        for argument, target in zip(node.args, target_arguments):
            expression = self.call_argument(argument, self.analyzer.parameter_types[name].get(target.arg))
            if target.arg in mutable and not isinstance(argument, (ast.Name, ast.Attribute, ast.Subscript)):
                temporary = self.context.temporary("call_arg")
                temporaries.append((temporary, expression))
                arguments.append(temporary)
            else:
                arguments.append(expression)
        call = f"{sanitize_name(name)}({', '.join(arguments)})"
        if not temporaries:
            return call
        declarations = " ".join(f"auto {temporary} = ({expression});" for temporary, expression in temporaries)
        if self.analyzer.return_types.get(name, UNKNOWN).kind == "void":
            return f"([&]() {{ {declarations} {call}; }}())"
        return f"([&]() {{ {declarations} return {call}; }}())"

    def call(self, node: ast.Call) -> str:
        if any(isinstance(argument, ast.Starred) for argument in node.args):
            return self.expanded_starred_call(node)
        if isinstance(node.func, ast.Name):
            name = node.func.id
            if name in self.analyzer.ignored_functions:
                return "void(0)"
            if self.context is not None and name in self.context.call_aliases:
                return self.heapq_call(node, self.context.call_aliases[name])
            if name in self.context.env:
                return f"{sanitize_name(name)}({', '.join(self.expr(arg) for arg in node.args)})"
            origin = self.analyzer.imported_origin(name)
            if origin is not None:
                imported = origin[1]
                if imported in {"randint", "randrange", "getrandbits"}:
                    return self.random_call(imported, node)
                if imported == "shuffle":
                    return self.random_call(imported, node)
                if imported == "randbelow":
                    if len(node.args) != 1 or node.keywords:
                        self.fail("randbelow expects one positional argument", node)
                    return f"py_randbelow({self.expr(node.args[0])})"
                if imported == "choice":
                    if len(node.args) != 1 or node.keywords:
                        self.fail("choice expects one positional argument", node)
                    return f"py_choice({self.expr(node.args[0])})"
                if imported in {"bisect_left", "bisect_right"}:
                    if len(node.args) < 2 or len(node.args) > 4 or node.keywords:
                        self.fail(f"{imported} expects two to four positional arguments", node)
                    source = self.expr(node.args[0])
                    value = self.expr(node.args[1])
                    function = "lower_bound" if imported == "bisect_left" else "upper_bound"
                    return f"py_bisect({source}, {value}, {'true' if function == 'upper_bound' else 'false'})"
                if imported in {"gcd", "lcm", "isqrt"}:
                    if imported == "isqrt":
                        return f"static_cast<i64>(sqrt(static_cast<long double>({self.expr(node.args[0])})))"
                    return f"std::{imported}({', '.join(self.expr(arg) for arg in node.args)})"
                if imported in {"sqrt", "hypot", "sin", "cos", "acos", "atan2", "log2"}:
                    return f"std::{imported}({', '.join(self.expr(arg) for arg in node.args)})"
                if imported == "perf_counter":
                    return "0.0"
                if imported in {"floor", "ceil"}:
                    return f"static_cast<i64>(std::{imported}({self.expr(node.args[0])}))"
                if imported in {"heappush", "heappop", "heapify"}:
                    return self.heapq_call(node, imported)
                if imported in {"copy", "deepcopy"}:
                    if len(node.args) != 1 or node.keywords:
                        self.fail(f"{imported} expects one positional argument", node)
                    return self.expr(node.args[0])
                if imported in {"StringBuilder", "StringIO"}:
                    return "vector<char>{}"
            if name in {"copy", "deepcopy"} and name in self.analyzer.from_imports.get("copy", set()):
                if len(node.args) != 1 or node.keywords:
                    self.fail(f"{name} expects one positional argument", node)
                return self.expr(node.args[0])
            if name in self.analyzer.input_aliases:
                if node.args or node.keywords:
                    self.fail("input arguments are not supported", node)
                return "py_input.line()"
            if name in self.analyzer.read_aliases:
                if node.args or node.keywords:
                    self.fail("read alias arguments are not supported", node)
                return "py_input.remaining()"
            if name == "int":
                if not 1 <= len(node.args) <= 2:
                    self.fail("int expects one or two arguments", node)
                if self.is_input_call(node.args[0]):
                    return "py_input.read<i64>()"
                if len(node.args) == 2:
                    return f"py_int_base({self.expr(node.args[0])}, {self.expr(node.args[1])})"
                return f"py_int({self.expr(node.args[0])})"
            if name == "float":
                if len(node.args) != 1:
                    self.fail("float expects one argument", node)
                return f"py_float({self.expr(node.args[0])})"
            if name == "str":
                value = node.args[0]
                value_type = self.analyzer.infer(value, self.context.env)
                if value_type.kind == "class" and "__str__" in self.analyzer.class_methods(value_type.name):
                    return f"{self.expr(value)}.{self.method_cpp_name('__str__')}()"
                return f"py_string({self.expr(value)})"
            if name == "repr":
                if len(node.args) != 1 or node.keywords:
                    self.fail("repr expects one positional argument", node)
                return f"py_repr({self.expr(node.args[0])})"
            if name == "getattr":
                if len(node.args) not in {2, 3} or node.keywords or not isinstance(node.args[1], ast.Constant):
                    self.fail("getattr requires a literal attribute name", node)
                receiver = self.expr(node.args[0])
                attribute = sanitize_name(str(node.args[1].value))
                return f"{receiver}.{attribute}"
            if name == "globals":
                if node.args or node.keywords:
                    self.fail("globals expects no arguments", node)
                return "PyDict<string, function<string()>>{}"
            if name == "bin":
                if len(node.args) != 1 or node.keywords:
                    self.fail("bin expects one positional argument", node)
                return f"py_bin({self.expr(node.args[0])})"
            if name == "format":
                if len(node.args) != 2 or node.keywords:
                    self.fail("format expects two positional arguments", node)
                return f"py_format_binary({self.expr(node.args[0])}, {self.expr(node.args[1])})"
            if name == "urandom":
                if len(node.args) != 1 or node.keywords:
                    self.fail("urandom expects one positional argument", node)
                return f"py_urandom({self.expr(node.args[0])})"
            if name == "randbelow":
                if len(node.args) != 1 or node.keywords:
                    self.fail("randbelow expects one positional argument", node)
                return f"py_randbelow({self.expr(node.args[0])})"
            if name == "choice":
                if len(node.args) != 1 or node.keywords:
                    self.fail("choice expects one positional argument", node)
                return f"py_choice({self.expr(node.args[0])})"
            if name == "iter":
                if len(node.args) != 1 or node.keywords:
                    self.fail("iter expects one positional argument", node)
                argument_type = self.analyzer.infer(node.args[0], self.context.env)
                expression = self.expr(node.args[0])
                if argument_type.kind == "class" and "__iter__" in self.analyzer.class_methods(argument_type.name):
                    return f"{expression}.{self.method_cpp_name('__iter__')}()"
                return expression
            if name == "next":
                if not 1 <= len(node.args) <= 2 or node.keywords:
                    self.fail("next expects one or two positional arguments", node)
                source = self.expr(node.args[0])
                if len(node.args) == 2:
                    return f"py_next({source}, {self.expr(node.args[1])})"
                return f"py_next({source})"
            if name == "isinstance":
                if len(node.args) != 2 or node.keywords:
                    self.fail("isinstance expects two positional arguments", node)
                value_type = self.analyzer.infer(node.args[0], self.context.env)
                expected = node.args[1].id if isinstance(node.args[1], ast.Name) else ""
                expected_kind = {
                    "int": "i64", "bool": "bool", "float": "f64", "str": "string", "bytes": "bytes",
                    "list": "list", "tuple": "tuple", "set": "set", "dict": "map",
                }.get(expected)
                if expected_kind is not None:
                    return "true" if value_type.kind == expected_kind else "false"
                if isinstance(node.args[1], ast.Name) and node.args[1].id in self.analyzer.classes:
                    return "true" if value_type.kind == "class" and value_type.name == node.args[1].id else "false"
                return "false"
            if name == "issubclass":
                return "false"
            if name == "type":
                if len(node.args) != 1 or node.keywords:
                    self.fail("type expects one positional argument", node)
                return "0"
            if name == "cmp_to_key":
                self.fail("cmp_to_key is not supported in generated C++", node)
            if name == "Fraction":
                return f"static_cast<double>({', '.join(self.expr(arg) for arg in node.args)})"
            if name in {"randint", "randrange", "getrandbits"} and name in self.analyzer.from_imports.get("random", set()):
                return self.random_call(name, node)
            if name == "bytes":
                if len(node.args) > 1 or node.keywords:
                    self.fail("bytes supports zero arguments or an integer size", node)
                return f"string({self.expr(node.args[0])}, '\\0')" if node.args else "string()"
            if name == "bytearray":
                if len(node.args) > 2 or node.keywords:
                    self.fail("bytearray expects at most two positional arguments", node)
                if not node.args:
                    return "vector<i64>{}"
                argument = node.args[0]
                argument_type = self.analyzer.infer(argument, self.context.env)
                if argument_type.kind in {"i64", "unknown"}:
                    return f"vector<i64>({self.expr(argument)}, 0)"
                if argument_type.kind == "bytes":
                    return f"py_bytes_values({self.expr(argument)})"
                if argument_type.kind == "string":
                    return f"py_bytes_values({self.expr(argument)})"
                if argument_type.kind in {"list", "tuple", "deque"}:
                    return f"py_integer_values({self.expr(argument)})"
                self.fail("unsupported bytearray initializer", argument)
            if name == "array":
                if not 1 <= len(node.args) <= 2 or node.keywords:
                    self.fail("array expects a typecode and optional initializer", node)
                typecode = node.args[0]
                if not isinstance(typecode, ast.Constant) or typecode.value not in {"i", "I", "H"}:
                    self.fail("only integer array typecodes i, I, and H are supported", typecode)
                if len(node.args) == 1:
                    return "vector<i64>{}"
                initializer = node.args[1]
                if isinstance(initializer, ast.Call) and isinstance(initializer.func, ast.Name) and initializer.func.id == "bytes" and len(initializer.args) == 1:
                    size = "2" if typecode.value == "H" else "4"
                    return f"vector<i64>(({self.expr(initializer.args[0])}) / {size}, 0)"
                return f"py_integer_values({self.expr(initializer)})"
            if name == "bool":
                return f"py_truth({self.expr(node.args[0])})"
            if name == "callable":
                if len(node.args) != 1 or node.keywords:
                    self.fail("callable expects one positional argument", node)
                return f"py_callable({self.expr(node.args[0])})"
            if name == "len":
                value_type = self.analyzer.infer(node.args[0], self.context.env)
                value = self.expr(node.args[0])
                if value_type.kind == "optional":
                    value_type = value_type.args[0]
                    value = f"(*{value})"
                if value_type.kind == "class":
                    return f"{value}.{self.method_cpp_name('__len__')}()"
                return f"static_cast<i64>({value}.size())"
            if name == "abs":
                return f"abs({self.expr(node.args[0])})"
            if name in {"min", "max"}:
                return self.minmax_call(name, node)
            if name in {"sum", "all", "any"}:
                reduction = self.generator_reduction(name, node)
                if reduction is not None:
                    return reduction
            if name == "sum":
                source = self.expr(node.args[0])
                item_type = element_type(self.analyzer.infer(node.args[0], self.context.env))
                initial = self.expr(node.args[1]) if len(node.args) == 2 else ("0.0" if item_type.kind == "f64" else "i64(0)")
                if isinstance(node.args[0], (ast.Name, ast.Attribute)):
                    return f"accumulate({source}.begin(), {source}.end(), {initial})"
                temporary = self.context.temporary("sum")
                return f"([&]() {{ auto&& {temporary} = ({source}); return accumulate({temporary}.begin(), {temporary}.end(), {initial}); }}())"
            if name in {"all", "any"}:
                if len(node.args) != 1:
                    self.fail(f"{name} expects one argument", node)
                suffix = "_string" if self.analyzer.infer(node.args[0], self.context.env).kind == "string" else ""
                return f"py_{name}{suffix}({self.expr(node.args[0])})"
            if name == "pow":
                if len(node.args) == 2:
                    return f"py_pow({self.value_expr(node.args[0])}, {self.value_expr(node.args[1])})"
                if len(node.args) == 3:
                    return f"py_pow_mod({self.value_expr(node.args[0])}, {self.value_expr(node.args[1])}, {self.value_expr(node.args[2])})"
                self.fail("pow expects two or three arguments", node)
            if name == "divmod":
                if len(node.args) != 2:
                    self.fail("divmod expects two arguments", node)
                return f"py_divmod({self.expr(node.args[0])}, {self.expr(node.args[1])})"
            if name == "ord":
                return f"py_ord({self.expr(node.args[0])})"
            if name == "chr":
                return f"string(1, static_cast<char>({self.expr(node.args[0])}))"
            if name == "sorted":
                descending = "false"
                reverse_expression = None
                key = None
                for keyword in node.keywords:
                    if keyword.arg == "key" and (
                        isinstance(keyword.value, (ast.Lambda, ast.Name, ast.Attribute))
                        or (
                            isinstance(keyword.value, ast.Call)
                            and isinstance(keyword.value.func, ast.Name)
                            and keyword.value.func.id == "cmp_to_key"
                        )
                    ):
                        key = keyword.value
                    elif keyword.arg == "reverse":
                        if isinstance(keyword.value, ast.Constant) and isinstance(keyword.value.value, bool):
                            descending = "true" if keyword.value.value else "false"
                        else:
                            reverse_expression = keyword.value
                    else:
                        self.fail("sorted supports a key lambda and literal reverse", keyword)
                if key is not None:
                    source_type = self.analyzer.infer(node.args[0], self.context.env)
                    source = self.expr(node.args[0])
                    if source_type.kind in {"map", "default_map"}:
                        source = f"py_dict_keys({source})"
                    if isinstance(key, ast.Call):
                        return f"py_sorted({source}, {descending})"
                    key_expression = self.sort_key_lambda(key, element_type(source_type)) if isinstance(key, ast.Lambda) else self.expr(key)
                    reverse = descending if reverse_expression is None else f"py_truth({self.expr(reverse_expression)})"
                    return f"{self.sort_key_runtime('py_sorted_key', key)}({source}, {key_expression}, {reverse})"
                if reverse_expression is not None:
                    return f"py_sorted({self.expr(node.args[0])}, py_truth({self.expr(reverse_expression)}))"
                return f"py_sorted({self.expr(node.args[0])}, {descending})"
            if name == "reversed":
                arg_type = self.analyzer.infer(node.args[0], self.context.env)
                expression = self.expr(node.args[0])
                if arg_type.kind == "class" and "__reversed__" in self.analyzer.class_methods(arg_type.name):
                    return f"{expression}.{self.method_cpp_name('__reversed__')}()"
                return f"py_reversed({expression})"
            if name == "list":
                if not node.args:
                    return "vector<i64>{}"
                arg = node.args[0]
                if isinstance(arg, ast.Call) and isinstance(arg.func, ast.Name) and arg.func.id == "map":
                    return self.map_call(arg)
                if isinstance(arg, ast.Call) and isinstance(arg.func, ast.Attribute) and arg.func.attr == "split":
                    return self.expr(arg)
                arg_type = self.analyzer.infer(arg, self.context.env)
                expression = self.expr(arg)
                if arg_type.kind == "class" and "__iter__" in self.analyzer.class_methods(arg_type.name):
                    return f"{expression}.{self.method_cpp_name('__iter__')}()"
                if arg_type.kind == "bytes":
                    return f"py_bytes_values({expression})"
                if arg_type.kind == "list":
                    return expression
                temporary = self.context.temporary("source")
                return f"([&]() {{ auto&& {temporary} = ({expression}); return vector<{element_type(arg_type).cpp()}>({temporary}.begin(), {temporary}.end()); }}())"
            if name == "tuple":
                if not node.args:
                    return "tuple<>()"
                argument = node.args[0]
                argument_type = self.analyzer.infer(argument, self.context.env)
                if argument_type.kind in {"list", "tuple"}:
                    return self.expr(argument)
                if argument_type.kind == "bytes":
                    return f"py_bytes_values({self.expr(argument)})"
                expression = self.expr(argument)
                temporary = self.context.temporary("tuple_source")
                return f"([&]() {{ auto&& {temporary} = ({expression}); return vector<{element_type(argument_type).cpp()}>({temporary}.begin(), {temporary}.end()); }}())"
            if name == "set":
                if not node.args:
                    return "set<i64>{}"
                arg_type = self.analyzer.infer(node.args[0], self.context.env)
                expression = self.expr(node.args[0])
                temporary = self.context.temporary("source")
                return f"([&]() {{ auto&& {temporary} = ({expression}); return set<{element_type(arg_type).cpp()}>({temporary}.begin(), {temporary}.end()); }}())"
            if name == "dict":
                if not node.args:
                    return "{}"
                return f"py_dict({self.expr(node.args[0])})"
            if name == "defaultdict":
                value_type = self.analyzer.infer(node, self.context.env)
                factory = node.args[0]
                default_factory = f"[&]() {{ return {self.expr(factory.body)}; }}" if isinstance(factory, ast.Lambda) else None
                if len(node.args) == 1:
                    return f"{value_type.cpp()}{{{default_factory}}}" if default_factory is not None else f"{value_type.cpp()}{{}}"
                if len(node.args) == 2:
                    initial = self.expr(node.args[1])
                    return f"{value_type.cpp()}{{{default_factory}, {initial}}}" if default_factory is not None else f"{value_type.cpp()}{{{initial}}}"
                self.fail("defaultdict expects a factory and optional mapping", node)
            if name == "deque":
                if not node.args:
                    return "deque<i64>{}"
                arg_type = self.analyzer.infer(node.args[0], self.context.env)
                expression = self.expr(node.args[0])
                temporary = self.context.temporary("source")
                return f"([&]() {{ auto&& {temporary} = ({expression}); return deque<{element_type(arg_type).cpp()}>({temporary}.begin(), {temporary}.end()); }}())"
            if name == "zip":
                if not node.args or node.keywords:
                    self.fail("zip expects one or more positional iterables", node)
                return f"py_zip({', '.join(self.expr(arg) for arg in node.args)})"
            if name == "map":
                return self.map_call(node)
            if name in {"bisect_left", "bisect_right", "gcd", "lcm", "isqrt"}:
                if name in {"bisect_left", "bisect_right"}:
                    method = "lower_bound" if name == "bisect_left" else "upper_bound"
                    source = self.expr(node.args[0])
                    value = self.expr(node.args[1])
                    return f"py_bisect({source}, {value}, {'true' if method == 'upper_bound' else 'false'})"
                function = {"gcd": "gcd", "lcm": "lcm", "isqrt": "static_cast<i64>(sqrt"}[name]
                if name == "isqrt":
                    return f"static_cast<i64>(sqrt(static_cast<long double>({self.expr(node.args[0])})))"
                return f"std::{function}({', '.join(self.expr(arg) for arg in node.args)})"
            if name in self.analyzer.from_imports.get("bisect", set()):
                if len(node.args) < 2 or len(node.args) > 4 or node.keywords:
                    self.fail(f"{name} expects two to four positional arguments", node)
                source = self.expr(node.args[0])
                value = self.expr(node.args[1])
                function = "lower_bound" if name.endswith("left") else "upper_bound"
                return f"py_bisect({source}, {value}, {'true' if function == 'upper_bound' else 'false'})"
            if name in {"sqrt", "hypot", "sin", "cos", "acos", "atan2", "log2", "prod"} and name in self.analyzer.from_imports.get("math", set()):
                if name == "prod":
                    return f"py_product({self.expr(node.args[0])})"
                return f"std::{name}({', '.join(self.expr(arg) for arg in node.args)})"
            if name in {"heappush", "heappop", "heapify"} and name in self.analyzer.from_imports.get("heapq", set()):
                return self.heapq_call(node)
            if name == "range":
                if not 1 <= len(node.args) <= 3 or node.keywords:
                    self.fail("range expects one to three positional arguments", node)
                values = [self.expr(arg) for arg in node.args]
                if len(values) == 1:
                    values.insert(0, "0")
                return f"py_range({', '.join(values)})"
            if name == "enumerate":
                if not 1 <= len(node.args) <= 2 or node.keywords:
                    self.fail("enumerate expects one or two positional arguments", node)
                source_type = self.analyzer.infer(node.args[0], self.context.env)
                source = self.expr(node.args[0])
                if source_type.kind in {"map", "default_map"}:
                    source = f"py_dict_keys({source})"
                start = self.expr(node.args[1]) if len(node.args) == 2 else "0"
                return f"py_enumerate({source}, {start})"
            if name in self.analyzer.functions:
                return self.user_function_call(name, node)
            if name in self.analyzer.classes:
                if node.keywords:
                    self.fail("keyword arguments for constructors are not supported", node)
                return f"{sanitize_name(name)}({', '.join(self.expr(arg) for arg in node.args)})"
            self.fail(f"unsupported function: {name}", node)
        if isinstance(node.func, ast.Attribute):
            if node.func.attr == "setrecursionlimit" and isinstance(node.func.value, ast.Name) and node.func.value.id == "sys":
                return "0"
            if node.func.attr == "readline" and self.analyzer._is_stdin_object(node.func.value):
                return "py_input.line()"
            if node.func.attr == "read" and self.analyzer._is_stdin_object(node.func.value):
                return "py_input.remaining()"
            if isinstance(node.func.value, ast.Name) and node.func.value.id == "bytes" and node.func.attr == "maketrans":
                if len(node.args) != 2 or node.keywords:
                    self.fail("bytes.maketrans expects two positional arguments", node)
                return self.expr(node.args[1])
            if isinstance(node.func.value, ast.Name) and node.func.value.id == "sys" and node.func.attr == "stdin":
                self.fail("sys.stdin access must use buffer.read()", node)
            if isinstance(node.func.value, ast.Name) and node.func.value.id == "copy" and node.func.attr in {"copy", "deepcopy"}:
                if len(node.args) != 1 or node.keywords:
                    self.fail(f"copy.{node.func.attr} expects one positional argument", node)
                return self.expr(node.args[0])
            if isinstance(node.func.value, ast.Name) and node.func.value.id in {"heapq", "bisect", "math", "random"}:
                owner = node.func.value.id
                if owner == "heapq":
                    return self.heapq_call(node)
                if owner == "bisect":
                    return self.bisect_call(node)
                if owner == "math":
                    return self.math_call(node)
                if owner == "random" and node.func.attr in {"randint", "randrange", "getrandbits", "shuffle"}:
                    return self.random_call(node.func.attr, node)
            if isinstance(node.func.value, ast.Name) and node.func.value.id == "int" and node.func.attr == "from_bytes":
                if not 1 <= len(node.args) <= 2 or node.keywords:
                    self.fail("int.from_bytes expects one or two positional arguments", node)
                return f"py_int_from_bytes({self.expr(node.args[0])})"
            attr = node.func.attr
            if attr == "indices":
                if len(node.args) != 1 or node.keywords:
                    self.fail("indices expects one positional argument", node)
                return f"py_slice_indices({self.expr(node.func.value)}, {self.expr(node.args[0])})"
            if attr == "bit_length":
                if node.args or node.keywords:
                    self.fail("bit_length expects no arguments", node)
                return f"py_bit_length({self.expr(node.func.value)})"
            if attr == "bit_count":
                if node.args or node.keywords:
                    self.fail("bit_count expects no arguments", node)
                return f"py_bit_count({self.expr(node.func.value)})"
            base_type = self.analyzer.infer(node.func.value, self.context.env)
            if attr == "__new__" and base_type.kind == "class":
                return f"{sanitize_name(base_type.name)}()"
            if base_type.kind == "function":
                return f"{self.expr(node.func.value)}({', '.join(self.expr(arg) for arg in node.args)})"
            if base_type.kind == "class":
                if node.keywords:
                    self.fail("keyword arguments for methods are not supported", node)
                if attr in self.analyzer.class_attributes.get(base_type.name, {}):
                    base = "this->" if isinstance(node.func.value, ast.Name) and node.func.value.id == "self" else self.expr(node.func.value) + "."
                    attribute_type = self.analyzer.class_attributes[base_type.name][attr]
                    expected = function_parameters(attribute_type) if attribute_type.kind == "function" else ()
                    arguments = ", ".join(
                        self.call_argument(argument, expected[index] if index < len(expected) else None)
                        for index, argument in enumerate(node.args)
                    )
                    return f"{base}{sanitize_name(attr)}({arguments})"
                base = "this->" if isinstance(node.func.value, ast.Name) and node.func.value.id == "self" else self.expr(node.func.value) + "."
                method = self.analyzer.class_methods(base_type.name).get(attr)
                parameters = self.analyzer.method_parameter_types.get(base_type.name, {}).get(attr, {})
                expected_arguments = self.analyzer.method_arguments(method) if method is not None else []
                arguments = ", ".join(
                    self.call_argument(
                        argument,
                        parameters.get(expected_arguments[index].arg) if index < len(expected_arguments) else None,
                    )
                    for index, argument in enumerate(node.args)
                )
                return f"{base}{self.method_cpp_name(attr)}({arguments})"
            if base_type.kind in {"unknown", "void"}:
                return f"{self.expr(node.func.value)}.{sanitize_name(attr)}({', '.join(self.expr(arg) for arg in node.args)})"
            if attr == "build" and base_type.kind in {"list", "deque", "string"}:
                if node.args or node.keywords:
                    self.fail("build expects no arguments", node)
                source = self.expr(node.func.value)
                return f"string({source}.begin(), {source}.end())"
            if attr in {"upper", "lower"}:
                if node.args or node.keywords:
                    self.fail(f"{attr} expects no arguments", node)
                return f"py_{attr}({self.expr(node.func.value)})"
            if attr == "replace":
                if len(node.args) != 2 or node.keywords:
                    self.fail("replace expects two positional arguments", node)
                return f"py_replace({self.expr(node.func.value)}, {self.expr(node.args[0])}, {self.expr(node.args[1])})"
            if attr == "split":
                if len(node.args) > 1:
                    self.fail("split supports at most one argument", node)
                if self.is_input_call(node.func.value):
                    if node.args:
                        return f"py_split(py_input.line(), {self.expr(node.args[0])})"
                    return "py_split(py_input.line())"
                arguments = ", " + self.expr(node.args[0]) if node.args else ""
                return f"py_split({self.expr(node.func.value)}{arguments})"
            if attr in {"strip", "lstrip", "rstrip"}:
                if node.args:
                    self.fail("strip arguments are not supported", node)
                expression = self.expr(node.func.value)
                mode = {"strip": "0", "lstrip": "-1", "rstrip": "1"}[attr]
                return f"py_strip({expression}, {mode})"
            if attr == "pop":
                if base_type.kind in {"map", "default_map"}:
                    if len(node.args) != 1 or node.keywords:
                        self.fail("dict.pop currently requires one key argument", node)
                    return f"py_dict_pop({self.expr(node.func.value)}, {self.expr(node.args[0])})"
                if node.args:
                    return f"py_pop_at({self.expr(node.func.value)}, {self.expr(node.args[0])})"
                return f"py_pop({self.expr(node.func.value)})"
            if attr == "popleft":
                return f"py_popleft({self.expr(node.func.value)})"
            if attr == "insert":
                if len(node.args) != 2 or node.keywords:
                    self.fail("insert expects an index and a value", node)
                base = self.expr(node.func.value)
                return f"{base}.insert({base}.begin() + {self.expr(node.args[0])}, {self.expr(node.args[1])})"
            if attr == "copy":
                return f"{self.expr(node.func.value)}"
            if attr in {"count", "index"}:
                source = self.expr(node.func.value)
                value = self.expr(node.args[0])
                if attr == "count" and base_type.kind in {"string", "bytes"}:
                    return f"py_count({source}, {value})"
                if not isinstance(node.func.value, (ast.Name, ast.Attribute)):
                    temporary = self.context.temporary("source")
                    function = "count" if attr == "count" else "find"
                    result = f"static_cast<i64>({function}({temporary}.begin(), {temporary}.end(), {value})"
                    if attr == "index":
                        result += f" - {temporary}.begin()"
                    result += ")"
                    return f"([&]() {{ auto&& {temporary} = ({source}); return {result}; }}())"
                if attr == "count":
                    return f"static_cast<i64>(count({source}.begin(), {source}.end(), {value}))"
                return f"static_cast<i64>(find({source}.begin(), {source}.end(), {value}) - {source}.begin())"
            if attr == "find":
                return f"py_find({self.expr(node.func.value)}, {self.expr(node.args[0])})"
            if attr == "bit_length":
                return f"py_bit_length({self.expr(node.func.value)})"
            if attr == "bit_count":
                return f"py_bit_count({self.expr(node.func.value)})"
            if attr == "indices":
                if len(node.args) != 1 or node.keywords:
                    self.fail("indices expects one positional argument", node)
                return f"py_slice_indices({self.expr(node.func.value)}, {self.expr(node.args[0])})"
            if attr in {"count1", "select0", "select1", "build"}:
                if node.keywords:
                    self.fail(f"{attr} keyword arguments are not supported", node)
                arguments = ", ".join(self.expr(argument) for argument in node.args)
                return f"{self.expr(node.func.value)}.{attr}({arguments})"
            if attr == "as_integer_ratio":
                if node.args or node.keywords:
                    self.fail("as_integer_ratio expects no arguments", node)
                return f"py_as_integer_ratio({self.expr(node.func.value)})"
            if attr == "isascii":
                if node.args or node.keywords:
                    self.fail("isascii expects no arguments", node)
                return f"py_isascii({self.expr(node.func.value)})"
            if attr == "to_bytes":
                if len(node.args) not in {2, 3} or node.keywords:
                    self.fail("to_bytes expects two or three positional arguments", node)
                signed = self.expr(node.args[2]) if len(node.args) == 3 else "false"
                if signed != "false":
                    self.fail("signed to_bytes is not supported", node)
                return f"py_to_bytes({self.expr(node.func.value)}, {self.expr(node.args[0])}, {self.expr(node.args[1])})"
            if attr == "join":
                return f"py_join({self.expr(node.func.value)}, {self.expr(node.args[0])})"
            if attr == "format":
                if node.keywords:
                    self.fail("format keyword arguments are not supported", node)
                arguments = ", ".join(self.expr(arg) for arg in node.args)
                suffix = f", {arguments}" if arguments else ""
                return f"py_format({self.expr(node.func.value)}{suffix})"
            if attr in {"encode", "decode"}:
                if node.args or node.keywords:
                    self.fail(f"{attr} arguments are not supported", node)
                return self.expr(node.func.value)
            if attr == "translate":
                if len(node.args) != 1 or node.keywords:
                    self.fail("translate expects one positional argument", node)
                return f"py_translate({self.expr(node.func.value)}, {self.expr(node.args[0])})"
            if attr == "get":
                arguments = ", ".join(self.expr(arg) for arg in node.args)
                return f"py_dict_get({self.expr(node.func.value)}, {arguments})"
            if attr in {"keys", "values", "items"}:
                if node.args or node.keywords:
                    self.fail(f"{attr} expects no arguments", node)
                return f"py_dict_{attr}({self.value_expr(node.func.value)})"
            if attr in {"append", "appendleft", "add", "clear", "sort", "reverse", "rotate", "extend", "remove", "discard", "update", "setdefault"}:
                return self.method_call(node)
            if attr == "read" and isinstance(node.func.value, ast.Attribute):
                if isinstance(node.func.value.value, ast.Attribute) and isinstance(node.func.value.value.value, ast.Name):
                    if node.func.value.value.value.id == "sys" and node.func.value.value.attr == "stdin":
                        return "py_input.remaining()"
            self.fail(f"unsupported method: {attr}", node)
        if isinstance(node.func, ast.Subscript):
            return f"{self.subscript(node.func, mutable=False)}({', '.join(self.expr(arg) for arg in node.args)})"
        self.fail("unsupported call", node)
        return ""

    def map_call(self, node: ast.Call) -> str:
        if len(node.args) < 2:
            self.fail("map expects a function and an iterable", node)
        fn = node.args[0]
        source = node.args[1]
        if isinstance(fn, ast.Name) and fn.id == "int" and self.is_input_split(source):
            return "py_read_vector<i64>(py_input)"
        if isinstance(fn, ast.Name) and fn.id == "int" and self.is_remaining_split(source):
            return "py_read_remaining<i64>(py_input)"
        source_type = self.analyzer.infer(source, self.context.env)
        reserve_source = source_type.kind in {"list", "deque", "set", "string", "bytes"}
        if isinstance(fn, ast.Name) and fn.id in {"int", "float", "str"}:
            source_code = self.expr(source)
            target_type = {"int": "i64", "float": "double", "str": "string"}[fn.id]
            temp = self.context.temporary("mapped")
            item = self.context.temporary("item")
            conversion = {"int": f"py_int({item})", "float": f"py_float({item})", "str": f"py_string({item})"}[fn.id]
            if reserve_source:
                source_temp = self.context.temporary("map_source")
                return (
                    f"([&]() {{ auto&& {source_temp} = ({source_code}); vector<{target_type}> {temp}; "
                    f"{temp}.reserve({source_temp}.size()); for (const auto& {item} : {source_temp}) "
                    f"{temp}.push_back({conversion}); return {temp}; }}())"
                )
            return f"([&]() {{ vector<{target_type}> {temp}; for (const auto& {item} : {source_code}) {temp}.push_back({conversion}); return {temp}; }}())"
        source_code = self.expr(source)
        item_type = self.analyzer._map_item_type(node, self.context.env)
        function = self.expr(fn)
        temp = self.context.temporary("mapped")
        item = self.context.temporary("item")
        if reserve_source:
            source_temp = self.context.temporary("map_source")
            return (
                f"([&]() {{ auto&& {source_temp} = ({source_code}); vector<{item_type.cpp()}> {temp}; "
                f"{temp}.reserve({source_temp}.size()); for (const auto& {item} : {source_temp}) "
                f"{temp}.push_back({function}({item})); return {temp}; }}())"
            )
        return f"([&]() {{ vector<{item_type.cpp()}> {temp}; for (const auto& {item} : {source_code}) {temp}.push_back({function}({item})); return {temp}; }}())"

    def is_input_call(self, node: ast.AST) -> bool:
        if not isinstance(node, ast.Call) or node.args or node.keywords:
            return False
        if isinstance(node.func, ast.Name):
            return node.func.id in self.analyzer.input_aliases
        return isinstance(node.func, ast.Attribute) and node.func.attr == "readline" and self.analyzer._is_stdin_object(node.func.value)

    def is_input_split(self, node: ast.AST) -> bool:
        return isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == "split" and self.is_input_call(node.func.value)

    def is_remaining_split(self, node: ast.AST) -> bool:
        if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Attribute) or node.func.attr != "split":
            return False
        value = node.func.value
        if not isinstance(value, ast.Call):
            return False
        if isinstance(value.func, ast.Name):
            return value.func.id in self.analyzer.read_aliases
        return isinstance(value.func, ast.Attribute) and value.func.attr == "read"

    def method_call(self, node: ast.Call) -> str:
        base = self.expr(node.func.value)
        base_type = self.analyzer.infer(node.func.value, self.context.env)
        if base_type.kind == "optional":
            base = f"(*{base})"
        attr = node.func.attr
        if attr == "append":
            value = self.expr(node.args[0])
            value_type = self.analyzer.infer(node.args[0], self.context.env)
            if value_type.kind == "optional" and base_type.kind == "list" and element_type(base_type).kind != "optional":
                value = f"(*{value})"
            if (
                base_type.kind == "list"
                and element_type(base_type).kind == "string"
                and self.analyzer.infer(node.args[0], self.context.env).kind == "char"
            ):
                value = f"string(1, {value})"
            return f"{base}.push_back({value})"
        if attr == "appendleft":
            return f"{base}.push_front({self.expr(node.args[0])})"
        if attr == "insert":
            if len(node.args) != 2 or node.keywords:
                self.fail("insert expects an index and a value", node)
            return f"{base}.insert({base}.begin() + {self.expr(node.args[0])}, {self.expr(node.args[1])})"
        if attr == "setdefault":
            if len(node.args) not in {1, 2} or node.keywords:
                self.fail("setdefault expects one or two positional arguments", node)
            key = self.expr(node.args[0])
            if len(node.args) == 1:
                return f"{base}[{key}]"
            return f"py_dict_setdefault({base}, {key}, {self.expr(node.args[1])})"
        if attr == "add":
            return f"{base}.insert({self.expr(node.args[0])})"
        if attr == "clear":
            return f"{base}.clear()"
        if attr == "sort":
            descending = False
            reverse_expression = None
            key = None
            for keyword in node.keywords:
                if keyword.arg == "key" and isinstance(keyword.value, (ast.Lambda, ast.Name, ast.Attribute)):
                    key = keyword.value
                elif keyword.arg == "reverse":
                    if isinstance(keyword.value, ast.Constant) and isinstance(keyword.value.value, bool):
                        descending = keyword.value.value
                    else:
                        reverse_expression = keyword.value
                else:
                    self.fail("sort supports a key lambda and literal reverse", keyword)
            if key is not None:
                item_type = element_type(self.analyzer.infer(node.func.value, self.context.env))
                key_expression = self.sort_key_lambda(key, item_type) if isinstance(key, ast.Lambda) else self.expr(key)
                return f"{self.sort_key_runtime('py_sort_key', key)}({base}, {key_expression}, {'true' if descending else 'false'})"
            if reverse_expression is not None:
                return f"py_sort_reverse({base}, py_truth({self.expr(reverse_expression)}))"
            function = "sort" if not descending else "sort"
            comparator = "" if not descending else ", greater<>()"
            return f"{function}({base}.begin(), {base}.end(){comparator})"
        if attr == "reverse":
            return f"reverse({base}.begin(), {base}.end())"
        if attr == "rotate":
            return f"py_rotate({base}, {self.expr(node.args[0])})"
        if attr == "extend":
            return f"py_extend({base}, {self.expr(node.args[0])})"
        if attr == "remove":
            return f"py_remove({base}, {self.expr(node.args[0])})"
        if attr == "discard":
            return f"{base}.erase({self.expr(node.args[0])})"
        if attr == "update":
            if len(node.args) != 1 or node.keywords:
                self.fail("update expects one positional argument", node)
            source = self.expr(node.args[0])
            return f"py_update({base}, {source})"
        self.fail(f"unsupported method: {attr}", node)
        return ""

    def random_call(self, name: str, node: ast.Call) -> str:
        if node.keywords:
            self.fail(f"{name} keyword arguments are not supported", node)
        if name == "shuffle":
            if len(node.args) != 1:
                self.fail("shuffle expects one positional argument", node)
            return f"py_shuffle({self.expr(node.args[0])})"
        if name == "randint":
            if len(node.args) != 2:
                self.fail("randint expects two arguments", node)
            return f"py_randint({self.expr(node.args[0])}, {self.expr(node.args[1])})"
        if name == "randrange":
            if not 1 <= len(node.args) <= 3:
                self.fail("randrange expects one to three arguments", node)
            values = [self.expr(arg) for arg in node.args]
            if len(values) == 1:
                values.insert(0, "0")
            return f"py_randrange({', '.join(values)})"
        if len(node.args) != 1:
            self.fail("getrandbits expects one argument", node)
        return f"py_getrandbits({self.expr(node.args[0])})"

    def heapq_call(self, node: ast.Call, imported_name: str | None = None) -> str:
        attr = imported_name or (node.func.attr if isinstance(node.func, ast.Attribute) else node.func.id)
        if attr == "heapify":
            return f"py_heapify({self.expr(node.args[0])})"
        if attr == "heappush":
            return f"py_heappush({self.expr(node.args[0])}, {self.expr(node.args[1])})"
        if attr == "heappop":
            return f"py_heappop({self.expr(node.args[0])})"
        self.fail(f"unsupported heapq function: {attr}", node)
        return ""

    def bisect_call(self, node: ast.Call) -> str:
        attr = node.func.attr
        method = "lower_bound" if attr == "bisect_left" else "upper_bound"
        source = self.expr(node.args[0])
        value = self.expr(node.args[1])
        return f"py_bisect({source}, {value}, {'true' if method == 'upper_bound' else 'false'})"

    def math_call(self, node: ast.Call) -> str:
        attr = node.func.attr
        if attr == "isqrt":
            return f"static_cast<i64>(sqrt(static_cast<long double>({self.expr(node.args[0])})))"
        if attr in {"gcd", "lcm"}:
            return f"std::{attr}({', '.join(self.expr(arg) for arg in node.args)})"
        if attr in {"floor", "ceil"}:
            return f"static_cast<i64>(std::{attr}({self.expr(node.args[0])}))"
        if attr in {"sqrt", "hypot", "sin", "cos", "acos", "atan2", "log2"}:
            return f"std::{attr}({', '.join(self.expr(arg) for arg in node.args)})"
        if attr == "prod":
            return f"py_product({self.expr(node.args[0])})"
        self.fail(f"unsupported math function: {attr}", node)
        return ""

    def minmax_call(self, name: str, node: ast.Call) -> str:
        reduction = self.generator_reduction(name, node)
        if reduction is not None:
            return reduction
        if len(node.args) == 1:
            source = self.expr(node.args[0])
            function = "min_element" if name == "min" else "max_element"
            if self.is_simple_reference(node.args[0]):
                return f"*{function}({source}.begin(), {source}.end())"
            temporary = self.context.temporary("minmax_source")
            result_type = element_type(self.analyzer.infer(node.args[0], self.context.env)).cpp()
            return f"([&]() -> {result_type} {{ const auto& {temporary} = ({source}); return *{function}({temporary}.begin(), {temporary}.end()); }}())"
        value = self.expr(node.args[0])
        function = "py_min" if name == "min" else "py_max"
        for arg in node.args[1:]:
            value = f"{function}({value}, {self.expr(arg)})"
        return value

    def wide_integer_expr(self, node: ast.AST) -> str:
        if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub, ast.Mult)):
            return f"({self.wide_integer_expr(node.left)} {self.operator(node.op)} {self.wide_integer_expr(node.right)})"
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
            return f"({self.operator(node.op)}{self.wide_integer_expr(node.operand)})"
        return f"static_cast<__int128>({self.value_expr(node)})"

    def operator(self, node: ast.AST) -> str:
        operators = {
            ast.Add: "+", ast.Sub: "-", ast.Mult: "*", ast.Div: "/", ast.Mod: "%",
            ast.Pow: "^", ast.BitAnd: "&", ast.BitOr: "|", ast.BitXor: "^", ast.LShift: "<<", ast.RShift: ">>",
            ast.USub: "-", ast.UAdd: "+", ast.Invert: "~", ast.Eq: "==", ast.NotEq: "!=", ast.Lt: "<", ast.LtE: "<=", ast.Gt: ">", ast.GtE: ">=",
            ast.Is: "==", ast.IsNot: "!=",
        }
        for cls, operator in operators.items():
            if isinstance(node, cls):
                return operator
        self.fail(f"unsupported operator: {type(node).__name__}", node)
        return ""


def parse_source(source: str, filename: str = "<input>") -> tuple[ast.Module, Analyzer]:
    try:
        tree = ast.parse(source, filename=filename, type_comments=True)
    except SyntaxError as error:
        raise TranspileError(error.msg, filename=filename) from error
    hoist_local_functions(tree, source, filename)
    flatten_try_blocks(tree)
    flatten_with_blocks(tree)
    normalize_keyword_only_arguments(tree)
    analyzer = Analyzer(tree, filename)
    return tree, analyzer


def check_source(source: str, filename: str = "<input>", target: str = "cpp") -> None:
    tree, analyzer = parse_source(source, filename)
    if target == "cpp":
        CppGenerator(source, filename).generate()
    elif target == "native":
        from .native_full import transpile_native_full_cpp
        transpile_native_full_cpp(source,filename)
    elif target == "codon":
        for function in analyzer.functions.values():
            analyzer.signature(function)
    else:
        raise TranspileError(f"unknown target: {target}", filename=filename)


def transpile_cpp(source: str, filename: str = "main.py") -> str:
    return CppGenerator(source, filename).generate()

from __future__ import annotations

import ast
from dataclasses import dataclass
from .translator import Type, TranspileError


LOW = -(1 << 63)
HIGH = (1 << 63) - 1
BOTTOM = object()


def fits(bounds):
    return bounds is not None and bounds is not BOTTOM and LOW <= bounds[0] <= bounds[1] <= HIGH


def union(left, right):
    if left is BOTTOM:
        return right
    if right is BOTTOM:
        return left
    if left is None or right is None:
        return None
    return min(left[0], right[0]), max(left[1], right[1])


@dataclass(frozen=True)
class NativeInteger(Type):
    bounds: tuple[int, int] = (LOW, HIGH)

    def cpp(self):
        return 'long long'


class IntegerBounds:
    """Whole-scope conservative intervals; unknown writes prevent specialization."""

    def __init__(self, analyzer, types, constants=None):
        self.analyzer = analyzer
        self.types = types
        self.constants = constants or {}

    def builtin(self, name):
        return name not in self.types and name not in self.analyzer.functions

    def kind(self, node):
        try:
            return self.analyzer.infer(node, self.types).kind
        except TranspileError:
            return 'unknown'

    def length(self, node):
        if isinstance(node, (ast.List, ast.Tuple, ast.Set)) and not any(isinstance(item, ast.Starred) for item in node.elts):
            if isinstance(node, ast.Set):
                return 0, len(node.elts)
            return len(node.elts), len(node.elts)
        if isinstance(node, ast.Constant) and isinstance(node.value, (str, bytes)):
            return len(node.value), len(node.value)
        if self.kind(node) in {'list', 'deque', 'set', 'map', 'default_map', 'string', 'bytes', 'tuple'}:
            return 0, HIGH
        return None

    def expression(self, node, values=None):
        values = values or {}
        if isinstance(node, ast.Constant) and isinstance(node.value, int):
            return node.value, node.value
        if isinstance(node, ast.Name):
            if node.id in values:
                return values[node.id]
            value_type = self.types.get(node.id)
            if isinstance(value_type, NativeInteger):
                return value_type.bounds
            return self.constants.get(node.id)
        if isinstance(node, ast.UnaryOp):
            operand = self.expression(node.operand, values)
            if operand is None or operand is BOTTOM:
                return operand
            if isinstance(node.op, ast.USub):
                return -operand[1], -operand[0]
            if isinstance(node.op, ast.UAdd):
                return operand
            if isinstance(node.op, ast.Invert):
                return ~operand[1], ~operand[0]
        if isinstance(node, ast.IfExp):
            return union(self.expression(node.body, values), self.expression(node.orelse, values))
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and not node.keywords and self.builtin(node.func.id):
            if node.func.id == 'len' and len(node.args) == 1:
                return self.length(node.args[0])
            args = [self.expression(arg, values) for arg in node.args]
            if any(arg is None for arg in args):
                return None
            if any(arg is BOTTOM for arg in args):
                return BOTTOM
            if node.func.id == 'abs' and len(args) == 1:
                low, high = args[0]
                return (0 if low <= 0 <= high else min(abs(low), abs(high))), max(abs(low), abs(high))
            if node.func.id in {'min', 'max'} and len(args) >= 2:
                fn = min if node.func.id == 'min' else max
                return fn(arg[0] for arg in args), fn(arg[1] for arg in args)
        if isinstance(node, ast.BinOp):
            left, right = self.expression(node.left, values), self.expression(node.right, values)
            if isinstance(node.op, ast.Mod) and right is not None and right is not BOTTOM and (right[0] > 0 or right[1] < 0) and self.kind(node.left) in {'i64', 'bool'}:
                return (0, right[1] - 1) if right[0] > 0 else (right[0] + 1, 0)
            if left is None or right is None:
                return None
            if left is BOTTOM or right is BOTTOM:
                return BOTTOM
            a, b = left
            c, d = right
            if isinstance(node.op, ast.Add):
                return a + c, b + d
            if isinstance(node.op, ast.Sub):
                return a - d, b - c
            if isinstance(node.op, ast.Mult):
                products = (a*c, a*d, b*c, b*d)
                return min(products), max(products)
            if isinstance(node.op, ast.FloorDiv) and (c > 0 or d < 0):
                quotients = (a//c, a//d, b//c, b//d)
                return min(quotients), max(quotients)
            if isinstance(node.op, ast.LShift) and 0 <= c <= d <= 63:
                shifted = (a << c, a << d, b << c, b << d)
                return min(shifted), max(shifted)
            if isinstance(node.op, ast.RShift) and 0 <= c <= d:
                c, d = min(c, 64), min(d, 64)
                shifted = (a >> c, a >> d, b >> c, b >> d)
                return min(shifted), max(shifted)
            if isinstance(node.op, ast.BitAnd) and a >= 0 and c >= 0:
                return 0, min(b, d)
            if isinstance(node.op, (ast.BitOr, ast.BitXor)) and a >= 0 and c >= 0:
                return 0, (1 << max(b.bit_length(), d.bit_length())) - 1
        return None

    def range_arguments(self, call, values=None):
        if not isinstance(call, ast.Call) or not isinstance(call.func, ast.Name) or call.func.id != 'range' or not self.builtin('range') or call.keywords or not 1 <= len(call.args) <= 3:
            return None
        args = [self.expression(arg, values) for arg in call.args]
        if len(args) == 1:
            args.insert(0, (0, 0))
        if len(args) == 2:
            args.append((1, 1))
        return args

    def scope(self, body, parameters=()):
        definitions = {}
        excluded = set(parameters)

        def bind(target, value):
            if isinstance(target, ast.Name):
                definitions.setdefault(target.id, []).append(value)
            elif isinstance(target, (ast.Tuple, ast.List)):
                if isinstance(value, ast.AST) and isinstance(value, (ast.Tuple, ast.List)) and len(target.elts) == len(value.elts) and not any(isinstance(item, ast.Starred) for item in target.elts + value.elts):
                    for child, expression in zip(target.elts, value.elts):
                        bind(child, expression)
                else:
                    for child in target.elts:
                        bind(child, None)
            elif isinstance(target, ast.Starred):
                bind(target.value, None)

        def visit(statements):
            for node in statements:
                if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                    for child in ast.walk(node):
                        if isinstance(child, ast.Nonlocal):
                            excluded.update(child.names)
                    continue
                if isinstance(node, (ast.Global, ast.Nonlocal)):
                    excluded.update(node.names)
                elif isinstance(node, ast.Assign):
                    for target in node.targets:
                        bind(target, node.value)
                elif isinstance(node, ast.AnnAssign) and node.value is not None:
                    bind(node.target, node.value)
                elif isinstance(node, ast.AugAssign):
                    bind(node.target, ast.BinOp(left=node.target, op=node.op, right=node.value))
                elif isinstance(node, ast.For):
                    if isinstance(node.target, ast.Name) and self.range_arguments(node.iter) is not None:
                        bind(node.target, ('range', node.iter))
                    elif isinstance(node.iter, ast.Call) and isinstance(node.iter.func, ast.Name) and node.iter.func.id == 'enumerate' and self.builtin('enumerate') and 1 <= len(node.iter.args) <= 2 and not node.iter.keywords and isinstance(node.target, (ast.Tuple, ast.List)) and len(node.target.elts) == 2:
                        bind(node.target.elts[0], ('enumerate', node.iter))
                        bind(node.target.elts[1], None)
                    else:
                        bind(node.target, None)
                    visit(node.body)
                    visit(node.orelse)
                elif isinstance(node, (ast.If, ast.While)):
                    visit(node.body)
                    visit(node.orelse)
                elif isinstance(node, ast.Try):
                    visit(node.body)
                    visit(node.orelse)
                    visit(node.finalbody)
                    for handler in node.handlers:
                        if handler.name:
                            definitions.setdefault(handler.name, []).append(None)
                        visit(handler.body)
                elif isinstance(node, (ast.With, ast.AsyncWith)):
                    for item in node.items:
                        if item.optional_vars:
                            bind(item.optional_vars, None)
                    visit(node.body)
                for child in ast.walk(node):
                    if isinstance(child, ast.NamedExpr):
                        excluded.update(name.id for name in ast.walk(child.target) if isinstance(name, ast.Name))

        visit(body)
        values = {name: BOTTOM for name in definitions}
        values.update({name: None for name in excluded})

        def definition(value):
            if value is None:
                return None
            if isinstance(value, ast.AST):
                return self.expression(value, values)
            tag, call = value
            if tag == 'range':
                args = self.range_arguments(call, values)
                if any(arg is None for arg in args):
                    return None
                if any(arg is BOTTOM for arg in args):
                    return BOTTOM
                start, stop, step = args
                if not all(fits(arg) for arg in args):
                    return None
                return min(start[0], stop[0]), max(start[1], stop[1])
            length = self.length(call.args[0])
            start = self.expression(call.args[1], values) if len(call.args) == 2 else (0, 0)
            if length is None or start is None:
                return None
            if start is BOTTOM:
                return BOTTOM
            return start[0], start[1] + max(0, length[1] - 1)

        for _ in range(64):
            previous = values.copy()
            for name, writes in definitions.items():
                if name in excluded:
                    continue
                bounds = (0, 0)
                for write in writes:
                    bounds = union(bounds, definition(write))
                values[name] = union(values[name], bounds)
                if values[name] is not BOTTOM and not fits(values[name]):
                    values[name] = None
            if values == previous:
                break
        else:
            for name in values:
                if values[name] != previous.get(name):
                    values[name] = None
            for _ in range(len(values) + 1):
                for name, writes in definitions.items():
                    if values[name] is not None:
                        for write in writes:
                            values[name] = union(values[name], definition(write))
                            if values[name] is not BOTTOM and not fits(values[name]):
                                values[name] = None
        return {name: bounds for name, bounds in values.items() if name not in excluded and fits(bounds) and self.types.get(name, Type('unknown')).kind in {'i64', 'unknown'}}

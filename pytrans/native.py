from __future__ import annotations

import ast
import json

from .translator import CPP_KEYWORDS, TranspileError


BUILTINS = {"input", "int", "str", "bool", "range", "print", "map", "len", "abs"}


class NativeEmitter:
    """Emit a statically typed scalar Python subset as standalone C++20."""

    def __init__(self, source: str, filename: str, int_type: str) -> None:
        if int_type not in {"int", "long long"}:
            raise ValueError("int_type must be 'int' or 'long long'")
        self.filename = filename
        self.int_type = int_type
        self.bits = 32 if int_type == "int" else 64
        try:
            self.tree = ast.parse(source, filename)
        except SyntaxError as error:
            raise TranspileError(str(error), filename=filename) from error
        self.types = {}
        self.names = {}
        self.used_names = {"cin", "cout", "std", "string", "getline", "abs", "ios", "numeric_limits", "streamsize", "runtime_error", "main"}
        self.source_names = {node.id for node in ast.walk(self.tree) if isinstance(node, ast.Name)}
        self.headers = {"iostream"}
        self.lines = []
        self.indent = 1
        self.serial = 0
        self.loops = []
        self.local_targets = {}
        self._collect(self.tree.body)
        for node in ast.walk(self.tree):
            if isinstance(node, ast.For):
                target = node.target.id
                inside = {id(child) for child in ast.walk(node)}
                elsewhere = any(
                    isinstance(child, ast.Name) and child.id == target and id(child) not in inside
                    for child in ast.walk(self.tree)
                )
                reassigned = target in self._assigned(node.body) or any(
                    isinstance(child, ast.Name) and child.id == target
                    for statement in node.orelse for child in ast.walk(statement)
                )
                if not elsewhere and not reassigned:
                    self.local_targets[id(node)] = target

    def _error(self, message, node):
        raise TranspileError(f"native: {message}", node, self.filename)

    def _name(self, name):
        if name not in self.names:
            if name.isascii() and name not in CPP_KEYWORDS and "__" not in name and not (name.startswith("_") and len(name) > 1 and name[1].isupper()):
                base = name
            else:
                base = name + "_" if name in CPP_KEYWORDS else "variable_" + name.strip("_").replace("__", "_") if name.isascii() else "variable"
            candidate = base
            while candidate in self.used_names or (candidate != name and candidate in self.source_names):
                candidate += "_"
            self.used_names.add(candidate)
            self.names[name] = candidate
        return self.names[name]

    def _temporary(self, purpose):
        self.serial += 1
        name = f"{purpose}_{self.serial}"
        while name in self.used_names:
            name += "_"
        self.used_names.add(name)
        return name

    def _line(self, text):
        self.lines.append("    " * self.indent + text)

    def _cpp_type(self, kind):
        if kind == "int":
            return self.int_type
        if kind == "str":
            self.headers.add("string")
            return "string"
        return "bool"

    def _bind(self, target, kind):
        if not isinstance(target, ast.Name):
            self._error("assignment targets must be scalar names", target)
        if target.id in BUILTINS:
            self._error("rebinding builtins is unsupported", target)
        previous = self.types.get(target.id)
        if previous is not None and previous != kind:
            self._error(f"variable {target.id!r} changes type from {previous} to {kind}", target)
        self.types[target.id] = kind
        self._name(target.id)

    def _call(self, node, name):
        return isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == name

    def _input(self, node):
        return self._call(node, "input") and not node.args and not node.keywords

    def _integer_input(self, node):
        return self._call(node, "int") and len(node.args) == 1 and not node.keywords and self._input(node.args[0])

    def _split_input(self, node):
        if not self._call(node, "map") or len(node.args) != 2 or node.keywords:
            return False
        cast, split = node.args
        return (
            isinstance(cast, ast.Name) and cast.id == "int"
            and isinstance(split, ast.Call) and not split.args and not split.keywords
            and isinstance(split.func, ast.Attribute) and split.func.attr == "split"
            and self._input(split.func.value)
        )

    def _assignment_types(self, target, value):
        if isinstance(target, (ast.Tuple, ast.List)):
            if not target.elts:
                self._error("empty unpacking targets are unsupported", target)
            if self._split_input(value):
                for item in target.elts:
                    self._bind(item, "int")
            elif isinstance(value, (ast.Tuple, ast.List)) and len(target.elts) == len(value.elts):
                kinds = [self._type(item) for item in value.elts]
                for item, kind in zip(target.elts, kinds):
                    self._bind(item, kind)
            else:
                self._error("unpacking requires map(int, input().split()) or a fixed scalar tuple", value)
        else:
            kind = "int" if self._integer_input(value) else "str" if self._input(value) else self._type(value)
            self._bind(target, kind)

    def _collect(self, body):
        for node in body:
            if isinstance(node, ast.Assign):
                for target in node.targets:
                    self._assignment_types(target, node.value)
            elif isinstance(node, ast.AnnAssign):
                if node.value is None or not isinstance(node.annotation, ast.Name) or node.annotation.id not in {"int", "str", "bool"}:
                    self._error("annotations require an initialized int, str, or bool", node)
                self._assignment_types(node.target, node.value)
                if self.types[node.target.id] != node.annotation.id:
                    self._error("annotation disagrees with inferred type", node)
            elif isinstance(node, ast.AugAssign):
                if not isinstance(node.target, ast.Name) or node.target.id not in self.types:
                    self._error("augmented assignment requires a previously typed variable", node)
                self._type(ast.copy_location(ast.BinOp(node.target, node.op, node.value), node))
            elif isinstance(node, ast.If):
                self._collect(node.body)
                self._collect(node.orelse)
            elif isinstance(node, ast.For):
                if not self._call(node.iter, "range") or node.iter.keywords or not 1 <= len(node.iter.args) <= 3:
                    self._error("for currently supports range(stop), range(start, stop[, step])", node)
                for arg in node.iter.args:
                    if self._type(arg) != "int":
                        self._error("range arguments must be integers", arg)
                self._bind(node.target, "int")
                self._collect(node.body)
                self._collect(node.orelse)
            elif isinstance(node, (ast.Expr, ast.Pass, ast.Break, ast.Continue)):
                pass
            else:
                self._error(f"unsupported statement {type(node).__name__}", node)

    def _type(self, node):
        if isinstance(node, ast.Constant):
            if type(node.value) is bool:
                return "bool"
            if type(node.value) is int:
                if not -(1 << (self.bits - 1)) <= node.value < (1 << (self.bits - 1)):
                    self._error(f"integer literal is outside signed {self.bits}-bit range", node)
                return "int"
            if type(node.value) is str:
                return "str"
        elif isinstance(node, ast.Name):
            if node.id not in self.types:
                self._error(f"cannot infer type of {node.id!r}", node)
            return self.types[node.id]
        elif isinstance(node, ast.UnaryOp):
            literal = self._literal_int(node)
            if literal is not None and -(1 << (self.bits - 1)) <= literal < (1 << (self.bits - 1)):
                return "int"
            kind = self._type(node.operand)
            if isinstance(node.op, ast.Not):
                return "bool"
            if kind == "int" and isinstance(node.op, (ast.UAdd, ast.USub, ast.Invert)):
                return "int"
        elif isinstance(node, ast.BinOp):
            left, right = self._type(node.left), self._type(node.right)
            if left == right == "str" and isinstance(node.op, ast.Add):
                return "str"
            if left == right == "int" and isinstance(node.op, (ast.Add, ast.Sub, ast.Mult, ast.FloorDiv, ast.Mod, ast.BitAnd, ast.BitOr, ast.BitXor)):
                return "int"
        elif isinstance(node, ast.Compare):
            kinds = [self._type(item) for item in [node.left, *node.comparators]]
            if len(set(kinds)) == 1 and all(isinstance(op, (ast.Eq, ast.NotEq, ast.Lt, ast.LtE, ast.Gt, ast.GtE)) for op in node.ops):
                return "bool"
        elif isinstance(node, ast.BoolOp):
            kinds = [self._type(item) for item in node.values]
            if len(set(kinds)) == 1:
                return kinds[0]
        elif isinstance(node, ast.IfExp):
            left, right = self._type(node.body), self._type(node.orelse)
            if left == right:
                return left
        elif isinstance(node, ast.Call) and not node.keywords and len(node.args) == 1:
            kind = self._type(node.args[0])
            if self._call(node, "bool"):
                return "bool"
            if self._call(node, "len") and kind == "str":
                return "int"
            if self._call(node, "abs") and kind == "int":
                return "int"
        self._error(f"unsupported expression {type(node).__name__}", node)

    def _assigned(self, body):
        return {child.id for node in body for child in ast.walk(node) if isinstance(child, ast.Name) and isinstance(child.ctx, ast.Store)}

    def _expr(self, node, initialized):
        kind = self._type(node)
        if isinstance(node, ast.Name):
            if node.id not in initialized:
                self._error(f"variable {node.id!r} may be uninitialized on this path", node)
            return self._name(node.id)
        if isinstance(node, ast.Constant):
            if kind == "bool":
                return "true" if node.value else "false"
            if kind == "int":
                return str(node.value) + ("LL" if self.int_type == "long long" else "")
            self.headers.add("string")
            length = f", {len(node.value.encode('utf-8'))}" if "\x00" in node.value else ""
            return "string(" + self._string(node.value) + length + ")"
        if isinstance(node, ast.UnaryOp):
            if self._literal_int(node) == -(1 << (self.bits - 1)):
                suffix = "LL" if self.int_type == "long long" else ""
                return f"(-{(1 << (self.bits - 1)) - 1}{suffix} - 1{suffix})"
            if isinstance(node.op, ast.Not):
                return f"(!{self._condition(node.operand, initialized)})"
            op = {ast.UAdd: "+", ast.USub: "-", ast.Invert: "~"}[type(node.op)]
            return f"({op}{self._expr(node.operand, initialized)})"
        if isinstance(node, ast.BinOp):
            left, right = self._expr(node.left, initialized), self._expr(node.right, initialized)
            if isinstance(node.op, (ast.FloorDiv, ast.Mod)):
                self.headers.add("stdexcept")
                a, b, r = (self._temporary(purpose) for purpose in ("dividend", "divisor", "remainder"))
                result = f"{a} / {b} - ({r} != 0 && (({r} < 0) != ({b} < 0)))" if isinstance(node.op, ast.FloorDiv) else f"{r} + (({r} != 0 && (({r} < 0) != ({b} < 0))) ? {b} : 0)"
                return f'([&]() {{ const {self.int_type} {a} = {left}; const {self.int_type} {b} = {right}; if ({b} == 0) throw runtime_error("division by zero"); const {self.int_type} {r} = {a} % {b}; return {result}; }}())'
            op = {ast.Add: "+", ast.Sub: "-", ast.Mult: "*", ast.BitAnd: "&", ast.BitOr: "|", ast.BitXor: "^"}[type(node.op)]
            return f"({left} {op} {right})"
        if isinstance(node, ast.Compare):
            values = [self._expr(item, initialized) for item in [node.left, *node.comparators]]
            ops = {ast.Eq: "==", ast.NotEq: "!=", ast.Lt: "<", ast.LtE: "<=", ast.Gt: ">", ast.GtE: ">="}
            return "(" + " && ".join(f"({a} {ops[type(op)]} {b})" for a, op, b in zip(values, node.ops, values[1:])) + ")"
        if isinstance(node, ast.BoolOp):
            if kind == "bool":
                return self._condition(node, initialized)
            first, *rest = node.values
            expression = self._expr(first, initialized)
            for item in rest:
                value = self._temporary("operand")
                truth = f"!{value}.empty()" if kind == "str" else value
                tail = self._expr(item, initialized)
                choices = f"{tail} : {value}" if isinstance(node.op, ast.And) else f"{value} : {tail}"
                expression = f"([&]() {{ const {self._cpp_type(kind)} {value} = {expression}; return {truth} ? {choices}; }}())"
            return expression
        if isinstance(node, ast.IfExp):
            return f"({self._condition(node.test, initialized)} ? {self._expr(node.body, initialized)} : {self._expr(node.orelse, initialized)})"
        if self._call(node, "bool"):
            return self._condition(node.args[0], initialized)
        if self._call(node, "len"):
            return f"static_cast<{self.int_type}>({self._expr(node.args[0], initialized)}.size())"
        if self._call(node, "abs"):
            self.headers.add("cstdlib")
            return f"abs({self._expr(node.args[0], initialized)})"
        self._error("unsupported expression", node)

    def _condition(self, node, initialized):
        if isinstance(node, ast.BoolOp):
            op = " && " if isinstance(node.op, ast.And) else " || "
            return "(" + op.join(self._condition(item, initialized) for item in node.values) + ")"
        value = self._expr(node, initialized)
        return f"(!{value}.empty())" if self._type(node) == "str" else value

    def _string(self, value):
        if any(0xD800 <= ord(char) <= 0xDFFF for char in value):
            self._error("surrogate characters are unsupported", self.tree)
        result = json.dumps(value, ensure_ascii=False)
        for number in range(32):
            result = result.replace(f"\\u{number:04x}", f"\\{number:03o}")
        return result

    def _discard_line(self):
        self.headers.add("limits")
        self._line("cin.ignore(numeric_limits<streamsize>::max(), '\\n');")

    def _assign(self, target, value, initialized):
        if isinstance(target, (ast.Tuple, ast.List)):
            if self._split_input(value):
                self._line("cin >> " + " >> ".join(self._name(item.id) for item in target.elts) + ";")
                self._discard_line()
            else:
                temporaries = []
                for item in value.elts:
                    name = self._temporary("unpack")
                    self._line(f"const {self._cpp_type(self._type(item))} {name} = {self._expr(item, initialized)};")
                    temporaries.append(name)
                for item, name in zip(target.elts, temporaries):
                    self._line(f"{self._name(item.id)} = {name};")
            initialized.update(item.id for item in target.elts)
        elif self._integer_input(value):
            self._line(f"cin >> {self._name(target.id)};")
            self._discard_line()
            initialized.add(target.id)
        elif self._input(value):
            name = self._name(target.id)
            self._line(f"getline(cin, {name});")
            initialized.add(target.id)
        else:
            self._line(f"{self._name(target.id)} = {self._expr(value, initialized)};")
            initialized.add(target.id)

    def _print(self, node, initialized):
        keywords = {}
        for keyword in node.keywords:
            if keyword.arg not in {"sep", "end", "flush"}:
                self._error("print supports sep, end, and flush", node)
            if keyword.arg in keywords:
                self._error("duplicate print keyword", node)
            keywords[keyword.arg] = keyword.value
        sep = keywords.get("sep", ast.Constant(" "))
        end = keywords.get("end", ast.Constant("\n"))
        if self._type(sep) != "str" or self._type(end) != "str":
            self._error("print sep and end must be strings", node)
        def output(item):
            if isinstance(item, ast.Constant) and isinstance(item.value, str) and "\x00" not in item.value:
                return self._string(item.value)
            expression = self._expr(item, initialized)
            if self._type(item) == "bool":
                return f'({expression} ? "True" : "False")'
            return expression
        values = []
        for index, item in enumerate(node.args):
            if index:
                values.append(output(sep))
            values.append(output(item))
        values.append(output(end))
        self._line("cout << " + " << ".join(values) + ";")
        if "flush" in keywords:
            self._line(f"if ({self._condition(keywords['flush'], initialized)}) cout.flush();")

    def _literal_int(self, node):
        if isinstance(node, ast.Constant) and type(node.value) is int:
            return node.value
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            value = self._literal_int(node.operand)
            if value is not None:
                return -value if isinstance(node.op, ast.USub) else value
        return None

    def _for(self, node, initialized):
        args = node.iter.args
        start = ast.Constant(0) if len(args) == 1 else args[0]
        stop = args[0] if len(args) == 1 else args[1]
        step = args[2] if len(args) == 3 else ast.Constant(1)
        modified = self._assigned(node.body) | {node.target.id}
        bounds = []
        for purpose, arg in zip(("range_start", "range_stop", "range_step"), (start, stop, step)):
            expression = self._expr(arg, initialized)
            literal = self._literal_int(arg)
            if literal is None and not (isinstance(arg, ast.Name) and arg.id not in modified):
                name = self._temporary(purpose)
                self._line(f"const {self.int_type} {name} = {expression};")
                expression = name
            bounds.append(expression)
        lower, upper, stride = bounds
        step_value = self._literal_int(step)
        if step_value == 0:
            self._error("range step cannot be zero", step)
        if step_value is None:
            self.headers.add("stdexcept")
            self._line(f'if ({stride} == 0) throw runtime_error("range step cannot be zero");')
        local = id(node) in self.local_targets
        counter = self._name(node.target.id) if local else self._temporary("range_index")
        test = f"{counter} < {upper}" if step_value is not None and step_value > 0 else f"{counter} > {upper}" if step_value is not None else f"({stride} > 0 ? {counter} < {upper} : {counter} > {upper})"
        increment = f"++{counter}" if step_value == 1 else f"--{counter}" if step_value == -1 else f"{counter} += {stride}"
        if step_value not in {1, -1}:
            self.headers.add("limits")
            positive = f"{counter} > numeric_limits<{self.int_type}>::max() - ({stride}) ? {upper} : {counter} + ({stride})"
            negative = f"{counter} < numeric_limits<{self.int_type}>::min() - ({stride}) ? {upper} : {counter} + ({stride})"
            update = positive if step_value is not None and step_value > 0 else negative if step_value is not None else f"{stride} > 0 ? ({positive}) : ({negative})"
            increment = f"{counter} = ({update})"
        broken = self._temporary("loop_broken") if node.orelse else None
        if broken:
            self._line(f"bool {broken} = false;")
        self._line(f"for ({self.int_type} {counter} = {lower}; {test}; {increment}) {{")
        self.indent += 1
        if not local:
            self._line(f"{self._name(node.target.id)} = {counter};")
        self.loops.append(broken)
        self._block(node.body, initialized | {node.target.id})
        self.loops.pop()
        self.indent -= 1
        self._line("}")
        result = set(initialized)
        literals = [self._literal_int(arg) for arg in (start, stop, step)]
        if all(value is not None for value in literals) and (literals[0] < literals[1] if literals[2] > 0 else literals[0] > literals[1]):
            result.add(node.target.id)
        if node.orelse:
            self._line(f"if (!{broken}) {{")
            self.indent += 1
            after_else = self._block(node.orelse, result)
            self.indent -= 1
            self._line("}")
            if not any(isinstance(child, ast.Break) for item in node.body for child in ast.walk(item)):
                result = after_else
        return result

    def _block(self, body, initialized):
        initialized = set(initialized)
        for node in body:
            if isinstance(node, (ast.Assign, ast.AnnAssign)):
                targets = node.targets if isinstance(node, ast.Assign) else [node.target]
                if len(targets) == 1:
                    self._assign(targets[0], node.value, initialized)
                else:
                    if self._input(node.value) or self._integer_input(node.value):
                        self._assign(targets[0], node.value, initialized)
                        expression = self._name(targets[0].id)
                    else:
                        expression = self._temporary("assigned_value")
                        self._line(f"const {self._cpp_type(self._type(node.value))} {expression} = {self._expr(node.value, initialized)};")
                    for target in targets:
                        if not isinstance(target, ast.Name):
                            self._error("chained assignments require scalar names", target)
                        self._line(f"{self._name(target.id)} = {expression};")
                        initialized.add(target.id)
            elif isinstance(node, ast.AugAssign):
                value = ast.copy_location(ast.BinOp(node.target, node.op, node.value), node)
                self._line(f"{self._name(node.target.id)} = {self._expr(value, initialized)};")
            elif isinstance(node, ast.If):
                self._line(f"if ({self._condition(node.test, initialized)}) {{")
                self.indent += 1
                left = self._block(node.body, initialized)
                self.indent -= 1
                if node.orelse:
                    self._line("} else {")
                    self.indent += 1
                    right = self._block(node.orelse, initialized)
                    self.indent -= 1
                else:
                    right = initialized
                self._line("}")
                initialized = left & right
            elif isinstance(node, ast.For):
                initialized = self._for(node, initialized)
            elif isinstance(node, ast.Expr):
                if self._call(node.value, "print"):
                    self._print(node.value, initialized)
                elif isinstance(node.value, ast.Constant) and isinstance(node.value.value, str):
                    pass
                else:
                    self._error("expression statements must be print calls", node)
            elif isinstance(node, (ast.Break, ast.Continue)):
                if not self.loops:
                    self._error("break/continue must be inside a loop", node)
                if isinstance(node, ast.Break) and self.loops[-1]:
                    self._line(f"{self.loops[-1]} = true;")
                self._line("break;" if isinstance(node, ast.Break) else "continue;")
            elif isinstance(node, ast.Pass):
                self._line(";")
        return initialized

    def emit(self):
        local = set(self.local_targets.values())
        declarations = [f"    {self._cpp_type(kind)} {self._name(name)};" for name, kind in self.types.items() if name not in local]
        self._block(self.tree.body, set())
        includes = "\n".join(f"#include <{header}>" for header in sorted(self.headers))
        body = "\n".join(declarations + self.lines)
        return f"{includes}\n\nusing namespace std;\n\nint main() {{\n    ios::sync_with_stdio(false);\n    cin.tie(nullptr);\n{body}\n    return 0;\n}}\n"


def transpile_native_cpp(source: str, filename: str = "main.py", *, int_type: str = "Int", registry=None) -> str:
    """Translate the existing typed Python subset to native-oriented C++.

    Int selects checked int64 arithmetic with Boost promotion. int64 selects
    full native syntax with fixed-width storage. Explicit int or long long
    retains the original scalar-only fixed-width emitter. The native
    library registry replaces exact AST fingerprints; unmatched definitions use
    ordinary rule-based generation.
    """
    if int_type == "Int":
        from .native_full import transpile_native_full_cpp
        return transpile_native_full_cpp(source,filename,registry=registry)
    if int_type == "int64":
        from .native_full import transpile_native_full_cpp
        return transpile_native_full_cpp(source,filename,registry=registry,integer_mode='fixed64')
    return NativeEmitter(source, filename, int_type).emit()

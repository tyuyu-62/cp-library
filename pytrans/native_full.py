from __future__ import annotations

import ast
import re
from dataclasses import replace
from .translator import CPP_KEYWORDS, CppGenerator, TranspileError, Type, sanitize_name, element_type, cpp_string, list_type
from .int_runtime import INT_RUNTIME
from .library_registry import LibraryRegistry
from .integer_bounds import IntegerBounds, NativeInteger, fits


_TOKEN = re.compile(r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'', re.S)


def mask_cpp(source):
    return _TOKEN.sub(lambda m: re.sub(r'[^\n]', ' ', m[0]), source)



def replace_identifiers(source, replacements):
    matches = list(re.finditer(r'\b\w+\b', mask_cpp(source)))
    for match in reversed(matches):
        replacement = replacements.get(match[0])
        if replacement is not None:
            source = source[:match.start()] + replacement + source[match.end():]
    return source

def runtime_units(source):
    source = re.sub(r'^#include[^\n]*\n|^using[^\n]*\n', '', source, flags=re.M)
    masked = mask_cpp(source)
    units = []
    start = 0
    depth = parens = 0
    for index, char in enumerate(masked):
        if char == '(':
            parens += 1
        elif char == ')':
            parens -= 1
        elif char == '{':
            depth += 1
        elif char == '}':
            depth -= 1
            if depth == 0 and parens == 0:
                next_text = masked[index + 1:].lstrip()
                if not next_text.startswith(';'):
                    units.append(source[start:index + 1].strip())
                    start = index + 1
        elif char == ';' and depth == 0 and parens == 0:
            units.append(source[start:index + 1].strip())
            start = index + 1
    if source[start:].strip():
        raise ValueError('incomplete C++ runtime unit')
    return [unit for unit in units if unit]


def unit_name(unit):
    clean = mask_cpp(unit)
    declaration = re.search(r'\b(?:struct|class)\s+(PyDict|PyDefaultDict|FastScanner|py_local_list)\b', clean)
    if declaration:
        return declaration[1]
    function = re.search(r'\b(py_\w+)\s*\(', clean)
    if function:
        return function[1]
    if re.search(r'\bpy_input\b', clean):
        return 'py_input'
    raise ValueError('unknown C++ runtime unit: ' + unit[:120])


def prune_runtime(runtime, body):
    units = runtime_units(runtime)
    groups = {}
    families = set(re.findall(r'\b(?:PyDict|PyDefaultDict|py_local_list|deque|vector|set|tuple)\b', mask_cpp(body)))
    if re.search(r'\bmake_tuple\s*\(|\bget\s*<', mask_cpp(body)):
        families.add('tuple')
    for unit in units:
        header = unit.split('{', 1)[0]
        name = unit_name(unit)
        overload = re.search(r'\b(py_\w+)\s*\(', header)
        if overload:
            parameters = header[overload.end():]
            unused = set(re.findall(r'\b(?:PyDict|PyDefaultDict|py_local_list|deque|set|tuple)\b', parameters)) - families
            if unused:
                continue
        groups.setdefault(name, []).append(unit)
    masked_body = mask_cpp(body)
    needed = set(re.findall(r'\b\w+\b', masked_body)) & groups.keys()
    needed = {name for name in needed if not name.startswith('py_') or name == 'py_input' or re.search(r'\b' + re.escape(name) + r'\s*(?:\(|<(?![<=])[^;{}]*>\s*\()', masked_body)}
    while True:
        dependencies = set()
        for name in needed:
            for unit in groups[name]:
                dependencies.update(re.findall(r'\b\w+\b', mask_cpp(unit)))
        expanded = needed | (dependencies & groups.keys())
        if expanded == needed:
            break
        needed = expanded
    selected = [unit for unit in units if unit_name(unit) in needed and unit in groups[unit_name(unit)]]
    headers = '\n'.join(re.findall(r'^#include[^\n]*', runtime, re.M))
    return headers + '\nusing namespace std;\nusing u64 = unsigned long long;\n' + '\n\n'.join(selected) + '\n'


class NativeCppGenerator(CppGenerator):
    integer_cpp = 'Int'
    def __init__(self, source, filename='main.py', registry=None):
        self.truth_values = {}
        self.emitting_main = False
        super().__init__(source, filename)
        global_analysis = self.integer_analysis(self.analyzer.global_types)
        bounds = global_analysis.scope(self.tree.body)
        written_globals = {name for node in ast.walk(self.tree) if isinstance(node, ast.Global) for name in node.names}
        for name, interval in bounds.items():
            if name not in written_globals:
                self.analyzer.global_types[name] = NativeInteger('i64', bounds=global_analysis.constants.get(name, interval))
        def globals_for(function):
            types = dict(function.value_types)
            local = self.analyzer.bound_names(function.node.body) | set(function.parameter_types)
            for name, value_type in self.analyzer.global_types.items():
                if name not in local and isinstance(value_type, NativeInteger):
                    types[name] = value_type
            return replace(function, value_types=types)
        self.program_ir = replace(self.program_ir, functions={name: globals_for(function) for name, function in self.program_ir.functions.items()}, methods={name: globals_for(function) for name, function in self.program_ir.methods.items()})
        original_infer = self.analyzer.infer
        def container_type(value):
            if isinstance(value, NativeInteger):
                return Type('i64')
            return replace(value, args=tuple(container_type(arg) for arg in value.args)) if value.args else value
        def native_infer(node, env):
            value = original_infer(node, env)
            return container_type(value) if value.args else value
        self.analyzer.infer = native_infer
        self.registry = registry if isinstance(registry, LibraryRegistry) else LibraryRegistry.load(registry)
        self.replacements = self.registry.match(ast.parse(source,filename=filename))
        self.nonnegative_names = set()
        self.cpp_libraries = {}

    def _function_ir(self, *args, **kwargs):
        value = super()._function_ir(*args, **kwargs)
        analysis = self.integer_analysis(value.value_types)
        bounds = analysis.scope(value.node.body, value.parameter_types)
        types = dict(value.value_types)
        for name, interval in bounds.items():
            if name not in value.rebound_names:
                types[name] = NativeInteger('i64', bounds=interval)
        local_names = self.analyzer.bound_names(value.node.body) | set(value.parameter_types)
        for name, interval in analysis.constants.items():
            if name not in local_names and name in types:
                types[name] = NativeInteger('i64', bounds=interval)
        locals_ir = {name: replace(local, value_type=types[name]) for name, local in value.locals.items()}
        return replace(value, value_types=types, locals=locals_ir, local_lists=frozenset())

    def integer_analysis(self, types=None):
        constants = {}
        for name, assignment in self.optimization_planner.constant_globals.items():
            value = ast.literal_eval(assignment.value)
            if isinstance(value, int) and not isinstance(value, bool):
                if types is not None or isinstance(self.context.env.get(name), NativeInteger):
                    constants[name] = (value, value)
        return IntegerBounds(self.analyzer, types if types is not None else self.context.env, constants)

    def native_integer(self, node):
        if isinstance(node, ast.Constant):
            return isinstance(node.value, int) and fits((node.value, node.value))
        if isinstance(node, ast.Name):
            return isinstance(self.context.env.get(node.id), NativeInteger)
        if isinstance(node, ast.UnaryOp):
            return self.native_integer(node.operand) and fits(self.integer_analysis().expression(node))
        if isinstance(node, ast.BinOp):
            return self.native_integer(node.left) and self.native_integer(node.right) and self.native_operation(node)
        if isinstance(node, ast.IfExp):
            return self.native_integer(node.body) and self.native_integer(node.orelse)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and self.integer_analysis().builtin(node.func.id):
            if node.func.id == 'len' and len(node.args) == 1:
                return self.integer_analysis().length(node.args[0]) is not None
            if node.func.id in {'min', 'max'} and len(node.args) == 2 and not node.keywords:
                return all(self.native_integer(arg) for arg in node.args)
            if node.func.id == 'abs' and len(node.args) == 1:
                return self.native_integer(node.args[0]) and fits(self.integer_analysis().expression(node))
        return False

    def native_operation(self, node):
        analysis = self.integer_analysis()
        if not fits(analysis.expression(node)):
            return False
        if isinstance(node.op, (ast.FloorDiv, ast.Mod)):
            left, right = analysis.expression(node.left), analysis.expression(node.right)
            return fits(left) and fits(right) and left[0] >= 0 and right[0] > 0
        if isinstance(node.op, (ast.LShift, ast.RShift)):
            left, right = analysis.expression(node.left), analysis.expression(node.right)
            return fits(left) and fits(right) and left[0] >= 0 and 0 <= right[0] <= right[1] < 64
        return isinstance(node.op, (ast.Add, ast.Sub, ast.Mult, ast.BitAnd, ast.BitOr, ast.BitXor))

    def emit(self, line=''):
        for value, truth in sorted(self.truth_values.items(), key=lambda item: -len(item[0])):
            line = line.replace(f'py_truth({value})', truth)
        super().emit(line)

    def truth_expression(self, node, value):
        value_type = self.analyzer.infer(node, self.context.env)
        if value_type.kind == 'bool':
            return value
        if value_type.kind in {'i64','f64','char'}:
            return f'({value} != 0)'
        if value_type.kind in {'list','deque','map','default_map','set','string','bytes'}:
            if isinstance(node, ast.Constant) and isinstance(node.value, (str, bytes)):
                return 'true' if node.value else 'false'
            return f'(!({value}).empty())'
        return f'py_truth({value})'

    def print_statement(self, node):
        if any(isinstance(arg, ast.Starred) for arg in node.args) or any(self.analyzer.infer(arg,self.context.env).kind in {'optional','class'} for arg in node.args):
            return super().print_statement(node)
        options = {keyword.arg: keyword.value for keyword in node.keywords}
        if any(key not in {'sep','end','flush'} or not isinstance(value,ast.Constant) for key,value in options.items()):
            return super().print_statement(node)
        sep = options['sep'].value if 'sep' in options else ' '
        end = options['end'].value if 'end' in options else '\n'
        if not isinstance(sep,str) or not isinstance(end,str):
            return super().print_statement(node)
        parts = []
        for index, arg in enumerate(node.args):
            if index:
                parts.append(cpp_string(sep))
            value = self.expr(arg)
            kind = self.analyzer.infer(arg,self.context.env).kind
            if kind == 'bool':
                value = f'({value} ? "True" : "False")'
            elif kind in {'list','deque','set','map','default_map','tuple'}:
                value = f'py_string({value})'
            if any(isinstance(child,ast.Call) for child in ast.walk(arg)):
                temporary = self.context.temporary('print_value')
                self.emit(f'auto {temporary} = {value};')
                value = temporary
            parts.append(value)
        parts.append(cpp_string(end))
        self.emit('cout << ' + ' << '.join(parts) + ';')
        if options.get('flush') and options['flush'].value:
            self.emit('cout.flush();')

    def _emit_runtime(self):
        before = len(self.lines)
        super()._emit_runtime()
        self.runtime = '\n'.join(self.lines[before:])
        del self.lines[before:]

    def _emit_globals(self):
        before = len(self.lines)
        super()._emit_globals()
        for index in range(before,len(self.lines)):
            if self.lines[index].startswith('constexpr i64 '):
                self.lines[index] = self.lines[index].replace('constexpr i64 ','constexpr long long ',1)

    def _emit_main(self):
        self.emitting_main = True
        super()._emit_main()
        self.emitting_main = False

    def statement(self, node):
        if self.emitting_main and self.context is not None:
            for name, value_type in self.analyzer.global_types.items():
                if isinstance(value_type, NativeInteger):
                    self.context.env[name] = value_type
        if isinstance(node, ast.AugAssign) and isinstance(node.target, ast.Name) and isinstance(self.context.env.get(node.target.id), NativeInteger):
            expression = ast.BinOp(left=node.target, op=node.op, right=node.value)
            self.assign_target(node.target, self.expr(expression))
            return
        super().statement(node)

    def _emit_class_declaration(self, name, node):
        entry = self.replacements.get(('class', name))
        if entry:
            token = f'CP_NATIVE_LIBRARY_{len(self.cpp_libraries)};'
            self.cpp_libraries[token] = entry.cpp_source
            self.lines.append(token)
        else:
            super()._emit_class_declaration(name, node)

    def _emit_class_methods(self, name, node):
        if ('class', name) not in self.replacements:
            super()._emit_class_methods(name, node)

    def _emit_prototype(self, name, node):
        if ('function', name) not in self.replacements:
            super()._emit_prototype(name, node)

    def _emit_function(self, name, node):
        entry = self.replacements.get(('function', name))
        if entry:
            token = f'CP_NATIVE_LIBRARY_{len(self.cpp_libraries)};'
            self.cpp_libraries[token] = entry.cpp_source
            self.lines.append(token)
        else:
            super()._emit_function(name, node)

    def nonnegative(self, node):
        if isinstance(node, ast.Constant) and isinstance(node.value, int):
            return node.value >= 0
        if isinstance(node, ast.Name):
            value_type = self.context.env.get(node.id)
            return node.id in self.nonnegative_names or isinstance(value_type, NativeInteger) and value_type.bounds[0] >= 0
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            return node.func.id in {'len', 'abs'}
        if isinstance(node, ast.BinOp):
            if isinstance(node.op, (ast.Add, ast.Mult, ast.BitAnd, ast.BitOr)):
                return self.nonnegative(node.left) and self.nonnegative(node.right)
            if isinstance(node.op, ast.Mod) and isinstance(node.right, ast.Constant) and isinstance(node.right.value, int):
                return node.right.value > 0
        return False

    def for_statement(self, node):
        previous = self.nonnegative_names
        self.nonnegative_names = set(previous)
        if isinstance(node.target, ast.Name) and isinstance(node.iter, ast.Call) and isinstance(node.iter.func, ast.Name):
            args = node.iter.args
            nonnegative = node.iter.func.id == 'range' and (len(args) == 1 or len(args) >= 2 and self.nonnegative(args[0]) and (len(args) == 2 or isinstance(args[2], ast.Constant) and args[2].value > 0))
            assigned = any(isinstance(child, ast.Name) and isinstance(child.ctx, ast.Store) and child.id == node.target.id for statement in node.body for child in ast.walk(statement))
            if nonnegative and not assigned:
                self.nonnegative_names.add(node.target.id)
        if not self.enumerate_for(node):
            super().for_statement(node)
        self.nonnegative_names = previous

    def enumerate_for(self, node):
        call = node.iter
        analysis = self.integer_analysis()
        if not isinstance(call, ast.Call) or not isinstance(call.func, ast.Name) or call.func.id != 'enumerate' or not analysis.builtin('enumerate') or call.keywords or not 1 <= len(call.args) <= 2 or not isinstance(node.target, (ast.Tuple, ast.List)) or len(node.target.elts) != 2:
            return False
        length = analysis.length(call.args[0])
        if length is None or analysis.kind(call.args[0]) not in {'list', 'deque', 'string', 'bytes'}:
            return False
        source = self.context.temporary('enumerate_source')
        self.emit(f'auto&& {source} = ({self.expr(call.args[0])});')
        start_node = call.args[1] if len(call.args) == 2 else ast.Constant(value=0)
        start_bounds = analysis.expression(start_node)
        index_bounds = (start_bounds[0], start_bounds[1] + max(0, length[1] - 1)) if fits(start_bounds) else None
        fixed = fits(index_bounds)
        start = self.context.temporary('enumerate_start')
        self.emit(f'{"long long" if fixed else "i64"} {start} = {self.expr(start_node)};')
        offset = self.context.temporary('offset')
        item_type = element_type(self.analyzer.infer(call.args[0], self.context.env))
        alias = isinstance(node.target.elts[1], ast.Name) and item_type.kind in {'list', 'deque', 'set', 'map', 'default_map', 'class'} and isinstance(call.args[0], (ast.Name, ast.Attribute))
        if alias:
            name = node.target.elts[1].id
            self.context.reference_expressions[name] = f'{source}[{offset}]'
        self.emit(f'for (size_t {offset} = 0; {offset} < {source}.size(); ++{offset}) {{')
        self.indent += 1
        index = f'({start} + static_cast<long long>({offset}))' if fixed else f'({self.integer_cpp}({start}) + static_cast<long long>({offset}))'
        self.assign_target(node.target.elts[0], index)
        if not alias:
            item = f'{source}[{offset}]'
            if analysis.kind(call.args[0]) == 'bytes':
                item = f'static_cast<long long>(static_cast<unsigned char>({item}))'
            self.assign_target(node.target.elts[1], item)
        for child in node.body:
            self.statement(child)
        self.indent -= 1
        self.emit('}')
        return True

    def range_for(self, node):
        arguments = self.integer_analysis().range_arguments(node.iter)
        if arguments is None or not all(fits(interval) for interval in arguments):
            return super().range_for(node)
        self.fixed_range_for(node, arguments)

    def fixed_range_for(self, node, arguments):
        if not isinstance(node.target, ast.Name):
            return super().range_for(node)
        expressions = [self.expr(arg) for arg in node.iter.args]
        if len(expressions) == 1:
            expressions.insert(0, '0LL')
        if len(expressions) == 2:
            expressions.append('1LL')
        start, stop, step = [self.context.temporary('range_' + kind) for kind in ('start', 'stop', 'step')]
        for name, expression in zip((start, stop, step), expressions):
            self.emit(f'long long {name} = static_cast<long long>({expression});')
        if arguments[2] in {(1, 1), (-1, -1)}:
            variable = self.context.temporary('range')
            positive = arguments[2] == (1, 1)
            self.emit(f'for (long long {variable} = {start}; {variable} {"<" if positive else ">"} {stop}; {"++" if positive else "--"}{variable}) {{')
        else:
            count = self.context.temporary('range_remaining')
            variable = self.context.temporary('range')
            self.emit(f'if ({step} == 0) throw runtime_error("range() arg 3 must not be zero");')
            distance = f'({step} > 0 ? static_cast<CP_NATIVE_WIDE>({stop}) - {start} : static_cast<CP_NATIVE_WIDE>({start}) - {stop})'
            width = f'({step} > 0 ? static_cast<CP_NATIVE_WIDE>({step}) : -static_cast<CP_NATIVE_WIDE>({step}))'
            delta = self.context.temporary('range_distance')
            stride = self.context.temporary('range_stride')
            self.emit(f'CP_NATIVE_WIDE {delta} = {distance}, {stride} = {width};')
            self.emit(f'unsigned long long {count} = {delta} > 0 ? static_cast<unsigned long long>(({delta} - 1) / {stride} + 1) : 0ULL;')
            self.emit(f'long long {variable} = {start};')
            self.emit(f'for (; {count} > 0; --{count}, ({count} > 0 ? {variable} += {step} : {variable})) {{')
        self.indent += 1
        self.assign_target(node.target, variable)
        for child in node.body:
            self.statement(child)
        self.indent -= 1
        self.emit('}')

    def subscript(self, node, mutable):
        base_type = self.analyzer.infer(node.value, self.context.env)
        if base_type.kind in {'list', 'deque', 'string', 'bytes'} and not isinstance(node.slice, ast.Slice) and self.nonnegative(node.slice):
            value = f'{self.expr(node.value)}.at(static_cast<size_t>({self.expr(node.slice)}))'
            return f'{self.integer_cpp}(static_cast<unsigned char>({value}))' if base_type.kind == 'bytes' else value
        return super().subscript(node, mutable)

    def expr(self, node):
        value = self.native_expr(node)
        if self.context is not None:
            self.truth_values[value] = self.truth_expression(node, value)
        return value

    def native_expr(self, node):
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd, ast.Invert)) and self.native_integer(node.operand):
            operator = '-' if isinstance(node.op, ast.USub) else '+' if isinstance(node.op, ast.UAdd) else '~'
            operand = self.expr(node.operand)
            if not fits(self.integer_analysis().expression(node)):
                operand = f'Int({operand})'
            return f'({operator}{operand})'
        if isinstance(node, ast.IfExp) and self.analyzer.infer(node,self.context.env).kind == 'i64':
            if self.native_integer(node):
                return f'({self.truth_expression(node.test,self.expr(node.test))} ? {self.expr(node.body)} : {self.expr(node.orelse)})'
            return f'({self.truth_expression(node.test,self.expr(node.test))} ? Int({self.expr(node.body)}) : Int({self.expr(node.orelse)}))'
        if isinstance(node, ast.Constant) and isinstance(node.value, int) and not isinstance(node.value, bool):
            value = node.value
            return f'{value}LL' if -(1 << 63) <= value < (1 << 63) else f'Int::parse("{value}")'
        if isinstance(node, ast.BinOp):
            left_type = self.analyzer.infer(node.left, self.context.env)
            right_type = self.analyzer.infer(node.right, self.context.env)
            if left_type.kind == right_type.kind == 'i64' and isinstance(node.op,ast.Mod) and isinstance(node.left,ast.BinOp) and isinstance(node.left.op,ast.Pow) and isinstance(node.left.right,ast.Constant) and isinstance(node.left.right.value,int) and node.left.right.value>=0:
                return f'py_pow_mod({self.value_expr(node.left.left)}, {self.value_expr(node.left.right)}, {self.value_expr(node.right)})'
            if 'f64' in {left_type.kind,right_type.kind}:
                left=f'static_cast<double>({self.value_expr(node.left)})'
                right=f'static_cast<double>({self.value_expr(node.right)})'
                if isinstance(node.op,ast.FloorDiv):
                    return f'std::floor({left} / {right})'
                if isinstance(node.op,ast.Mod):
                    return f'py_float_mod({left}, {right})'
                if isinstance(node.op,ast.Pow):
                    return f'std::pow({left}, {right})'
            if left_type.kind in {'i64','bool'} and right_type.kind in {'i64','bool'} and not isinstance(node.op, ast.Pow):
                if isinstance(node.op, ast.Div):
                    return f'(static_cast<double>({self.value_expr(node.left)}) / static_cast<double>({self.value_expr(node.right)}))'
                operator = '/' if isinstance(node.op, ast.FloorDiv) else self.operator(node.op)
                left = self.value_expr(node.left)
                right = self.value_expr(node.right)
                if self.native_integer(node.left) and self.native_integer(node.right):
                    if self.native_operation(node):
                        return f'({left} {operator} {right})'
                    left = f'Int({left})'
                return f'({left} {operator} {right})'
        return super().expr(node)

    def fixed_integer(self, node):
        return self.native_integer(node)

    def call(self, node):
        if isinstance(node.func, ast.Name) and node.func.id == 'len' and len(node.args) == 1 and self.integer_analysis().builtin('len') and self.integer_analysis().length(node.args[0]) is not None:
            value = super().call(node)
            return value.replace('static_cast<i64>', 'static_cast<long long>')
        if isinstance(node.func, ast.Name) and node.func.id == 'abs' and len(node.args) == 1 and self.integer_analysis().builtin('abs') and self.native_integer(node.args[0]):
            value = self.expr(node.args[0])
            return f'std::abs({value})' if fits(self.integer_analysis().expression(node)) else f'abs(Int({value}))'
        if isinstance(node.func,ast.Attribute) and node.func.attr == 'from_bytes' and isinstance(node.func.value,ast.Name) and node.func.value.id == 'int':
            if not 1 <= len(node.args) <= 2:
                self.fail('int.from_bytes expects one or two positional arguments',node)
            options={keyword.arg:keyword.value for keyword in node.keywords}
            if set(options)-{'signed','byteorder'}:
                self.fail('unsupported int.from_bytes keyword',node)
            order=self.expr(node.args[1]) if len(node.args)==2 else self.expr(options['byteorder']) if 'byteorder' in options else 'string("big")'
            signed=self.expr(options['signed']) if 'signed' in options else 'false'
            return f'py_int_from_bytes({self.expr(node.args[0])}, {order}, {signed})'
        if isinstance(node.func,ast.Attribute) and node.func.attr == 'to_bytes' and self.analyzer.infer(node.func.value,self.context.env).kind == 'i64':
            if not 1 <= len(node.args) <= 2:
                self.fail('to_bytes expects a length and optional byteorder',node)
            options={keyword.arg:keyword.value for keyword in node.keywords}
            if set(options)-{'signed','byteorder'}:
                self.fail('unsupported to_bytes keyword',node)
            order=self.expr(node.args[1]) if len(node.args)==2 else self.expr(options['byteorder']) if 'byteorder' in options else 'string("big")'
            signed=self.expr(options['signed']) if 'signed' in options else 'false'
            return f'py_to_bytes({self.value_expr(node.func.value)}, {self.expr(node.args[0])}, {order}, {signed})'
        function = None
        if isinstance(node.func,ast.Attribute) and isinstance(node.func.value,ast.Name) and node.func.value.id == 'math':
            function = node.func.attr
        elif isinstance(node.func,ast.Name) and node.func.id not in self.context.env:
            origin = self.analyzer.imported_origin(node.func.id)
            if origin and origin[0] == 'math':
                function = origin[1]
        if function == 'isqrt' and len(node.args)==1:
            return f'Int({self.value_expr(node.args[0])}).isqrt()'
        if function in {'floor','ceil'} and len(node.args)==1 and self.analyzer.infer(node.args[0],self.context.env).kind=='i64':
            return f'Int({self.value_expr(node.args[0])})'
        if function in {'sqrt','hypot','sin','cos','acos','atan2','log2','floor','ceil'}:
            arguments = ', '.join(f'static_cast<double>({self.value_expr(arg)})' for arg in node.args)
            value = f'std::{function}({arguments})'
            return f'Int({value})' if function in {'floor','ceil'} else value
        if isinstance(node.func,ast.Name) and node.func.id=='pow' and len(node.args)==2 and any(self.analyzer.infer(arg,self.context.env).kind=='f64' for arg in node.args):
            return 'std::pow(' + ', '.join(f'static_cast<double>({self.value_expr(arg)})' for arg in node.args) + ')'

        if isinstance(node.func, ast.Name) and node.func.id in {'min','max'} and len(node.args) == 2 and not node.keywords:
            types = [self.analyzer.infer(arg, self.context.env) for arg in node.args]
            if all(value.kind == 'i64' for value in types):
                if all(self.native_integer(arg) for arg in node.args):
                    return f'std::{node.func.id}({self.value_expr(node.args[0])}, {self.value_expr(node.args[1])})'
                return f'std::{node.func.id}(Int({self.value_expr(node.args[0])}), Int({self.value_expr(node.args[1])}))'
        if isinstance(node.func, ast.Attribute) and node.func.attr in {'bit_length','bit_count'} and not node.args:
            if self.analyzer.infer(node.func.value, self.context.env).kind == 'i64':
                return f'Int({self.value_expr(node.func.value)}).{node.func.attr}()'
        value = super().call(node)
        if isinstance(node.func,ast.Name) and node.func.id in self.analyzer.functions and node.func.id not in self.context.declared and node.func.id not in self.context.env:
            name = sanitize_name(node.func.id)
            if value.startswith(name + '('):
                value = '::' + value
        return value.replace('std::gcd(', 'integer_gcd(').replace('std::lcm(', 'integer_lcm(')

    def emit_starred_unpack(self, target, value):
        stars = [i for i, child in enumerate(target.elts) if isinstance(child, ast.Starred)]
        if len(stars) != 1:
            self.fail('starred assignment requires exactly one starred target', target)
        value_type = self.analyzer.infer(value, self.context.env)
        if value_type.kind not in {'list','deque'}:
            return super().emit_starred_unpack(target, value)
        index = stars[0]
        suffix = len(target.elts) - index - 1
        temporary = self.context.temporary('unpacked')
        self.emit(f'auto {temporary} = {self.expr(value)};')
        self.emit(f'if ({temporary}.size() < {len(target.elts) - 1}) throw runtime_error("not enough values to unpack");')
        for position, child in enumerate(target.elts[:index]):
            self.assign_unpack_target(child, f'{temporary}[{position}]', element_type(value_type))
        self.assign_unpack_target(target.elts[index].value, f'vector<{element_type(value_type).cpp()}>({temporary}.begin() + {index}, {temporary}.end() - {suffix})', list_type(element_type(value_type)))
        for offset, child in enumerate(target.elts[index+1:]):
            self.assign_unpack_target(child, f'{temporary}[{temporary}.size() - {suffix-offset}]', element_type(value_type))
        return True

    def finish_body(self, body):
        source_names = {node.id for node in ast.walk(self.tree) if isinstance(node, ast.Name)}
        source_names.update(node.arg for node in ast.walk(self.tree) if isinstance(node, ast.arg))
        source_names.update(node.attr for node in ast.walk(self.tree) if isinstance(node, ast.Attribute))
        source_names.update(node.name for node in ast.walk(self.tree) if isinstance(node, (ast.FunctionDef, ast.ClassDef)))
        excluded = set(self.analyzer.global_types) | set(self.analyzer.functions) | set(self.analyzer.classes)
        masked = mask_cpp(body)
        originals = {}
        for name in source_names - excluded - CPP_KEYWORDS:
            escaped = sanitize_name(name)
            if escaped == name or escaped in source_names or not name.isascii():
                continue
            if re.search(r'\b' + re.escape(escaped) + r'\s*\(', masked):
                continue
            if re.search(r'(?<![\w.:>])' + re.escape(name) + r'\b', masked):
                continue
            originals[escaped] = name
        body = replace_identifiers(body, originals)
        identifiers = set(re.findall(r'\b\w+\b', mask_cpp(body)))
        renames = {}
        for name in sorted(identifiers):
            match = re.fullmatch(r'__py_rebound_(.+)_(\d+)', name)
            if match:
                base, number = match.groups()
                base = originals.get(base, base)
                candidate = f'{base}_{number}'
                while candidate in identifiers:
                    number = str(int(number)+1)
                    candidate = f'{base}_{number}'
                identifiers.add(candidate)
                renames[name] = candidate
        body = replace_identifiers(body, renames)
        for token, cpp in self.cpp_libraries.items():
            body = body.replace(token, cpp)
        return body

    def generate(self):
        body = super().generate()
        body = replace_identifiers(body, {'i64':'Int', '__int128':'Int', 'CP_NATIVE_WIDE':'__int128'})
        body = re.sub(r'\bconstexpr Int\b', 'const Int', body)
        body = self.finish_body(body)
        runtime = adapt_runtime(self.runtime)
        runtime += '\ndouble py_float_mod(double a,double b) { if(b==0) throw domain_error("float modulo by zero"); double r=fmod(a,b); if(r!=0 && ((r<0)!=(b<0))) r+=b; return r; }\nInt integer_gcd(Int a, Int b) { a=abs(a); b=abs(b); while(b!=0) { Int r=a%b; a=b; b=r; } return a; }\nInt integer_lcm(Int a, Int b) { return a==0 || b==0 ? Int(0) : abs(a / integer_gcd(a,b) * b); }\n'
        runtime = runtime.replace('integer_gcd', 'py_native_gcd').replace('integer_lcm', 'py_native_lcm')
        body = body.replace('integer_gcd', 'py_native_gcd').replace('integer_lcm', 'py_native_lcm')
        runtime = prune_runtime(runtime, body)
        return (INT_RUNTIME + '\n' if re.search(r'\bInt\b',mask_cpp(runtime+body)) else '') + runtime + '\n' + body


def adapt_runtime(runtime):
    runtime = replace_identifiers(runtime, {'i64':'Int','__int128':'Int'})
    units = runtime_units(runtime)
    replacements = {
        'py_int_from_bytes': """template<class C>
Int py_int_from_bytes(const C& source,const string& byteorder="big",bool signed_value=false) {
    if(byteorder!="big" && byteorder!="little") throw invalid_argument("invalid byteorder");
    vector<unsigned char> bytes; bytes.reserve(source.size());
    for(const auto& value:source) bytes.push_back(static_cast<unsigned char>(value));
    Int::Big result;
    boost::multiprecision::import_bits(result,bytes.begin(),bytes.end(),8,byteorder=="big");
    if(signed_value && !bytes.empty() && ((byteorder=="big" ? bytes.front() : bytes.back()) & 128)) result -= Int::Big(1) << (bytes.size()*8);
    return Int::from_big(std::move(result));
}""",
        'py_to_bytes': """string py_to_bytes(const Int& value,const Int& length,const string& byteorder,bool signed_value=false) {
    if(byteorder!="big" && byteorder!="little") throw invalid_argument("invalid byteorder");
    size_t size=static_cast<size_t>(length);
    if(size>numeric_limits<size_t>::max()/8) throw length_error("byte count too large");
    if(size==0) { if(value==0 || signed_value && value==-1) return string(); throw overflow_error("integer does not fit byte count"); }
    if(!signed_value && value<0) throw overflow_error("negative integer cannot use unsigned bytes");
    Int needed=(value<0 ? ~value : value).bit_length();
    if(needed > Int(size)*8-(signed_value ? 1 : 0)) throw overflow_error("integer does not fit byte count");
    Int::Big bits=value.big_value(); if(bits<0) bits+=Int::Big(1) << (size*8);
    vector<unsigned char> bytes; boost::multiprecision::export_bits(bits,back_inserter(bytes),8,false);
    string result(size,char(0));
    for(size_t i=0;i<bytes.size();++i) result[byteorder=="big" ? size-1-i : i]=static_cast<char>(bytes[i]);
    return result;
}""",
        'py_floor_div': 'Int py_floor_div(Int a, Int b) { return a / b; }',
        'py_mod': 'Int py_mod(Int a, Int b) { return a % b; }',
        'py_pow': '''Int py_pow(Int a, Int b) {
    if(b<0) throw domain_error("negative integer exponent requires a floating result");
    Int result=1; while(b>0) { if(b&1) result*=a; b>>=1; if(b) a*=a; } return result;
}''',
        'py_mod_wide': 'Int py_mod_wide(Int a, Int b) { return a % b; }',
        'py_bit_length': 'Int py_bit_length(Int value) { return value.bit_length(); }',
        'py_bit_count': 'Int py_bit_count(Int value) { return value.bit_count(); }',
        'py_pow_mod': '''Int py_pow_mod(Int a, Int b, Int mod) {
    if(mod==0) throw domain_error("pow() modulus cannot be zero");
    a %= mod;
    if(b<0) {
        Int r=a, next=abs(mod), x=1, y=0;
        while(next!=0) { Int q=r/next, remainder=r-q*next, z=x-q*y; r=next; next=remainder; x=y; y=z; }
        if(r!=1) throw domain_error("base is not invertible for the given modulus");
        a=x%mod; b=-b;
    }
    Int result=Int(1)%mod;
    while(b>0) { if(b&1) result=result*a%mod; b>>=1; if(b) a=a*a%mod; }
    return result;
}''',
    }
    adapted = []
    for unit in units:
        name = unit_name(unit)
        if name in replacements:
            unit = replacements[name]
        if name == 'FastScanner':
            begin = unit.index('        u64 result = 0;')
            end = unit.index('\n    }', begin)
            unit = unit[:begin] + """        u64 small = 0;
        const u64 limit = negative ? (u64(1) << 63) : u64(9223372036854775807LL);
        while (pos < data.size() && data[pos] >= '0' && data[pos] <= '9') {
            unsigned digit = data[pos++] - '0';
            if (small <= (limit - digit) / 10) { small = small * 10 + digit; continue; }
            Int result = Int(small) * 10 + digit;
            while (pos < data.size() && data[pos] >= '0' && data[pos] <= '9') result = result * 10 + (data[pos++] - '0');
            return negative ? -result : result;
        }
        if (negative && small == (u64(1) << 63)) return Int(-9223372036854775807LL - 1);
        return negative ? Int(-static_cast<long long>(small)) : Int(small);""" + unit[end:]
        if name == 'py_at' and 'template <class C>' in unit:
            const = 'const ' if 'py_at(const C&' in unit else ''
            unit = f"""template <class C>
decltype(auto) py_at({const}C& container, const Int& index) {{
    long long position = static_cast<long long>(index);
    auto size = container.size();
    if (position < 0) position += static_cast<long long>(size);
    if (position < 0 || static_cast<unsigned long long>(position) >= size) throw out_of_range("index out of range");
    return container[static_cast<size_t>(position)];
}}"""
            if not const:
                unit += '''
template <class C> requires (requires(C& container) { container.size(); container[0]; })
decltype(auto) py_at(C& container, long long index) {
    auto size = container.size();
    if (index < 0) index += static_cast<long long>(size);
    if (index < 0 || static_cast<unsigned long long>(index) >= size) throw out_of_range("index out of range");
    return container[static_cast<size_t>(index)];
}'''
        if name == 'py_truth':
            unit = unit.replace('py_truth(Int value)', 'py_truth(const Int& value)')
        if name == 'py_int':
            unit = unit.replace('return stoll(value);', 'return Int::parse(value);')
            unit = unit.replace("return value >= '0' && value <= '9' ? value - '0' : stoll(string(1, value));", 'return Int::parse(string(1, value));')
        if name == 'py_int_base':
            unit = unit.replace('stoll(value, nullptr, static_cast<int>(base))', 'Int::parse(value, static_cast<int>(base))')
        if name == 'py_find':
            unit = unit.replace('? -1 :', '? Int(-1) :')
        if name in {'py_randint', 'py_randrange','py_randbelow'}:
            unit = unit.replace('uniform_int_distribution<Int>', 'uniform_int_distribution<long long>')
        if name == 'py_bin':
            unit = '''string py_bin(Int value) {
    bool negative = value < 0;
    value = abs(value);
    string result;
    do { result.push_back((value & 1) ? '1' : '0'); value >>= 1; } while (value != 0);
    result += "b0"; if(negative) result.push_back('-'); reverse(result.begin(),result.end()); return result;
}'''
        adapted.append(unit)
    headers = '\n'.join(re.findall(r'^#include[^\n]*', runtime, re.M))
    return headers + '\nusing namespace std;\nusing u64 = unsigned long long;\n' + '\n\n'.join(adapted)


def transpile_native_full_cpp(source, filename='main.py', *, registry=None, integer_mode='adaptive'):
    if integer_mode == 'fixed64':
        from .native_fixed64 import Fixed64CppGenerator
        return Fixed64CppGenerator(source,filename,registry).generate()
    if integer_mode != 'adaptive':
        raise ValueError('integer_mode must be adaptive or fixed64')
    return NativeCppGenerator(source,filename,registry).generate()

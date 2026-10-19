from __future__ import annotations

import ast
import re
from .translator import CppGenerator, sanitize_name
from .native_full import NativeCppGenerator, replace_identifiers, prune_runtime, runtime_units, unit_name, mask_cpp
from .integer_bounds import LOW, HIGH, fits
from .fixed64_analysis import fixed_callbacks, bounded_range_accesses
from .fixed64_input import Fixed64InputMixin, INPUT_RUNTIME


class Fixed64CppGenerator(Fixed64InputMixin, NativeCppGenerator):
    """Full native syntax with fixed-width storage and no BigInt promotion."""

    integer_cpp = 'i64'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.replacements = {key: entry for key, entry in self.replacements.items() if not re.search(r'\bInt\b', mask_cpp(entry.cpp_source))}
        self.static_callbacks = fixed_callbacks(self)
        self.callback_aliases = {}
        self.guarded_indices = set()
        self.prepare_input_optimizations()

    def fixed_bounds(self, node):
        bounds=self.integer_analysis().expression(node)
        if bounds is not None:
            return bounds
        if isinstance(node,ast.BinOp):
            left,right=self.fixed_bounds(node.left),self.fixed_bounds(node.right)
            if left is None or right is None:return None
            a,b=left;c,d=right
            if isinstance(node.op,ast.Add):return a+c,b+d
            if isinstance(node.op,ast.Sub):return a-d,b-c
            if isinstance(node.op,ast.Mult):
                values=(a*c,a*d,b*c,b*d);return min(values),max(values)
            if isinstance(node.op,ast.BitAnd):
                if c>=0:return 0,d
                if a>=0:return 0,b
                return LOW,HIGH
            if isinstance(node.op,(ast.BitOr,ast.BitXor)):return LOW,HIGH
            return None
        if isinstance(node,ast.UnaryOp):
            value=self.fixed_bounds(node.operand)
            if value is None:return None
            if isinstance(node.op,ast.USub):return -value[1],-value[0]
            if isinstance(node.op,ast.UAdd):return value
        if isinstance(node,(ast.Name,ast.Attribute,ast.Subscript,ast.Call)) and self.analyzer.infer(node,self.context.env).kind in {'i64','unknown'}:
            return LOW,HIGH
        return None

    def assignment(self, node):
        if self.emit_scalar_list_unpack(node):return
        if len(node.targets)==1 and isinstance(node.targets[0],ast.Name) and node.targets[0].id in self.context.native_callable_aliases and isinstance(node.value,ast.Attribute):
            target=self.static_callbacks.get((self.context.class_name,node.value.attr))
            if target:
                entry=self.callback_aliases.setdefault(id(self.context),(self.context,{}))
                entry[1][node.targets[0].id]=target
        return super().assignment(node)

    def fixed_expression_safe(self, node):
        if not fits(self.fixed_bounds(node)):
            return False
        if isinstance(node, ast.BinOp):
            return isinstance(node.op,(ast.Add,ast.Sub,ast.Mult,ast.BitAnd,ast.BitOr,ast.BitXor,ast.Mod)) and self.fixed_expression_safe(node.left) and self.fixed_expression_safe(node.right)
        if isinstance(node, ast.UnaryOp):
            return isinstance(node.op,(ast.UAdd,ast.USub,ast.Invert)) and self.fixed_expression_safe(node.operand)
        return True

    def emit_unpack_read(self, target, value):
        if self.emit_input_unpack(target,value):return True
        if isinstance(value,ast.Call) and isinstance(value.func,ast.Name) and value.func.id=='map' and len(value.args)==2 and not value.keywords and isinstance(value.args[0],ast.Name) and value.args[0].id=='int' and self.is_input_split(value.args[1]) and not value.args[1].args and not value.args[1].keywords and self.integer_analysis().builtin('map') and self.integer_analysis().builtin('int'):
            temporary=self.context.temporary('values')
            self.emit(f'auto {temporary} = py_input.read_i64_fixed_line<{len(target.elts)}>();')
            for index,child in enumerate(target.elts):self.assign_target(child,f'{temporary}[{index}]')
            return True
        return super().emit_unpack_read(target,value)

    def native_expr(self, node):
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub) and isinstance(node.operand, ast.Constant) and isinstance(node.operand.value, int) and node.operand.value == 1 << 63:
            return '(-9223372036854775807LL - 1LL)'
        if isinstance(node, ast.Constant) and isinstance(node.value, int) and not isinstance(node.value, bool):
            if not LOW <= node.value <= HIGH:
                self.fail('integer literal is outside the fixed64 range; select adaptive integers', node)
            return f'{node.value}LL'
        if isinstance(node, ast.BinOp):
            kinds = {self.analyzer.infer(child, self.context.env).kind for child in (node.left, node.right)}
            if isinstance(node.op, ast.Div) or 'f64' in kinds:
                return NativeCppGenerator.native_expr(self, node)
            if isinstance(node.op,ast.Mod) and not (isinstance(node.left,ast.BinOp) and isinstance(node.left.op,ast.Pow)):
                bounds=self.fixed_bounds(node.left)
                if self.fixed_expression_safe(node.left) and kinds <= {'i64','unknown','bool'}:
                    left,right=self.value_expr(node.left),self.value_expr(node.right)
                    divisor=self.fixed_bounds(node.right)
                    if fits(divisor) and bounds[0]>=0 and divisor[0]>0:return f'({left} % {right})'
                    return f'py_mod({left}, {right})'
        return CppGenerator.expr(self, node)

    def call(self, node):
        target=None
        if not node.keywords:
            if isinstance(node.func,ast.Name) and node.func.id in self.context.native_callable_aliases:
                entry=self.callback_aliases.get(id(self.context))
                if entry and entry[0] is self.context:target=entry[1].get(node.func.id)
            elif isinstance(node.func,ast.Attribute) and isinstance(node.func.value,ast.Name):
                owner=self.analyzer.infer(node.func.value,self.context.env)
                if owner.kind=='class':target=self.static_callbacks.get((owner.name,node.func.attr))
        if target:
            parameters=self.program_ir.functions[target].parameter_types
            expected=list(parameters.values())
            if len(node.args)==len(expected):
                arguments=', '.join(self.call_argument(arg,expected[index]) for index,arg in enumerate(node.args))
                return f'::{sanitize_name(target)}({arguments})'
        if isinstance(node.func, ast.Attribute) and node.func.attr in {'from_bytes', 'to_bytes'}:
            return NativeCppGenerator.call(self, node)
        origin = self.analyzer.imported_origin(node.func.id) if isinstance(node.func, ast.Name) else None
        function = origin[1] if origin and origin[0] == 'math' else node.func.attr if isinstance(node.func, ast.Attribute) and isinstance(node.func.value, ast.Name) and node.func.value.id == 'math' else None
        if function in {'floor','ceil'} and len(node.args) == 1 and self.analyzer.infer(node.args[0], self.context.env).kind == 'i64':
            return self.value_expr(node.args[0])
        is_isqrt = isinstance(node.func, ast.Attribute) and isinstance(node.func.value, ast.Name) and node.func.value.id == 'math' and node.func.attr == 'isqrt'
        if (is_isqrt or origin == ('math', 'isqrt')) and len(node.args) == 1:
            return f'py_fixed_isqrt({self.value_expr(node.args[0])})'
        value = CppGenerator.call(self, node)
        if isinstance(node.func, ast.Name) and node.func.id in self.analyzer.functions and node.func.id not in self.context.declared and node.func.id not in self.context.env:
            name = sanitize_name(node.func.id)
            if value.startswith(name + '('):
                value = '::' + value
        return value

    def range_for(self, node):
        arguments = self.integer_analysis().range_arguments(node.iter)
        if arguments is None:
            return CppGenerator.range_for(self, node)
        bounds=[value if fits(value) else (LOW,HIGH) for value in arguments]
        accesses=bounded_range_accesses(self,node)
        if not accesses:
            return self.fixed_range_for(node,bounds)
        expressions=[self.expr(arg) for arg in node.iter.args]
        if len(expressions)==1:expressions.insert(0,'0LL')
        start,stop=[self.context.temporary('bounded_'+kind) for kind in ('start','stop')]
        self.emit(f'const i64 {start}={expressions[0]}, {stop}={expressions[1]};')
        positive=len(expressions)<3 or isinstance(node.iter.args[2],ast.Constant) and node.iter.args[2].value==1
        lower=f'static_cast<__int128>({start})' if positive else f'(static_cast<__int128>({stop})+1)'
        upper=f'(static_cast<__int128>({stop})-1)' if positive else f'static_cast<__int128>({start})'
        guards=[]
        for base,index,offset in accesses:
            delta=self.wide_integer_expr(offset) if offset is not None else '0LL'
            guards.append(f'({lower}+({delta})>=0 && {upper}+({delta})<static_cast<__int128>({self.expr(base)}.size()))')
        self.emit('if ('+' && '.join(guards)+') {')
        self.indent+=1
        previous=self.guarded_indices
        self.guarded_indices=previous|{(base.id,ast.dump(index)) for base,index,_ in accesses}
        self.fixed_range_for(node,bounds)
        self.guarded_indices=previous
        self.indent-=1
        self.emit('} else {')
        self.indent+=1
        self.fixed_range_for(node,bounds)
        self.indent-=1
        self.emit('}')

    def subscript(self,node,mutable):
        if isinstance(node.value,ast.Name) and (node.value.id,ast.dump(node.slice)) in self.guarded_indices:
            return f'{self.expr(node.value)}[static_cast<size_t>({self.expr(node.slice)})]'
        return super().subscript(node,mutable)

    def generate(self):
        body = self.finish_body(replace_identifiers(CppGenerator.generate(self), {'CP_NATIVE_WIDE':'__int128'}))
        runtime = fixed_runtime(self.runtime)
        runtime = prune_runtime(runtime, body)
        return runtime.replace('using u64 = unsigned long long;', 'using u64 = unsigned long long;\nusing i64 = long long;', 1) + '\n' + body


def fixed_runtime(runtime):
    replacements = {
        'py_floor_div': '''i64 py_floor_div(i64 a,i64 b) {
    if(b==0) throw domain_error("integer division by zero");
    if(a==(-9223372036854775807LL-1) && b==-1) throw overflow_error("fixed64 division overflow");
    i64 q=a/b,r=a%b; if(r!=0 && ((r<0)!=(b<0))) --q; return q;
}''',
        'py_mod': '''i64 py_mod(i64 a,i64 b) {
    if(b==0) throw domain_error("integer modulo by zero");
    i64 r=b==-1 ? 0 : a%b; if(r!=0 && ((r<0)!=(b<0))) r+=b; return r;
}''',
        'py_to_bytes': '''string py_to_bytes(i64 value,i64 length,const string& byteorder="big",bool signed_value=false) {
    if(byteorder!="big" && byteorder!="little") throw invalid_argument("invalid byteorder");
    if(length<0) throw domain_error("negative byte count");
    if(!signed_value && value<0) throw overflow_error("negative integer cannot use unsigned bytes");
    if(length==0) { if(value==0 || signed_value && value==-1) return string(); throw overflow_error("integer does not fit byte count"); }
    if(length<8) {
        unsigned bits=static_cast<unsigned>(length*8);
        if(signed_value ? value < -(1LL<<(bits-1)) || value > (1LL<<(bits-1))-1 : static_cast<u64>(value) >= (1ULL<<bits)) throw overflow_error("integer does not fit byte count");
    }
    string result(static_cast<size_t>(length),value<0 ? char(255) : char(0));
    u64 bits=static_cast<u64>(value);
    for(i64 i=0;i<min(length,8LL);++i) result[static_cast<size_t>(byteorder=="big" ? length-1-i : i)]=static_cast<char>((bits>>(i*8))&255);
    return result;
}''',
        'py_int_from_bytes': '''template<class C>
i64 py_int_from_bytes(const C& source,const string& byteorder="big",bool signed_value=false) {
    if(byteorder!="big" && byteorder!="little") throw invalid_argument("invalid byteorder");
    if(source.empty()) return 0;
    auto high=byteorder=="big" ? source.front() : source.back();
    bool negative=signed_value && (static_cast<unsigned char>(high)&128);
    u64 bits=0;
    for(size_t i=0;i<source.size();++i) {
        auto byte=static_cast<unsigned char>(source[byteorder=="big" ? source.size()-1-i : i]);
        if(i<8) bits |= static_cast<u64>(byte) << (i*8);
        else if(byte!=(negative ? 255 : 0)) throw overflow_error("integer does not fit fixed64");
    }
    if(negative && source.size()<8) bits |= ~0ULL << (source.size()*8);
    if((bits>>63)!=(negative ? 1 : 0)) throw overflow_error("integer does not fit fixed64");
    return static_cast<i64>(bits);
}''',
    }
    units = []
    for unit in runtime_units(runtime):
        unit = replacements.get(unit_name(unit),unit)
        if unit_name(unit) == 'FastScanner':
            unit = unit.replace('return negative ? -static_cast<i64>(result) : static_cast<i64>(result);', 'if(negative && result==(u64(1)<<63)) return (-9223372036854775807LL-1);\n        return negative ? -static_cast<i64>(result) : static_cast<i64>(result);')
            method='''
    template<size_t N> array<i64,N> read_i64_fixed_line() {
        array<i64,N> result{};size_t count=0;
        while(pos<data.size()) {
            while(pos<data.size() && data[pos]!='\\n' && py_space(static_cast<unsigned char>(data[pos]))) ++pos;
            if(pos==data.size()) break;
            if(data[pos]=='\\n') { ++pos;break; }
            if(count==N) throw runtime_error("too many values to unpack");
            result[count++]=read_i64_token();
        }
        if(count!=N) throw runtime_error("not enough values to unpack");
        return result;
    }
'''
            unit=unit[:-2]+method+INPUT_RUNTIME+'};'
        units.append(unit)
    headers = '\n'.join(re.findall(r'^#include[^\n]*', runtime, re.M))
    return headers + '\n#include <array>\nusing namespace std;\nusing i64 = long long;\nusing u64 = unsigned long long;\n' + '\n\n'.join(units) + '''
double py_float_mod(double a,double b) {
    if(b==0) throw domain_error("float modulo by zero");
    double r=fmod(a,b); if(r!=0 && ((r<0)!=(b<0))) r+=b; return r;
}
i64 py_fixed_isqrt(i64 value) {
    if(value<0) throw domain_error("isqrt() argument must be nonnegative");
    u64 root=static_cast<u64>(sqrt(static_cast<long double>(value)));
    while(static_cast<__int128>(root)*root>value) --root;
    while(static_cast<__int128>(root+1)*(root+1)<=value) ++root;
    return static_cast<i64>(root);
}
'''

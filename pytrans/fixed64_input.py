from __future__ import annotations

import ast
from .translator import I64, list_type


INPUT_RUNTIME = r'''
    template<size_t N,class F> array<i64,N> read_i64_fixed_line_map(F function) {
        auto result=read_i64_fixed_line<N>();
        for(auto& value:result) value=function(value);
        return result;
    }
    template<size_t P,size_t S,class F>
    array<i64,P+S> read_i64_star_line(vector<i64>& rest,F function) {
        array<i64,P+S> edges{};
        rest.clear();size_t count=0;
        if(rest.capacity()==0) {
            size_t end=data.find('\n',pos);
            if(end==string::npos) end=data.size();
            rest.reserve(min<size_t>((end-pos+1)/2,1<<20));
        }
        while(pos<data.size()) {
            while(pos<data.size() && data[pos]!='\n' && py_space(static_cast<unsigned char>(data[pos]))) ++pos;
            if(pos==data.size()) break;
            if(data[pos]=='\n') { ++pos;break; }
            i64 value=function(read_i64_token());
            if(count<P) edges[count]=value;
            else rest.push_back(value);
            ++count;
        }
        if(count<P+S) throw runtime_error("not enough values to unpack");
        for(size_t i=0;i<S;++i) edges[P+i]=rest[rest.size()-S+i];
        rest.resize(rest.size()-S);
        return edges;
    }
'''


class Fixed64InputMixin:
    """Specialize pure integer-line readers and preserve observable list ownership."""

    def prepare_input_optimizations(self):
        self.input_parents={child:parent for parent in ast.walk(self.tree) for child in ast.iter_child_nodes(parent)}
        self.input_readers={}
        self.input_global_writes=set()
        self.input_reflection=any(isinstance(node,ast.Call) and isinstance(node.func,ast.Name) and node.func.id in {'locals','globals','vars','eval','exec'} for node in ast.walk(self.tree))
        for statement in self.tree.body:
            if not isinstance(statement,(ast.FunctionDef,ast.ClassDef)):
                self.input_global_writes.update(child.id for child in ast.walk(statement) if isinstance(child,ast.Name) and isinstance(child.ctx,ast.Store))
        for node in ast.walk(self.tree):
            if isinstance(node,ast.Global):self.input_global_writes.update(node.names)
        if self.input_reflection:return
        for name,function in self.analyzer.functions.items():
            if name in self.input_global_writes or self.analyzer.hoisted_captures.get(name) or ('function',name) in self.replacements:continue
            args=function.args
            if function.decorator_list or args.vararg or args.kwarg or args.kwonlyargs:continue
            body=[statement for statement in function.body if not (isinstance(statement,ast.Expr) and isinstance(statement.value,ast.Constant) and isinstance(statement.value.value,str))]
            if len(body)!=1 or not isinstance(body[0],ast.Return):continue
            parameters=[arg.arg for arg in args.posonlyargs+args.args]
            value=body[0].value
            expression=None
            if isinstance(value,ast.ListComp) and len(value.generators)==1:
                loop=value.generators[0]
                if loop.ifs or loop.is_async or not isinstance(loop.target,ast.Name) or not self.plain_integer_input_split(loop.iter,parameters):continue
                variable=loop.target.id;expression=value.elt
            elif isinstance(value,ast.Call) and isinstance(value.func,ast.Name) and value.func.id=='list' and len(value.args)==1 and not value.keywords and self.direct_integer_map(value.args[0],parameters):
                if 'list' in parameters or 'list' in self.input_global_writes or 'list' in self.analyzer.functions:continue
                variable='__input_value';expression=ast.Call(func=ast.Name(id='int',ctx=ast.Load()),args=[ast.Name(id=variable,ctx=ast.Load())],keywords=[])
            else:continue
            def pure(node):
                if isinstance(node,ast.Constant):return isinstance(node.value,int) and not isinstance(node.value,bool)
                if isinstance(node,ast.Name):return node.id in parameters and node.id!=variable
                if isinstance(node,ast.BinOp) and isinstance(node.op,(ast.Add,ast.Sub)):return pure(node.left) and pure(node.right)
                if isinstance(node,ast.UnaryOp) and isinstance(node.op,(ast.UAdd,ast.USub)):return pure(node.operand)
                return isinstance(node,ast.Call) and isinstance(node.func,ast.Name) and node.func.id=='int' and len(node.args)==1 and not node.keywords and isinstance(node.args[0],ast.Name) and node.args[0].id==variable
            converted=any(isinstance(node,ast.Call) and isinstance(node.func,ast.Name) and node.func.id=='int' for node in ast.walk(expression))
            if not converted or variable=='int' or not pure(expression) or 'int' in parameters or 'int' in self.input_global_writes or 'int' in self.analyzer.functions:continue
            defaults={}
            valid=True
            for parameter,default in zip(parameters[len(parameters)-len(args.defaults):],args.defaults):
                try:value=ast.literal_eval(default)
                except (ValueError,TypeError):valid=False;break
                if not isinstance(value,int) or isinstance(value,bool):valid=False;break
                defaults[parameter]=ast.Constant(value=value)
            if valid:self.input_readers[name]=(parameters,defaults,expression,variable)

    def plain_integer_input_split(self,node,parameters=()):
        if not self.is_input_split(node) or node.args or node.keywords:return False
        reader=node.func.value.func
        loaded={child.id for child in ast.walk(reader) if isinstance(child,ast.Name)}
        if loaded&set(parameters) or self.context is not None and loaded&self.context.env.keys():return False
        if isinstance(reader,ast.Name):
            if reader.id in parameters or reader.id in self.analyzer.functions:return False
            for statement in ast.walk(self.tree):
                if isinstance(statement,(ast.Assign,ast.AnnAssign)):
                    targets=statement.targets if isinstance(statement,ast.Assign) else [statement.target]
                    if any(isinstance(target,ast.Name) and target.id==reader.id for target in targets):
                        if not isinstance(statement,ast.Assign) or not self.analyzer.is_input_alias_assignment(statement):return False
        return True

    def direct_integer_map(self,node,parameters=()):
        return isinstance(node,ast.Call) and isinstance(node.func,ast.Name) and node.func.id=='map' and len(node.args)==2 and not node.keywords and isinstance(node.args[0],ast.Name) and node.args[0].id=='int' and not {'map','int'}&set(parameters) and not {'map','int'}&self.input_global_writes and not {'map','int'}&self.analyzer.functions.keys() and self.plain_integer_input_split(node.args[1],parameters)

    def integer_line_mapper(self,node):
        if self.input_reflection:return None
        if self.direct_integer_map(node) and self.integer_analysis().builtin('map') and self.integer_analysis().builtin('int'):
            return '[](i64 input_value) { return input_value; }'
        if not isinstance(node,ast.Call) or not isinstance(node.func,ast.Name) or node.func.id not in self.input_readers or node.keywords or node.func.id in self.context.env:return None
        parameters,defaults,expression,variable=self.input_readers[node.func.id]
        if len(node.args)>len(parameters) or any(isinstance(arg,ast.Starred) for arg in node.args):return None
        values={}
        for index,name in enumerate(parameters):
            value=node.args[index] if index<len(node.args) else defaults.get(name)
            if value is None or self.analyzer.infer(value,self.context.env).kind not in {'i64','bool'}:return None
            values[name]=value
        bindings={}
        for name,value in values.items():
            temporary=self.context.temporary('input_argument')
            self.emit(f'const i64 {temporary} = {self.expr(value)};')
            bindings[name]=temporary
        def render(value):
            if isinstance(value,ast.Call):return 'input_value'
            if isinstance(value,ast.Name):return bindings[value.id]
            if isinstance(value,ast.Constant):return self.expr(value)
            if isinstance(value,ast.UnaryOp):return f'({self.operator(value.op)}{render(value.operand)})'
            return f'({render(value.left)} {self.operator(value.op)} {render(value.right)})'
        return '[&](i64 input_value) { return '+render(expression)+'; }'

    def input_scalar_targets(self,targets,source=None):
        return all(isinstance(child,ast.Name) and child.id!=source and child.id not in self.context.reference_expressions and self.analyzer.infer(child,self.context.env).kind in {'i64','unknown','bool'} for child in targets)

    def emit_input_unpack(self,target,value):
        if not self.input_scalar_targets(target.elts) or len(target.elts)>64:return False
        if self.direct_integer_map(value):return False
        mapper=self.integer_line_mapper(value)
        if mapper is None:return False
        temporary=self.context.temporary('input_values')
        self.emit(f'auto {temporary} = py_input.read_i64_fixed_line_map<{len(target.elts)}>({mapper});')
        for index,child in enumerate(target.elts):self.assign_unpack_target(child,f'{temporary}[{index}]',I64)
        return True

    def emit_scalar_list_unpack(self,node):
        if self.input_reflection:return False
        if len(node.targets)!=1 or not isinstance(node.targets[0],(ast.Tuple,ast.List)) or not isinstance(node.value,ast.Name):return False
        target=node.targets[0];kind=self.analyzer.infer(node.value,self.context.env)
        if kind.kind!='list' or kind.args[0].kind!='i64' or not self.input_scalar_targets(target.elts,node.value.id):return False
        temporary=self.context.temporary('unpack_source')
        self.emit(f'const auto& {temporary} = {self.expr(node.value)};')
        self.emit(f'if ({temporary}.size() != {len(target.elts)}) throw runtime_error("wrong number of values to unpack");')
        for index,child in enumerate(target.elts):self.assign_unpack_target(child,f'{temporary}[{index}]',I64)
        return True

    def reusable_input_rest(self,name,target):
        if self.input_reflection or any(name in captures for captures in self.analyzer.hoisted_captures.values()):return False
        scope=self.input_parents.get(target)
        while scope is not None and not isinstance(scope,ast.FunctionDef):scope=self.input_parents.get(scope)
        if scope is None or name in {arg.arg for arg in scope.args.posonlyargs+scope.args.args+scope.args.kwonlyargs}:return False
        if any(isinstance(node,(ast.Try,ast.Global,ast.Nonlocal,ast.Lambda)) or isinstance(node,ast.FunctionDef) and node is not scope for node in ast.walk(scope)):return False
        uses=[node for node in ast.walk(scope) if isinstance(node,ast.Name) and node.id==name]
        if sum(isinstance(node.ctx,ast.Store) for node in uses)!=1 or not any(isinstance(node.ctx,ast.Load) for node in uses):return False
        for node in uses:
            if not isinstance(node.ctx,ast.Load):continue
            parent=self.input_parents.get(node)
            if not isinstance(parent,ast.Assign) or parent.value is not node or len(parent.targets)!=1 or not isinstance(parent.targets[0],(ast.Tuple,ast.List)) or not self.input_scalar_targets(parent.targets[0].elts,name):return False
        return True

    def emit_starred_unpack(self,target,value):
        stars=[index for index,child in enumerate(target.elts) if isinstance(child,ast.Starred)]
        if len(stars)!=1:return super().emit_starred_unpack(target,value)
        index=stars[0];star=target.elts[index].value;suffix=len(target.elts)-index-1
        scalars=target.elts[:index]+target.elts[index+1:]
        if not isinstance(star,ast.Name) or not self.input_scalar_targets(scalars,star.id) or len({child.id for child in scalars})!=len(scalars) or len(scalars)>64:return super().emit_starred_unpack(target,value)
        kind=self.analyzer.infer(star,self.context.env)
        if kind.kind!='list' or kind.args[0].kind!='i64':return super().emit_starred_unpack(target,value)
        mapper=self.integer_line_mapper(value)
        if mapper is None:return super().emit_starred_unpack(target,value)
        reuse=self.reusable_input_rest(star.id,target)
        rest=self.expr(star) if reuse else self.context.temporary('input_rest')
        if not reuse:self.emit(f'vector<i64> {rest};')
        edges=self.context.temporary('input_edges')
        self.emit(f'auto {edges} = py_input.read_i64_star_line<{index},{suffix}>({rest}, {mapper});')
        for position,child in enumerate(target.elts[:index]):self.assign_unpack_target(child,f'{edges}[{position}]',I64)
        if not reuse:self.assign_unpack_target(star,f'std::move({rest})',list_type(I64))
        for position,child in enumerate(target.elts[index+1:]):self.assign_unpack_target(child,f'{edges}[{index+position}]',I64)
        return True

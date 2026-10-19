from __future__ import annotations

import ast


def bounded_range_accesses(generator,loop):
    """Plan guarded direct accesses for scalar-only loops without mutation hazards."""
    if not isinstance(loop.target,ast.Name) or loop.orelse:
        return []
    call=loop.iter
    if len(call.args)==3:
        step=call.args[2]
        if not (isinstance(step,ast.Constant) and step.value==1 or isinstance(step,ast.UnaryOp) and isinstance(step.op,ast.USub) and isinstance(step.operand,ast.Constant) and step.operand.value==1):return []
    if any(not isinstance(statement,(ast.Assign,ast.AugAssign)) for statement in loop.body):return []
    analysis=generator.integer_analysis()
    for node in ast.walk(loop):
        if isinstance(node,(ast.NamedExpr,ast.Starred,ast.Lambda,ast.ListComp,ast.GeneratorExp)):return []
        if isinstance(node,ast.Call) and node is not call:
            if node.keywords or not isinstance(node.func,ast.Name) or node.func.id not in {'min','max','len','abs','int'} or not analysis.builtin(node.func.id):return []
            if any(generator.analyzer.infer(arg,generator.context.env).kind in {'class','function','optional'} for arg in node.args):return []
        if isinstance(node,ast.Attribute):return []
    body_nodes=[node for statement in loop.body for node in ast.walk(statement)]
    stores={node.id for node in body_nodes if isinstance(node,ast.Name) and isinstance(node.ctx,ast.Store)}
    if loop.target.id in stores:return []
    result={}
    for node in body_nodes:
        if not isinstance(node,ast.Subscript):continue
        if not isinstance(node.value,ast.Name) or node.value.id in stores:return []
        kind=generator.analyzer.infer(node.value,generator.context.env)
        if kind.kind!='list' or kind.args[0].kind not in {'i64','f64','bool'}:return []
        index=node.slice;offset=None
        if isinstance(index,ast.Name) and index.id==loop.target.id:pass
        elif isinstance(index,ast.BinOp) and isinstance(index.left,ast.Name) and index.left.id==loop.target.id and isinstance(index.op,(ast.Add,ast.Sub)) and isinstance(index.right,(ast.Name,ast.Constant)):
            if isinstance(index.right,ast.Name) and (index.right.id in stores or index.right.id==loop.target.id):return []
            if generator.analyzer.infer(index.right,generator.context.env).kind not in {'i64','bool'}:return []
            offset=index.right if isinstance(index.op,ast.Add) else ast.UnaryOp(op=ast.USub(),operand=index.right)
        else:return []
        result[(node.value.id,ast.dump(index))]=(node.value,index,offset)
    bound_names={node.id for arg in call.args for node in ast.walk(arg) if isinstance(node,ast.Name) and isinstance(node.ctx,ast.Load)}
    if bound_names&stores:return []
    return list(result.values())


def fixed_callbacks(generator):
    """Bind immutable constructor callbacks only when every call has one identity."""
    if generator.replacements:
        return {}
    if any(isinstance(node,ast.Call) and isinstance(node.func,ast.Attribute) and node.func.attr in {'__new__','__init__'} for node in ast.walk(generator.tree)):
        return {}
    analyzer=generator.analyzer
    parents={child:parent for parent in ast.walk(generator.tree) for child in ast.iter_child_nodes(parent)}
    global_stores=set()
    def module_stores(nodes):
        for node in nodes:
            if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef,ast.ClassDef)):
                continue
            global_stores.update(child.id for child in ast.walk(node) if isinstance(child,ast.Name) and isinstance(child.ctx,ast.Store))
    module_stores(generator.tree.body)
    for scope in ast.walk(generator.tree):
        if isinstance(scope,(ast.FunctionDef,ast.AsyncFunctionDef)):
            global_stores.update(name for statement in ast.walk(scope) if isinstance(statement,ast.Global) for name in statement.names)
    def global_function(node,call):
        if not isinstance(node,ast.Name) or node.id not in analyzer.functions or node.id in global_stores or analyzer.hoisted_captures.get(node.id):return None
        scope=parents.get(call)
        while scope is not None and not isinstance(scope,(ast.FunctionDef,ast.AsyncFunctionDef)):
            scope=parents.get(scope)
        if scope is not None:
            bound=analyzer.bound_names(scope.body)|{arg.arg for arg in scope.args.posonlyargs+scope.args.args}
            bound-=set(name for statement in scope.body if isinstance(statement,ast.Global) for name in statement.names)
            if node.id in bound:return None
        if generator.program_ir.functions[node.id].mutable_parameters:return None
        return node.id
    result={}
    for name,attributes in generator.optimization_planner.immutable_callables.items():
        if analyzer.class_bases.get(name) or any(name in bases for bases in analyzer.class_bases.values()):continue
        initializer=analyzer.class_methods(name).get('__init__')
        if initializer is None or any(isinstance(child,(ast.Return,ast.Yield,ast.YieldFrom)) for child in ast.walk(initializer)):continue
        arguments=[arg.arg for arg in initializer.args.posonlyargs+initializer.args.args if arg.arg!='self']
        defaults=dict(zip(arguments[len(arguments)-len(initializer.args.defaults):],initializer.args.defaults))
        constructors=[node for node in ast.walk(generator.tree) if isinstance(node,ast.Call) and isinstance(node.func,ast.Name) and node.func.id==name]
        if not constructors:continue
        for attribute in attributes:
            if any(isinstance(child,ast.Call) and isinstance(child.func,ast.Attribute) and child.func.attr==attribute for child in ast.walk(initializer)):continue
            assignments=[node for node in initializer.body if isinstance(node,ast.Assign) and len(node.targets)==1 and isinstance(node.targets[0],ast.Attribute) and isinstance(node.targets[0].value,ast.Name) and node.targets[0].value.id=='self' and node.targets[0].attr==attribute]
            all_writes=[node for node in ast.walk(initializer) if isinstance(node,ast.Attribute) and isinstance(node.ctx,ast.Store) and node.attr==attribute]
            if len(assignments)!=1 or len(all_writes)!=1:continue
            value=assignments[0].value
            if not isinstance(value,ast.Name) or value.id not in arguments:continue
            if any(isinstance(child,ast.Name) and isinstance(child.ctx,ast.Store) and child.id==value.id for child in ast.walk(initializer)):continue
            position=arguments.index(value.id);identities=[]
            for call in constructors:
                if any(isinstance(arg,ast.Starred) for arg in call.args) or any(keyword.arg is None for keyword in call.keywords):identities.append(None);continue
                argument=call.args[position] if position<len(call.args) else next((keyword.value for keyword in call.keywords if keyword.arg==value.id),defaults.get(value.id))
                identities.append(global_function(argument,call))
            if identities and identities[0] is not None and all(identity==identities[0] for identity in identities):result[(name,attribute)]=identities[0]
    return result

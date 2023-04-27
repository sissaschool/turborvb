import csnake as cs

def add_this(cw, d):

    """
    this function adds wrapper for one cudalib API call

    it is a mess and needs a refactor, but for time being its fine
    """

    args = []
    offload_args = []
    called_args = []
    return_value = ""

    for x in d["args"]:
        qualifiers = ["const"]
        if "const" in x:
            if not x["const"]:
                qualifiers.remove("const")
        name = x["name"]
        type_ = x["type"] + " *"

        v = cs.Variable(name, type_, qualifiers = qualifiers)
        entry_offload = {}
        entry_offload["type"] = v.primitive
        entry_offload["name"] = v.name

        if "recast_offload" in x:
            entry_offload["type"] = x["recast_offload"]

        args.append(v)
        if x["offload"]:
            offload_args.append(entry_offload)

        try:
            if x["return"]:
                return_value = x["name"]
        except:
            call_name = x["name"]
            if x["offload"]:
                call_name = call_name + "_"
            elif "recast" in x:
                call_name = "* " + call_name + "__"
            else:
                call_name = "* " + call_name
            if "norefpass" in x:
                if x["norefpass"]:
                    call_name = call_name[2:]
            called_args.append(call_name)


    function_name = f"{d['library'].lower()}_{d['name']}_offload_"
    if "suffix" in d:
        if not d["suffix"]:
            function_name = function_name[:-8]
    fun = cs.Function(function_name,
                      "void",
                      arguments = args)

    n_name = d["library"].lower() + d["name"][0].upper() + d["name"][1:].lower()
    if "new_name" in d:
        n_name = d["new_name"]

    for x in d["args"]:
        try:
            new_type = x["recast"]["to"]
            fun.add_code((f"{new_type} * {x['name']}__ = ({new_type} *) *{x['name']};",))
        except:
            pass
    if len(offload_args) > 0:
        fun.add_code((f"#pragma omp target data use_device_ptr({', '.join([ x['name'] for x in offload_args ])})",))
    fun.add_code(("{",))
    for oa in offload_args:
        fun.add_code((f"{oa['type']} * {oa['name']}_ = ({oa['type']} *) {oa['name']};",))
    if "versions" in d:
        for v in d["versions"]:
            called_args_ = list(called_args)
            for vv in v:
                called_args_ = [ vv["pass"] if x == f"* {vv['variable']}" else x for x in called_args_ ]
                called_args_string = ", ".join(called_args_)
            condition = [ f"(toupper(*{vv['variable']}) == *\"{vv['value']}\")" for vv in v ]
            condition = " && ".join(condition)
            fun.add_code(f"if ({condition}) " + "{")
            fun.add_code((f"{('*'+return_value+' = ') if return_value != '' else ''}{n_name}({called_args_string});", ))
            fun.add_code("}")
    else:
        called_args_string = ", ".join(called_args)
        fun.add_code((f"{('*'+return_value+' = ') if return_value != '' else ''}{n_name}({called_args_string});", ))
    fun.add_code(("}",))

    cw.start_if_def(f"_{d['library'].upper()}")
    cw.add_function_definition(fun)
    cw.end_if_def()

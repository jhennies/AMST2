
def set_resources(
        sn_args, rule_names, mem_defaults, runtime_defaults,
        mem_args=None, runtime_args=None
):

    assert len(rule_names) == len(mem_defaults) == len(runtime_defaults)

    def _split_args(args):
        result_dict = dict()
        if args is not None:
            for entry in args:
                rule_name, value = str.split(entry, ':')
                result_dict[rule_name] = int(value)
        return result_dict

    mem_dict = _split_args(mem_args)
    runtime_dict = _split_args(runtime_args)

    sn_args.set_resources = {
        rule_name: dict(
            mem_mb=mem_defaults[idx] if rule_name not in mem_dict else mem_dict[rule_name],
            runtime=runtime_defaults[idx] if rule_name not in runtime_dict else runtime_dict[rule_name]
        )
        for idx, rule_name in enumerate(rule_names)
    }

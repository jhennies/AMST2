
def get_cluster_settings(sn_args, config_filepath):

    import json
    with open(config_filepath, mode='r') as f:
        config = json.load(f)

    import grp
    import os
    group_name = grp.getgrgid(os.getgid()).gr_name

    sn_args.executor = config['executor']

    from snakemake.resources import DefaultResources
    default_resources = DefaultResources()
    default_resources.set_resource('partition', 'htc-el8')
    default_resources.set_resource('mem_mb', 1000)
    default_resources.set_resource('runtime', 10)
    default_resources.set_resource('slurm_account', group_name)
    default_resources.set_resource('name', '{rule}-{wildcards}')
    sn_args.default_resources = default_resources

    sn_args.restart_times = 1
    sn_args.max_jobs_per_second = 10
    sn_args.max_status_checks_per_second = 1
    sn_args.local_cores = 1
    sn_args.latency_wait = 3
    sn_args.jobs = 500
    sn_args.keep_going = True
    sn_args.rerun_incomplete = True
    sn_args.printshellcmds = True
    sn_args.scheduler = 'greedy'
    sn_args.use_conda = True
    sn_args.verbose = False
    sn_args.quiet = False

    return sn_args

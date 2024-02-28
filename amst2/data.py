
import os
import numpy as np


def snk_stack_to_ome_zarr():

    # ----------------------------------------------------
    import argparse

    parser = argparse.ArgumentParser(
        description='Convert a dataset within a h5 container or a tif stack to ome.zarr',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('stack_path', type=str,
                        help='Input h5 container or tif stack')
    parser.add_argument('target_dirpath', type=str,
                        help='Output director which will contain the cache directory and, by default, also the '
                             'ome.zarr output')
    parser.add_argument('--ome_zarr_filename', type=str, default=None,
                        help='Filename of the ome.zarr output')
    parser.add_argument('--ome_zarr_filepath', type=str, default=None,
                        help='Full filepath for the ome.zarr output. Over-writes --ome_zarr_filename')
    parser.add_argument('--stack_pattern', type=str, default='*.tif',
                        help='File pattern for globbing the input stack; default="*.tif"')
    parser.add_argument('--stack_key', type=str, default='data',
                        help='Path within input h5 file; default="data"')
    parser.add_argument('--resolution', type=float, nargs=3, default=(1., 1., 1.),
                        help='Resolution of input data; default=(1., 1., 1.)')
    parser.add_argument('--unit', type=str, default='pixel',
                        help='Unit of input resolution; default="pixel"')
    parser.add_argument('--downsample_type', type=str, default='Average',
                        help='Downsample type used to create the resolution pyramid; default="Average"')
    parser.add_argument('--downsample_factors', type=int, nargs='+', default=(2, 2, 2),
                        help='Downsample factors used to create the resolution pyramid; default=(2, 2, 2)')
    parser.add_argument('--name', type=str, default=None,
                        help='Name of the dataset; defaults to the filename (without ome.zarr extension)')
    parser.add_argument('--chunk_size', type=int, nargs=3, default=[1, 512, 512],
                        help='Chunk size of the ome-zarr dataset; default=[1, 512, 512]')
    parser.add_argument('--save_bounds', action='store_true',
                        help='Saves a json file alongside that contains the bounds of the non-zero area of each slice')
    parser.add_argument('--cores', type=int, default=1,
                        help='Maximum number of available cores')
    parser.add_argument('--batch_size', type=int, default=16,
                        help='Number of slices forming one batch; default=16')
    parser.add_argument('--continue_run', action='store_true',
                        help='Needed to continue a run')
    parser.add_argument('--cluster', type=str, default=None,
                        help='Use "--cluster slurm" to run on a slurm cluster (only tested on EMBL resources)')
    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()
    stack_path = args.stack_path
    target_dirpath = args.target_dirpath
    ome_zarr_filename = args.ome_zarr_filename
    ome_zarr_filepath = args.ome_zarr_filepath
    stack_pattern = args.stack_pattern
    stack_key = args.stack_key
    resolution = args.resolution
    unit = args.unit
    downsample_type = args.downsample_type
    downsample_factors = args.downsample_factors
    name = args.name
    chunk_size = args.chunk_size
    save_bounds = args.save_bounds
    cores = args.cores
    batch_size = args.batch_size
    continue_run = args.continue_run
    cluster = args.cluster
    verbose = args.verbose

    # import time
    import json
    from amst2.library.data import estimate_mem_mb

    if ome_zarr_filepath is None:
        if ome_zarr_filename is not None:
            ome_zarr_filepath = os.path.join(target_dirpath, ome_zarr_filename)
        else:
            ome_zarr_filepath = os.path.join(
                target_dirpath,
                os.path.split(stack_path)[1].replace('.h5', '') + '.ome.zarr'
            )
    if ome_zarr_filename is None:
        ome_zarr_filename = os.path.split(ome_zarr_filepath)[1]

    target_dirpath = os.path.abspath(target_dirpath)
    ome_zarr_filepath = os.path.abspath(ome_zarr_filepath)
    stack_path = os.path.abspath(stack_path)

    if not os.path.exists(target_dirpath):
        os.mkdir(target_dirpath)
    cache_dirpath = os.path.join(target_dirpath, 'snk_cache')
    if not os.path.exists(cache_dirpath):
        os.mkdir(cache_dirpath)
    this_cache_dirpath = os.path.join(cache_dirpath, f'stack_to_ome_zarr_{ome_zarr_filename.replace(".ome.zarr", "")}')
    if not os.path.exists(this_cache_dirpath):
        os.mkdir(this_cache_dirpath)
    elif os.path.exists(this_cache_dirpath) and not continue_run:
        raise RuntimeError('Cache directory exists. If you want to continue use --continue_run')
    log_dirpath = os.path.join(this_cache_dirpath, 'logs')
    if not os.path.exists(log_dirpath):
        os.mkdir(log_dirpath)

    from squirrel.library.io import load_data_handle
    data_h, shape_h = load_data_handle(stack_path, key=stack_key, pattern=stack_pattern)
    batch_ids = [x for x in range(0, shape_h[0], batch_size)]

    src_dirpath = os.path.dirname(os.path.realpath(__file__))

    run_info = dict(
        stack_path=stack_path,
        target_dirpath=target_dirpath,
        ome_zarr_filename=ome_zarr_filename,
        ome_zarr_filepath=ome_zarr_filepath,
        stack_pattern=stack_pattern,
        stack_key=stack_key,
        resolution=resolution,
        unit=unit,
        downsample_type=downsample_type,
        downsample_factors=downsample_factors,
        name=name,
        chunk_size=chunk_size,
        save_bounds=save_bounds,
        cores=cores,
        verbose=verbose,
        batch_size=batch_size,
        batch_ids=batch_ids,
        stack_shape=shape_h,
        continue_run=continue_run,
        src_dirpath=src_dirpath
    )

    with open(os.path.join(this_cache_dirpath, 'run.json'), mode='w') as f:
        json.dump(run_info, f, indent=2)

    from snakemake.cli import parse_args, args_to_api
    from pathlib import Path
    parser, sn_args = parse_args({})
    sn_args.dryrun = False
    sn_args.unlock = False
    sn_args.directory = Path(this_cache_dirpath)
    sn_args.snakefile = Path(os.path.join(src_dirpath, 'snakemake_workflows/stack_to_ome_zarr.snk'))
    sn_args.cores = cores
    sn_args.set_threads = dict(batch_to_ome_zarr=min(batch_size, cores))
    if cluster is not None:
        sn_args.set_resources = dict(
            batch_to_ome_zarr=dict(
                # I'm purposely not using snakemake's functionality here to determine mem_mb on-the-fly
                mem_mb=int(np.ceil(estimate_mem_mb(data_h) * batch_size * 8)),
                time=10,
                partition='htc-el8'
            )
            # default=dict(
            #     partition='htc-el8',
            #     mem_mb=1000,
            #     time_min=10
            # )
        )

        # sn_args.profile = 'cluster_slurm'  # os.path.join(src_dirpath, 'cluster', cluster, 'config.yaml')
        sn_args.executor = cluster
        # sn_args.resources = dict(
        #     partition='htc-el8',
        #
        # )
        # sn_args.configfile = Path(os.path.join(src_dirpath, 'cluster', cluster, 'config.yaml'))

        jobscript_filepath = os.path.join(this_cache_dirpath, 'jobscript.sh')
        with open(jobscript_filepath, mode='w') as f:
            f.write(
                f'mkdir -p {os.path.join(this_cache_dirpath, "logs", "{rule}")} && '
                'sbatch '
                '--partition={resources.partition} '
                '--qos={resources.qos} '
                '--cpus-per-task={threads} '
                '--mem={resources.mem_mb} '
                '--job-name=smk-{rule}-{wildcards} '
                f'--output={os.path.join(log_dirpath, "{rule}-{wildcards}-%j.out")} '
                '--account=mattei'
            )
        sn_args.jobscript = jobscript_filepath
        # sn_args.jobscript = f'mkdir -p {os.path.join(this_cache_dirpath, "logs", "{rule}")} && ' \
        #                     'sbatch ' \
        #                     '--partition={resources.partition} ' \
        #                     '--qos={resources.qos} ' \
        #                     '--cpus-per-task={threads} ' \
        #                     '--mem={resources.mem_mb} ' \
        #                     '--job-name=smk-{rule}-{wildcards} ' \
        #                     f'--output={os.path.join(log_dirpath, "{rule}-{wildcards}-%j.out")}'
        from snakemake.resources import DefaultResources
        default_resources = DefaultResources()
        default_resources.set_resource('partition', 'htc-el8')
        default_resources.set_resource('mem_mb', 1000)
        default_resources.set_resource('time', 10)
        sn_args.default_resources = default_resources
        # sn_args.default_resources = dict(
        #     partition='htc-el8',
        #     mem_mb=1000,
        #     time_min=10
        # )
        sn_args.restart_times = 1
        sn_args.max_jobs_per_second = 10
        sn_args.max_status_checks_per_second = 1
        sn_args.local_cores = 1
        sn_args.latency_wait = 60
        sn_args.jobs = 500
        sn_args.keep_going = True
        sn_args.rerun_incomplete = True
        sn_args.printshellcmds = True
        sn_args.scheduler = 'greedy'
        sn_args.use_conda = True
    args_to_api(sn_args, parser)


if __name__ == '__main__':
    snk_stack_to_ome_zarr()

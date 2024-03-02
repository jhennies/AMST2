
def add_common_arguments_to_parser(parser):

    parser.add_argument('--cores', type=int, default=1,
                        help='Maximum number of available cores')
    parser.add_argument('--batch_size', type=int, default=16,
                        help='Number of slices forming one batch; default=16')
    parser.add_argument('--continue_run', action='store_true',
                        help='Needed to continue a run')
    parser.add_argument('--cluster', type=str, default=None,
                        help='Use "--cluster slurm" to run on a slurm cluster (only tested on EMBL resources)')
    parser.add_argument('-v', '--verbose', action='store_true')


def common_args_to_dict(args):

    return dict(
        cores=args.cores,
        batch_size=args.batch_size,
        continue_run=args.continue_run,
        cluster=args.cluster,
        verbose=args.verbose
    )


def add_ome_zarr_arguments_to_parser(parser):

    parser.add_argument('--resolution', type=float, nargs=3, default=(1., 1., 1.),
                        help='Resolution of input data; default=(1., 1., 1.)')
    parser.add_argument('--unit', type=str, default='pixel',
                        help='Unit of input resolution; default="pixel"')
    parser.add_argument('--downsample_type', type=str, default='Average',
                        help='Downsample type used to create the resolution pyramid; default="Average"')
    parser.add_argument('--downsample_factors', type=int, nargs='+', default=(2, 2, 2),
                        help='Downsample factors used to create the resolution pyramid; default=(2, 2, 2)')
    parser.add_argument('--chunk_size', type=int, nargs=3, default=[1, 512, 512],
                        help='Chunk size of the ome-zarr dataset; default=[1, 512, 512]')
    parser.add_argument('--name', type=str, default=None,
                        help='Name of the dataset; defaults to the filename (without ome.zarr extension)')


def ome_zarr_args_to_dict(args):

    return dict(
        resolution=args.resolution,
        unit=args.unit,
        downsample_type=args.downsample_type,
        downsample_factors=args.downsample_factors,
        name=args.name,
        chunk_size=args.chunk_size
    )


def add_snakemake_arguments_to_parser(parser):

    parser.add_argument('--dryrun', action='store_true',
                        help='Snakemake argument')
    parser.add_argument('--unlock', action='store_true',
                        help='Snakemake argument')


def args_to_snakemake_arguments(args, sn_args):

    sn_args.dryrun = args.dryrun
    sn_args.unlock = args.unlock
    sn_args.cores = args.cores


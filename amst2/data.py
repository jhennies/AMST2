
import os


def snk_stack_to_ome_zarr():

    # ----------------------------------------------------
    import argparse
    from amst2.library.argument_parsing import (
        add_common_arguments_to_parser, common_args_to_dict,
        add_ome_zarr_arguments_to_parser, ome_zarr_args_to_dict,
        add_snakemake_arguments_to_parser, args_to_snakemake_arguments
    )

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
    parser.add_argument('--save_bounds', action='store_true',
                        help='Saves a json file alongside that contains the bounds of the non-zero area of each slice')

    add_common_arguments_to_parser(parser)
    add_ome_zarr_arguments_to_parser(parser)
    add_snakemake_arguments_to_parser(parser)

    args = parser.parse_args()

    stack_path = os.path.abspath(args.stack_path)
    target_dirpath = args.target_dirpath
    ome_zarr_filename = args.ome_zarr_filename
    ome_zarr_filepath = args.ome_zarr_filepath
    stack_pattern = args.stack_pattern
    stack_key = args.stack_key
    save_bounds = args.save_bounds

    ome_zarr_args = ome_zarr_args_to_dict(args)
    common_args = common_args_to_dict(args)

    import json
    from amst2.library.input_file_and_dirpaths import make_cache_folder_structure, solve_output_path_and_name

    ome_zarr_filepath, ome_zarr_filename = solve_output_path_and_name(
        ome_zarr_filepath,
        ome_zarr_filename,
        target_dirpath,
        os.path.split(stack_path)[1].replace('.h5', '') + '.ome.zarr'
    )

    cache_dirpath, this_cache_dirpath = make_cache_folder_structure(
        target_dirpath,
        f'stack_to_ome_zarr_{ome_zarr_filename.replace(".ome.zarr", "")}',
        continue_run=args.continue_run
    )

    from squirrel.library.io import load_data_handle
    data_h, shape_h = load_data_handle(stack_path, key=stack_key, pattern=stack_pattern)
    batch_ids = [x for x in range(0, shape_h[0], args.batch_size)]

    src_dirpath = os.path.dirname(os.path.realpath(__file__))

    run_info = dict(
        stack_path=stack_path,
        target_dirpath=target_dirpath,
        ome_zarr_filename=ome_zarr_filename,
        ome_zarr_filepath=ome_zarr_filepath,
        stack_pattern=stack_pattern,
        stack_key=stack_key,
        save_bounds=save_bounds,
        **ome_zarr_args,
        **common_args,
        batch_ids=batch_ids,
        stack_shape=shape_h,
        src_dirpath=src_dirpath
    )

    with open(os.path.join(this_cache_dirpath, 'run.json'), mode='w') as f:
        json.dump(run_info, f, indent=2)

    from snakemake.cli import parse_args, args_to_api
    parser, sn_args = parse_args({})
    args_to_snakemake_arguments(args, sn_args)

    from pathlib import Path
    sn_args.snakefile = Path(os.path.join(src_dirpath, 'snakemake_workflows/stack_to_ome_zarr.snk'))
    sn_args.set_threads = dict(batch_to_ome_zarr=min(args.batch_size, args.cores))
    sn_args.directory = Path(this_cache_dirpath)

    if args.cluster is not None:
        # from amst2.library.data import estimate_mem_mb
        from amst2.cluster.slurm import get_cluster_settings

        sn_args.set_resources = dict(
            batch_to_ome_zarr=dict(
                # I'm purposely not using snakemake's functionality here to determine mem_mb on-the-fly
                mem_mb=1600,  #int(np.ceil(estimate_mem_mb(data_h) * batch_size * 8)),
                runtime=10
            ),
            create_ome_zarr=dict(
                mem_mb=100,
                runtime=5
            )
        )

        sn_args = get_cluster_settings(sn_args, os.path.join(src_dirpath, 'cluster', 'embl.json'))

    args_to_api(sn_args, parser)


if __name__ == '__main__':
    snk_stack_to_ome_zarr()


import os


def snk_stack_to_ome_zarr():

    # Argument parsing -----------------------------------

    import argparse
    from amst2.library.argument_parsing import (
        add_common_arguments_to_parser,
        add_output_locations_to_parser, output_locations_to_dict,
        add_ome_zarr_arguments_to_parser,
        add_snakemake_arguments_to_parser, args_to_snakemake_arguments,
        args_to_dict,
        write_run_json
    )

    parser = argparse.ArgumentParser(
        description='Convert a dataset within a h5 container or a tif stack to ome.zarr',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('stack_path', type=str,
                        help='Input h5 container or tif stack')
    parser.add_argument('--stack_pattern', type=str, default='*.tif',
                        help='File pattern for globbing the input stack; default="*.tif"')
    parser.add_argument('--stack_key', type=str, default='data',
                        help='Path within input h5 file; default="data"')
    parser.add_argument('--save_bounds', action='store_true',
                        help='Saves a json file alongside that contains the bounds of the non-zero area of each slice')

    common_arg_fields = add_common_arguments_to_parser(parser)
    add_output_locations_to_parser(parser)
    ome_zarr_fields = add_ome_zarr_arguments_to_parser(parser)
    add_snakemake_arguments_to_parser(parser)

    args = parser.parse_args()

    stack_path = os.path.abspath(args.stack_path)
    stack_pattern = args.stack_pattern
    stack_key = args.stack_key
    save_bounds = args.save_bounds

    ome_zarr_args = args_to_dict(args, ome_zarr_fields)
    common_args = args_to_dict(args, common_arg_fields)
    output_location_args = output_locations_to_dict(
        args,
        workflow_name='stack_to_ome_zarr',
        fallback_output_name=os.path.split(stack_path)[1].replace('.h5', '') + '.ome.zarr'
    )

    assert common_args['batch_size'] in [4, 8, 16, 32, 64], 'Only allowing batch sizes of [4, 8, 16, 32, 64]!'
    assert common_args['batch_size'] % ome_zarr_args['chunk_size'][0] == 0
    chunk_size = [ome_zarr_args['chunk_size']]
    downsample_factors = ome_zarr_args['downsample_factors']
    for ds_factor in downsample_factors:
        z_chunk = int(chunk_size[-1][0] / ds_factor)
        if z_chunk == 0:
            z_chunk = 1
        chunk_size.append([z_chunk, chunk_size[0][1], chunk_size[0][1]])
    ome_zarr_args['chunk_size'] = chunk_size

    # Generate run.json ----------------------------------

    from squirrel.library.io import load_data_handle
    data_h, shape_h = load_data_handle(stack_path, key=stack_key, pattern=stack_pattern)
    batch_ids = [x for x in range(0, shape_h[0], args.batch_size)]

    src_dirpath = os.path.dirname(os.path.realpath(__file__))
    dtype = str(data_h[0].dtype)

    run_info = dict(
        stack_path=stack_path,
        stack_pattern=stack_pattern,
        stack_key=stack_key,
        save_bounds=save_bounds,
        batch_ids=batch_ids,
        stack_shape=shape_h,
        src_dirpath=src_dirpath,
        dtype=dtype,
        **output_location_args,
        **ome_zarr_args,
        **common_args
    )

    write_run_json(run_info, output_location_args)

    # Snakemake stuff ------------------------------------

    from snakemake.cli import parse_args, args_to_api
    from pathlib import Path

    parser, sn_args = parse_args({})
    args_to_snakemake_arguments(args, sn_args, output_location_args=output_location_args)
    sn_args.snakefile = Path(os.path.join(src_dirpath, 'snakemake_workflows/stack_to_ome_zarr.snk'))
    sn_args.set_threads = dict(batch_to_ome_zarr=min(args.batch_size, args.cores))

    if args.cluster is not None:
        # from amst2.library.data import estimate_mem_mb
        from amst2.cluster.slurm import get_cluster_settings

        sn_args.set_resources = dict(
            batch_to_ome_zarr=dict(
                # I'm purposely not using snakemake's functionality here to determine mem_mb on-the-fly
                mem_mb=8000,  # int(np.ceil(estimate_mem_mb(data_h) * batch_size * 8)),
                runtime=10
            ),
            create_ome_zarr=dict(
                mem_mb=1024,
                runtime=5
            )
        )

        sn_args = get_cluster_settings(sn_args, os.path.join(src_dirpath, 'cluster', 'embl.json'))

    args_to_api(sn_args, parser)


def snk_ome_zarr_to_stack():

    # Argument parsing -----------------------------------

    import argparse
    from amst2.library.argument_parsing import (
        add_common_arguments_to_parser,
        add_snakemake_arguments_to_parser, args_to_snakemake_arguments,
        args_to_dict,
        write_run_json
    )

    parser = argparse.ArgumentParser(
        description='Convert a dataset within a h5 container or a tif stack to ome.zarr',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('ome_zarr_filepath', type=str,
                        help='Ome-zarr dataset which will be converted to tif stack')
    parser.add_argument('target_dirpath', type=str,
                        help='Output directory which will contain the tif slices.')
    parser.add_argument('stack_name', type=str,
                        help='Name of the tif stack (will be a subfolder of the target dirpath)')
    parser.add_argument('--ome_zarr_key', type=str, default='s0',
                        help='Path within input ome-zarr dataset; default="s0"')

    common_arg_fields = add_common_arguments_to_parser(parser)
    add_snakemake_arguments_to_parser(parser)

    args = parser.parse_args()

    ome_zarr_filepath = args.ome_zarr_filepath
    target_dirpath = args.target_dirpath
    stack_name = args.stack_name
    ome_zarr_key = args.ome_zarr_key

    common_args = args_to_dict(args, common_arg_fields)

    from amst2.library.input_file_and_dirpaths import make_cache_folder_structure
    cache_dirpath, this_cache_dirpath = make_cache_folder_structure(
        target_dirpath,
        f'ome_zarr_to_stack_{stack_name}',
        continue_run=common_args['continue_run']
    )

    output_location_args = dict(
        cache_dirpath=cache_dirpath,
        this_cache_dirpath=this_cache_dirpath
    )

    tif_stack_dirpath = os.path.join(target_dirpath, stack_name)

    # Generate run.json ----------------------------------

    from squirrel.library.ome_zarr import get_ome_zarr_handle
    data_h = get_ome_zarr_handle(ome_zarr_filepath, key=ome_zarr_key, mode='r')
    batch_ids = [x for x in range(0, data_h.shape[0], args.batch_size)]

    src_dirpath = os.path.dirname(os.path.realpath(__file__))
    dtype = str(data_h[0].dtype)

    run_info = dict(
        ome_zarr_filepath=ome_zarr_filepath,
        ome_zarr_key=ome_zarr_key,
        target_dirpath=target_dirpath,
        batch_ids=batch_ids,
        stack_shape=data_h.shape,
        src_dirpath=src_dirpath,
        dtype=dtype,
        cache_dirpath=cache_dirpath,
        this_cache_dirpath=this_cache_dirpath,
        tif_stack_dirpath=tif_stack_dirpath,
        **common_args
    )

    write_run_json(run_info, output_location_args)

    # Snakemake stuff ------------------------------------

    from snakemake.cli import parse_args, args_to_api
    from pathlib import Path

    parser, sn_args = parse_args({})
    args_to_snakemake_arguments(args, sn_args, output_location_args=output_location_args)
    sn_args.snakefile = Path(os.path.join(src_dirpath, 'snakemake_workflows/ome_zarr_to_stack.snk'))
    sn_args.set_threads = dict(ome_zarr_batch_to_stack=min(args.batch_size, args.cores, args.max_cores_per_task))

    if args.cluster is not None:
        from amst2.cluster.slurm import get_cluster_settings

        sn_args.set_resources = dict(
            ome_zarr_batch_to_stack=dict(
                # I'm purposely not using snakemake's functionality here to determine mem_mb on-the-fly
                mem_mb=8000,  # int(np.ceil(estimate_mem_mb(data_h) * batch_size * 8)),
                runtime=10
            )
        )

        sn_args = get_cluster_settings(sn_args, os.path.join(src_dirpath, 'cluster', 'embl.json'))

    args_to_api(sn_args, parser)


if __name__ == '__main__':
    snk_stack_to_ome_zarr()

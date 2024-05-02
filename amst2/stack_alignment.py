
import os


def snk_default_amst_pre_alignment():

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
        description='Runs the default amst pre-alignment workflow',
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('input_ome_zarr_filepath', type=str,
                        help='Input ome zarr filepath')
    parser.add_argument('--local_auto_mask', action='store_true',
                        help='Generates a mask using image > 0')
    parser.add_argument('--no_local_alignment', action='store_true',
                        help='Use this flag to turn off the local (slice-to-slice) alignment')
    parser.add_argument('--template_roi', type=float, nargs=5, default=None,
                        metavar=('x-min', 'y-min', 'x-max', 'y-max', 'z'),
                        help='If supplied template matching will be performed with the here defined template\n'
                             'Note: the values for x, y and z are in dataset units!')
    parser.add_argument('--tm_search_roi', type=float, nargs=4, default=None,
                        metavar=('x-min', 'y-min', 'x-max', 'y-max'),
                        help='The template will only be matched in this area.\n'
                             'Note: the values for x, y and z are in dataset units!')
    parser.add_argument('-cm', '--combine_median', type=int, default=8,
                        help='Median smoothing of offsets when combining local and TM')
    parser.add_argument('-cs', '--combine_sigma', type=float, default=8.,
                        help='Gaussian smoothing of offsets when combining local and TM')
    parser.add_argument('--no_previews', action='store_true',
                        help='No preview outputs are generated for the alignement steps')
    parser.add_argument('--preview_downsample_level', type=int, default=2,
                        help='Downsample level of preview volumes; default=2')

    common_arg_fields = add_common_arguments_to_parser(parser)
    add_output_locations_to_parser(parser)
    ome_zarr_fields = add_ome_zarr_arguments_to_parser(parser, omit=['resolution'])
    add_snakemake_arguments_to_parser(parser)

    args = parser.parse_args()

    input_ome_zarr_filepath = args.input_ome_zarr_filepath
    local_auto_mask = args.local_auto_mask
    no_local_alignment = args.no_local_alignment
    template_roi = args.template_roi
    combine_median = args.combine_median
    combine_sigma = args.combine_sigma
    no_previews = args.no_previews
    tm_search_roi = args.tm_search_roi
    preview_downsample_level = args.preview_downsample_level

    common_args = args_to_dict(args, common_arg_fields)
    output_location_args = output_locations_to_dict(
        args,
        workflow_name='stack_alignment',
        fallback_output_name='pre-align.ome.zarr'
    )
    ome_zarr_args = args_to_dict(args, ome_zarr_fields)
    meta_dirpath = os.path.join(
        output_location_args['target_dirpath'],
        f'{output_location_args["output_ome_zarr_filename"].replace(".ome.zarr", "")}.meta'
    )
    if not os.path.exists(meta_dirpath):
        os.mkdir(meta_dirpath)

    # Generate run.json ----------------------------------

    from squirrel.library.io import load_data_handle
    from squirrel.library.ome_zarr import (
        get_scale_of_downsample_level, get_ome_zarr_handle, get_unit_of_dataset
    )
    data_h, shape_h = load_data_handle(input_ome_zarr_filepath, key='s0', pattern=None)
    batch_ids = [x for x in range(0, shape_h[0], common_args['batch_size'])]
    ome_zarr_h = get_ome_zarr_handle(input_ome_zarr_filepath, key=None, mode='r')
    resolution = get_scale_of_downsample_level(ome_zarr_h, 0)
    unit = get_unit_of_dataset(ome_zarr_h)
    dtype = str(data_h.dtype)

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

    src_dirpath = os.path.dirname(os.path.realpath(__file__))

    run_info = dict(
        input_ome_zarr_filepath=input_ome_zarr_filepath,
        no_local_alignment=no_local_alignment,
        use_template_matching=template_roi is not None,
        combine_median=combine_median,
        combine_sigma=combine_sigma,
        no_previews=no_previews,
        batch_ids=batch_ids,
        stack_shape=shape_h,
        src_dirpath=src_dirpath,
        preview_downsample_level=preview_downsample_level,
        output_dtype=dtype,
        resolution=resolution,
        unit=unit,
        meta_dirpath=meta_dirpath,
        local_alignment_params=dict(
            auto_mask=local_auto_mask,
            transform='translation',
            key='s0',
            pre_fix_big_jumps=True  # Hard-coding because if slices are off it really doesn't work otherwise!
        ),
        template_matching_params=dict(
            template_roi=template_roi,
            search_roi=tm_search_roi,
            resolution=resolution,
            key='s0',
            save_template=True
        ),
        **output_location_args,
        **common_args,
        **ome_zarr_args
    )

    write_run_json(run_info, output_location_args)

    # Snakemake stuff ------------------------------------

    from snakemake.cli import parse_args, args_to_api
    from pathlib import Path

    parser, sn_args = parse_args({})
    args_to_snakemake_arguments(args, sn_args, output_location_args=output_location_args)
    sn_args.snakefile = Path(os.path.join(src_dirpath, 'snakemake_workflows/default_amst_pre_alignment.snk'))
    sn_args.set_threads = dict(
        elastix_stack_alignment=min(args.batch_size, args.cores, args.max_cores_per_task),
        template_matching=min(args.batch_size, args.cores, args.max_cores_per_task),
        apply_final_transformations=min(args.batch_size, args.cores, args.max_cores_per_task)
    )

    if args.cluster is not None:

        from amst2.cluster.slurm import get_cluster_settings

        sn_args.set_resources = dict(
            elastix_stack_alignment=dict(
                mem_mb=16000,
                runtime=30,
            ),
            local_alignment_preview=dict(
                mem_mb=8000,
                runtime=30
            ),
            template_matching=dict(
                mem_mb=24000,
                runtime=30
            ),
            template_matching_preview=dict(
                mem_mb=8000,
                runtime=30
            ),
            finalize_and_join_transforms=dict(
                mem_mb=100,
                runtime=5
            ),
            final_preview=dict(
                mem_mb=8000,
                runtime=30
            ),
            create_ome_zarr=dict(
                mem_mb=1024,
                runtime=5
            ),
            apply_final_transformations=dict(
                mem_mb=8000,
                runtime=30
            )
        )

        sn_args = get_cluster_settings(sn_args, os.path.join(src_dirpath, 'cluster', 'embl.json'))

    args_to_api(sn_args, parser)


if __name__ == '__main__':
    snk_default_amst_pre_alignment()

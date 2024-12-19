
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
    parser.add_argument('--local_gaussian_sigma', type=float, default=0.,
                        help='Perform a gaussian smoothing before registration')
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
    parser.add_argument('--mem', type=str, nargs='+', default=None,
                        help='Cluster job memory amounts for the snakemake rules.\n'
                             'For each rule define like so:\n'
                             '"rule_name:16000"')
    parser.add_argument('--runtime', type=str, nargs='+', default=None,
                        help='Cluster job runtimes for the snakemake rules. \n'
                             'For each rule define like so:\n'
                             '"rule_name:30"')
    parser.add_argument('--preview_downsample_level', type=int, default=2,
                        help='Downsample level of preview volumes; default=2')

    common_arg_fields = add_common_arguments_to_parser(parser)
    add_output_locations_to_parser(parser)
    ome_zarr_fields = add_ome_zarr_arguments_to_parser(parser, omit=['resolution'])
    add_snakemake_arguments_to_parser(parser)

    args = parser.parse_args()

    input_ome_zarr_filepath = os.path.abspath(args.input_ome_zarr_filepath)
    local_auto_mask = args.local_auto_mask
    local_gaussian_sigma = args.local_gaussian_sigma
    no_local_alignment = args.no_local_alignment
    template_roi = args.template_roi
    combine_median = args.combine_median
    combine_sigma = args.combine_sigma
    no_previews = args.no_previews
    tm_search_roi = args.tm_search_roi
    preview_downsample_level = args.preview_downsample_level
    mem = args.mem
    runtime = args.runtime

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
        elastix_stack_alignment_workflow_params=dict(
            auto_mask=local_auto_mask,
            gaussian_sigma=local_gaussian_sigma,
            transform='translation',
            key='s0',
            pre_fix_big_jumps=True,  # Hard-coding because if slices are off it really doesn't work otherwise!
            pre_fix_iou_thresh=0.9
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
        from amst2.cluster.general import set_resources

        set_resources(
            sn_args,
            [
                'elastix_stack_alignment',
                'local_alignment_preview',
                'template_matching',
                'template_matching_preview',
                'finalize_and_join_transforms',
                'final_preview',
                'create_ome_zarr',
                'apply_final_transformations'
            ],
            [16000, 8000, 24000, 8000, 100, 8000, 1024, 8000],
            [30, 30, 30, 30, 5, 30, 5, 30],
            mem_args=mem,
            runtime_args=runtime
        )

        sn_args = get_cluster_settings(sn_args, os.path.join(src_dirpath, 'cluster', 'embl.json'))

    args_to_api(sn_args, parser)


def snk_elastix_stack_alignment():

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
        description='Runs an elastix stack alignment workflow',
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('input_ome_zarr_filepath', type=str,
                        help='Input ome zarr filepath')
    parser.add_argument('--transform', type=str, default='translation',
                        help='The transformation that is applied to the images; default="translation"')
    parser.add_argument('--auto_mask', action='store_true',
                        help='Generates a mask using image > 0')
    parser.add_argument('--gaussian_sigma', type=float, default=0.,
                        help='Perform a gaussian smoothing before registration')
    parser.add_argument('--use_edges', action='store_true',
                        help='Computes edges before registration using a sobel filter')
    parser.add_argument('--z_step', type=int, default=1,
                        help='Uses only every nth data slice; good to determine long range correspondence; default=1')
    parser.add_argument('--average_for_z_step', action='store_true',
                        help='Use the average of n=z_step slices for fixed and moving images')
    parser.add_argument('--auto_pad', action='store_true',
                        help='Compute the information necessary for auto-padding and apply if requested')
    parser.add_argument('--determine_bounds', action='store_true',
                        help='Determine the bounds of each slice; Does not work with z_step > 1; Implicit for auto_pad')
    parser.add_argument('--apply_final', action='store_true',
                        help='Apply the final transformation and create result ome-zarr')
    parser.add_argument('--no_preview', action='store_true',
                        help='No preview outputs are generated for the alignement steps')
    parser.add_argument('--no_fixing_of_big_jumps', action='store_true',
                        help='Switches off usage of cross-correlation to fix big initial jumps\n'
                             'Only use this if you are sure that there are no big jumps in the input stack!')
    parser.add_argument('--pre_fix_iou_thresh', type=float, default=0.9,
                        help='IoU-threshold below which the big jumps are pre-aligned using cross-correlation')
    parser.add_argument('--parameter_map', type=str, default=None,
                        help='Elastix parameter map file')
    parser.add_argument('--mem', type=str, nargs='+', default=None,
                        help='Cluster job memory amounts for the snakemake rules.\n'
                             'For each rule define like so:\n'
                             '"rule_name:16000"')
    parser.add_argument('--runtime', type=str, nargs='+', default=None,
                        help='Cluster job runtimes for the snakemake rules. \n'
                             'For each rule define like so:\n'
                             '"rule_name:30"')
    parser.add_argument('--preview_downsample_level', type=int, default=2,
                        help='Downsample level of preview volumes; default=2')

    common_arg_fields = add_common_arguments_to_parser(parser)
    add_output_locations_to_parser(parser)
    ome_zarr_fields = add_ome_zarr_arguments_to_parser(parser, omit=['resolution'])
    add_snakemake_arguments_to_parser(parser)

    args = parser.parse_args()

    input_ome_zarr_filepath = os.path.abspath(args.input_ome_zarr_filepath)
    transform = args.transform
    auto_mask = args.auto_mask
    gaussian_sigma = args.gaussian_sigma
    use_edges = args.use_edges
    z_step = args.z_step
    average_for_z_step = args.average_for_z_step
    auto_pad = args.auto_pad
    determine_bounds = args.determine_bounds
    apply_final = args.apply_final
    no_preview = args.no_preview
    no_fixing_of_big_jumps = args.no_fixing_of_big_jumps
    pre_fix_iou_thresh = args.pre_fix_iou_thresh
    parameter_map = os.path.abspath(args.parameter_map) if args.parameter_map is not None else None
    preview_downsample_level = args.preview_downsample_level
    mem = args.mem
    runtime = args.runtime

    common_args = args_to_dict(args, common_arg_fields)
    output_location_args = output_locations_to_dict(
        args,
        workflow_name='stack_alignment',
        fallback_output_name='elastix.ome.zarr'
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
    print(f'input_ome_zarr_filepath = {input_ome_zarr_filepath}')
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
        no_preview=no_preview,
        batch_ids=batch_ids,
        stack_shape=shape_h,
        src_dirpath=src_dirpath,
        preview_downsample_level=preview_downsample_level,
        output_dtype=dtype,
        resolution=resolution,
        unit=unit,
        meta_dirpath=meta_dirpath,
        auto_pad=auto_pad,
        determine_bounds=determine_bounds,
        apply_final=apply_final,
        elastix_stack_alignment_workflow_params=dict(
            auto_mask=auto_mask,
            gaussian_sigma=gaussian_sigma,
            use_edges=use_edges,
            z_step=z_step,
            average_for_z_step=average_for_z_step,
            transform=transform,
            parameter_map=parameter_map,
            key='s0',
            pre_fix_big_jumps=not no_fixing_of_big_jumps,
            pre_fix_iou_thresh=pre_fix_iou_thresh
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
    sn_args.snakefile = Path(os.path.join(src_dirpath, 'snakemake_workflows/elastix_stack_alignment.snk'))
    sn_args.set_threads = dict(
        elastix_stack_alignment=min(args.batch_size, args.cores, args.max_cores_per_task),
        apply_final_transformations=min(args.batch_size, args.cores, args.max_cores_per_task)
    )

    if args.cluster is not None:

        from amst2.cluster.slurm import get_cluster_settings
        from amst2.cluster.general import set_resources

        set_resources(
            sn_args,
            [
                'elastix_stack_alignment',
                'alignment_preview',
                'create_ome_zarr',
                'apply_final_transformations'
            ],
            [16000, 8000, 1024, 8000],
            [30, 30, 5, 30],
            mem_args=mem,
            runtime_args=runtime
        )

        sn_args = get_cluster_settings(sn_args, os.path.join(src_dirpath, 'cluster', 'embl.json'))

    args_to_api(sn_args, parser)


def snk_apply_transformation():

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
        description='Applies a stack of transformations to an image stack',
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('input_ome_zarr_filepath', type=str,
                        help='Input ome zarr filepath')
    parser.add_argument('input_transforms_filepaths', type=str, nargs='+',
                        help='Filepaths to transformations file (*.json) or directory')
    parser.add_argument('--no_autopad', action='store_true',
                        help='Switches off auto-padding')
    parser.add_argument('--mem', type=str, nargs='+', default=None,
                        help='Cluster job memory amounts for the snakemake rules.\n'
                             'For each rule define like so:\n'
                             '"rule_name:16000"')
    parser.add_argument('--runtime', type=str, nargs='+', default=None,
                        help='Cluster job runtimes for the snakemake rules. \n'
                             'For each rule define like so:\n'
                             '"rule_name:30"')

    common_arg_fields = add_common_arguments_to_parser(parser)
    add_output_locations_to_parser(parser)
    ome_zarr_fields = add_ome_zarr_arguments_to_parser(parser, omit=['resolution'])
    add_snakemake_arguments_to_parser(parser)

    args = parser.parse_args()

    input_ome_zarr_filepath = os.path.abspath(args.input_ome_zarr_filepath)
    input_transforms_filepaths = [os.path.abspath(x) for x in args.input_transforms_filepaths]
    no_autopad = args.no_autopad
    mem = args.mem
    runtime = args.runtime

    common_args = args_to_dict(args, common_arg_fields)
    output_location_args = output_locations_to_dict(
        args,
        workflow_name='apply_transforms',
        fallback_output_name=f'{os.path.splitext(os.path.split(input_transforms_filepaths[0])[1])[0]}.ome.zarr'
    )
    ome_zarr_args = args_to_dict(args, ome_zarr_fields)

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

    shapes = []
    for input_transforms_filepath in input_transforms_filepaths:
        if os.path.isdir(input_transforms_filepath):
            from squirrel.library.elastix import ElastixStack
            this_stack = ElastixStack(dirpath=input_transforms_filepath)
            if this_stack.image_shape() is not None:
                shapes.append(this_stack.image_shape())
        else:
            from squirrel.library.affine_matrices import AffineStack
            this_stack = AffineStack(filepath=input_transforms_filepath)
            if this_stack.exists_meta('stack_shape') and this_stack.get_meta('stack_shape') is not None:
                shapes.append(this_stack.get_meta('stack_shape')[1:])
    if len(shapes) > 0:
        import numpy as np
        shapes = np.array(shapes)
        assert np.all(shapes == shapes[0], axis=0).all(), 'Only allowed for transform stacks with equal image shapes'
        stack_shape = np.array([shape_h[0], shapes[0][0], shapes[0][1]]).tolist()
    else:
        stack_shape = shape_h

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
        input_transforms_filepaths=input_transforms_filepaths,
        batch_ids=batch_ids,
        stack_shape=stack_shape,
        src_dirpath=src_dirpath,
        output_dtype=dtype,
        resolution=resolution,
        unit=unit,
        no_autopad=no_autopad,
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
    sn_args.snakefile = Path(os.path.join(src_dirpath, 'snakemake_workflows/apply_transformation.snk'))
    sn_args.set_threads = dict(
        apply_final_transformations=min(args.batch_size, args.cores, args.max_cores_per_task)
    )

    if args.cluster is not None:

        from amst2.cluster.slurm import get_cluster_settings
        from amst2.cluster.general import set_resources

        set_resources(
            sn_args,
            [
                'create_ome_zarr',
                'apply_final_transformations'
            ],
            [1024, 8000],
            [5, 30],
            mem_args=mem,
            runtime_args=runtime
        )

        sn_args = get_cluster_settings(sn_args, os.path.join(src_dirpath, 'cluster', 'embl.json'))

    args_to_api(sn_args, parser)


def snk_template_matching_stack_alignment():
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
        description='Runs a template matching stack alignment workflow',
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('input_ome_zarr_filepath', type=str,
                        help='Input ome zarr filepath')
    parser.add_argument('--template_roi', type=float, nargs=5, default=None,
                        metavar=('x-min', 'y-min', 'x-max', 'y-max', 'z'),
                        help='This defines the template\n'
                             'Note: the values for x, y and z are in dataset units!')
    parser.add_argument('--tm_search_roi', type=float, nargs=4, default=None,
                        metavar=('x-min', 'y-min', 'x-max', 'y-max'),
                        help='The template will only be matched in this area.\n'
                             'Note: the values for x, y and z are in dataset units!')
    parser.add_argument('-cm', '--combine_median', type=int, default=8,
                        help='Median smoothing of offsets when combining local and TM')
    parser.add_argument('-cs', '--combine_sigma', type=float, default=8.,
                        help='Gaussian smoothing of offsets when combining local and TM')
    parser.add_argument('--auto_pad', action='store_true',
                        help='Compute the information necessary for auto-padding and apply if requested')
    parser.add_argument('--determine_bounds', action='store_true',
                        help='Determine the bounds of each slice; Implicit for auto_pad')
    parser.add_argument('--apply_final', action='store_true',
                        help='Apply the final transformation and create result ome-zarr')
    parser.add_argument('--no_preview', action='store_true',
                        help='No preview outputs are generated for the alignement steps')
    parser.add_argument('--mem', type=str, nargs='+', default=None,
                        help='Cluster job memory amounts for the snakemake rules.\n'
                             'For each rule define like so:\n'
                             '"rule_name:16000"')
    parser.add_argument('--runtime', type=str, nargs='+', default=None,
                        help='Cluster job runtimes for the snakemake rules. \n'
                             'For each rule define like so:\n'
                             '"rule_name:30"')
    parser.add_argument('--preview_downsample_level', type=int, default=2,
                        help='Downsample level of preview volumes; default=2')

    common_arg_fields = add_common_arguments_to_parser(parser)
    add_output_locations_to_parser(parser)
    ome_zarr_fields = add_ome_zarr_arguments_to_parser(parser, omit=['resolution'])
    add_snakemake_arguments_to_parser(parser)

    args = parser.parse_args()

    input_ome_zarr_filepath = os.path.abspath(args.input_ome_zarr_filepath)
    template_roi = args.template_roi
    tm_search_roi = args.tm_search_roi
    combine_median = args.combine_median
    combine_sigma = args.combine_sigma
    auto_pad = args.auto_pad
    determine_bounds = args.determine_bounds
    apply_final = args.apply_final
    no_preview = args.no_preview
    preview_downsample_level = args.preview_downsample_level
    mem = args.mem
    runtime = args.runtime

    common_args = args_to_dict(args, common_arg_fields)
    output_location_args = output_locations_to_dict(
        args,
        workflow_name='stack_alignment',
        fallback_output_name='elastix.ome.zarr'
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
    print(f'input_ome_zarr_filepath = {input_ome_zarr_filepath}')
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
        no_preview=no_preview,
        batch_ids=batch_ids,
        stack_shape=shape_h,
        src_dirpath=src_dirpath,
        preview_downsample_level=preview_downsample_level,
        output_dtype=dtype,
        resolution=resolution,
        unit=unit,
        meta_dirpath=meta_dirpath,
        auto_pad=auto_pad,
        determine_bounds=determine_bounds,
        apply_final=apply_final,
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
    sn_args.snakefile = Path(os.path.join(src_dirpath, 'snakemake_workflows/template_matching_alignment.snk'))
    sn_args.set_threads = dict(
        template_matching=min(args.batch_size, args.cores, args.max_cores_per_task),
        apply_final_transformations=min(args.batch_size, args.cores, args.max_cores_per_task)
    )

    if args.cluster is not None:

        from amst2.cluster.slurm import get_cluster_settings
        from amst2.cluster.general import set_resources

        set_resources(
            sn_args,
            [
                'template_matching',
                'template_matching_preview',
                'create_ome_zarr',
                'apply_final_transformations'
            ],
            [16000, 8000, 1024, 8000],
            [30, 30, 5, 30],
            mem_args=mem,
            runtime_args=runtime
        )

        sn_args = get_cluster_settings(sn_args, os.path.join(src_dirpath, 'cluster', 'embl.json'))

    args_to_api(sn_args, parser)


if __name__ == '__main__':
    snk_default_amst_pre_alignment()

import os


def _run_nsbs_pre_align(parameter_yaml, verbose=False):

    from .lib import load_parameter_yaml
    parameter_dict = load_parameter_yaml(parameter_yaml)

    output_dirpath = parameter_dict['general']['output_dirpath']

    from .lib import run_stack_to_ome_zarr
    run_stack_to_ome_zarr(parameter_dict, parameter_key='stack_to_ome_zarr', verbose=verbose)

    if not os.path.exists(os.path.join(output_dirpath, 'stack_to_ome_zarr', 'stack_to_ome_zarr.done')):
        return

    from .lib import run_nsbs_alignment
    parameter_dict['sbs_alignment']['apply_final'] = True
    run_nsbs_alignment(
        parameter_dict, parameter_key='sbs_alignment',
        input_dirpath=os.path.join(output_dirpath, 'stack_to_ome_zarr/input-raw.ome.zarr'),
        verbose=verbose
    )

    if not os.path.exists(os.path.join(output_dirpath, 'sbs_alignment', 'sbs_alignment.done')):
        return

    from .lib import run_nsbs_alignment
    run_nsbs_alignment(
        parameter_dict, parameter_key='nsbs_alignment',
        input_dirpath=os.path.join(output_dirpath, 'sbs_alignment/sbs.ome.zarr'),
        verbose=verbose
    )

    if not os.path.exists(os.path.join(output_dirpath, 'nsbs_alignment', 'nsbs_alignment.done')):
        return

    from squirrel.workflows.transformation import dot_product_on_affines_workflow
    dot_product_on_affines_workflow(
        [
            os.path.join(output_dirpath, 'sbs_alignment/sbs.meta/elastix.json'),
            os.path.join(output_dirpath, 'nsbs_alignment/nsbs.meta/elastix.json')
        ],
        os.path.join(output_dirpath, 'nsbs_alignment/nsbs.meta/combined.json'),
        keep_meta=0,
        verbose=verbose
    )
    if not os.path.exists(os.path.join(output_dirpath, 'nsbs_alignment/nsbs.meta/combined.json')):
        return
    from squirrel.workflows.transformation import apply_auto_pad_workflow
    apply_auto_pad_workflow(
        os.path.join(output_dirpath, 'nsbs_alignment/nsbs.meta/combined.json'),
        os.path.join(output_dirpath, 'nsbs-pre-align.json'),
        verbose=verbose
    )
    if not os.path.exists(os.path.join(output_dirpath, 'nsbs-pre-align.json')):
        return

    from .lib import run_apply_transformation
    run_apply_transformation(
        parameter_dict,
        parameter_key='apply_transformation',
        transforms_filepath=os.path.join(output_dirpath, 'nsbs-pre-align.json'),
        input_dirpath=os.path.join(output_dirpath, 'stack_to_ome_zarr/input-raw.ome.zarr'),
        output_filename='nsbs-pre-align.ome.zarr',
        verbose=verbose
    )

    if not os.path.exists(os.path.join(output_dirpath, 'apply_transformation', 'apply_transformation.done')):
        return

    if 'result_to_tif_stack' in parameter_dict and parameter_dict['result_to_tif_stack']['active']:
        from .lib import run_ome_zarr_to_stack
        run_ome_zarr_to_stack(
            parameter_dict,
            parameter_key='result_to_tif_stack',
            input_dirpath=os.path.join(output_dirpath, 'apply_transformation/nsbs-pre-align.ome.zarr'),
            output_dirname='nsbs-pre-align',
            verbose=verbose
        )


def nsbs_pre_align():

    # Argument parsing -----------------------------------

    import argparse

    parser = argparse.ArgumentParser(
        description='Runs a pre-alignment with nth-slice-by-slice alignment for long-range registration',
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('parameter_yaml', type=str)
    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()

    parameter_yaml = args.parameter_yaml
    verbose = args.verbose

    if verbose:
        print(f'parameter_yaml = {parameter_yaml}')

    _run_nsbs_pre_align(parameter_yaml, verbose=verbose)

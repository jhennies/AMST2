
import os


def _run_amst(parameter_yaml, verbose=False):

    from .lib import load_parameter_yaml
    parameter_dict = load_parameter_yaml(parameter_yaml)

    output_dirpath = parameter_dict['general']['output_dirpath']

    from .lib import run_amst
    run_amst(
        parameter_dict,
        parameter_key='amst',
        pre_align_dirpath=os.path.join(output_dirpath, 'apply_pre_align', 'nsbs-pre-align.ome.zarr'),
        verbose=verbose
    )

    if not os.path.exists(os.path.join(output_dirpath, 'amst', 'amst.done')):
        return

    from .lib import run_apply_transformation
    apply_amst_param_key = 'apply_amst'
    pre_align_transforms = os.path.join(output_dirpath, 'nsbs-pre-align.json')
    pre_align_transforms = parameter_dict['general']['pre_align_transforms'] if 'pre_align_transforms' in parameter_dict['general'] else pre_align_transforms
    pre_align_transforms = parameter_dict[apply_amst_param_key]['pre_align_transforms'] if 'pre_align_transforms' in parameter_dict[apply_amst_param_key] else pre_align_transforms
    run_apply_transformation(
        parameter_dict,
        parameter_key=apply_amst_param_key,
        transforms_filepath=[
            pre_align_transforms,
            os.path.join(output_dirpath, 'amst', 'amst.meta', 'amst')
        ],
        input_dirpath=os.path.join(output_dirpath, 'stack_to_ome_zarr', 'input-raw.ome.zarr'),
        output_filename='amst.ome.zarr',
        verbose=verbose
    )

    if not os.path.exists(os.path.join(output_dirpath, apply_amst_param_key, f'{apply_amst_param_key}.done')):
        return

    if 'amst_to_tif_stack' in parameter_dict and parameter_dict['amst_to_tif_stack']['active']:
        from .lib import run_ome_zarr_to_stack
        run_ome_zarr_to_stack(
            parameter_dict,
            parameter_key='amst_to_tif_stack',
            input_dirpath=os.path.join(output_dirpath, apply_amst_param_key, 'amst.ome.zarr'),
            output_dirname='amst',
            verbose=verbose
        )


def amst():

    # Argument parsing -----------------------------------

    import argparse

    parser = argparse.ArgumentParser(
        description='Runs AMST',
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('parameter_yaml', type=str)
    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()

    parameter_yaml = args.parameter_yaml
    verbose = args.verbose

    if verbose:
        print(f'parameter_yaml = {parameter_yaml}')

    _run_amst(parameter_yaml, verbose=verbose)


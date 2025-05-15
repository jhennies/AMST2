import os


def _run_nsbs_pre_align(parameter_yaml, verbose=False):

    from .lib import load_parameter_yaml
    parameter_dict = load_parameter_yaml(parameter_yaml)
    if verbose:
        print(f'parameter_dict = {parameter_dict}')

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
        parameter_key='apply_pre_align',
        transforms_filepath=os.path.join(output_dirpath, 'nsbs-pre-align.json'),
        input_dirpath=os.path.join(output_dirpath, 'stack_to_ome_zarr/input-raw.ome.zarr'),
        output_filename='nsbs-pre-align.ome.zarr',
        verbose=verbose
    )

    if not os.path.exists(os.path.join(output_dirpath, 'apply_pre_align', 'apply_pre_align.done')):
        return

    if 'pre_align_to_tif_stack' in parameter_dict and parameter_dict['pre_align_to_tif_stack']['active']:
        from .lib import run_ome_zarr_to_stack
        run_ome_zarr_to_stack(
            parameter_dict,
            parameter_key='pre_align_to_tif_stack',
            input_dirpath=os.path.join(output_dirpath, 'apply_pre_align/nsbs-pre-align.ome.zarr'),
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


def _run_cleanup_nsbs_pre_align(parameter_yaml, verbose=False, dryrun=False):

    import shutil
    from glob import glob

    from .lib import load_parameter_yaml
    parameter_dict = load_parameter_yaml(parameter_yaml)

    output_dirpath = parameter_dict['general']['output_dirpath']

    # rm -r */snk_cache
    def remove_snk_caches(project_dir):
        dirpaths = glob(os.path.join(project_dir, '*', 'snk_cache'))
        for dirpath in dirpaths:
            print(f'rm -r {dirpath}')
            if not dryrun:
                shutil.rmtree(dirpath)
    remove_snk_caches(output_dirpath)

    # rm */*.done
    def remove_done_files(project_dir):
        files = glob(os.path.join(project_dir, '*', '*.done'))
        for file in files:
            print(f'rm {file}')
            if not dryrun:
                os.remove(file)
    remove_done_files(output_dirpath)

    # # mv apply_pre_align/nsbs-pre-align.ome.zarr/ .

    # print('rmdir apply_pre_align')
    # if not dryrun:
    #     os.rmdir(os.path.join(output_dirpath, 'apply_pre_align'))

    print('rm -r nsbs_alignment')
    if not dryrun:
        shutil.rmtree(os.path.join(output_dirpath, 'nsbs_alignment'))

    print('rm -r sbs_alignment')
    if not dryrun:
        shutil.rmtree(os.path.join(output_dirpath, 'sbs_alignment'))

    print('mv pre_align_to_tif_stack/nsbs-pre-align/ .')
    if not dryrun:
        shutil.move(os.path.join(output_dirpath, 'pre_align_to_tif_stack', 'nsbs-pre-align'), output_dirpath)

    print('rmdir pre_align_to_tif_stack')
    if not dryrun:
        os.rmdir(os.path.join(output_dirpath, 'pre_align_to_tif_stack'))


def cleanup_nsbs_pre_align():

    # Argument parsing -----------------------------------

    import argparse

    parser = argparse.ArgumentParser(
        description='Cleans up a pre-alignment with nth-slice-by-slice alignment for long-range registration',
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('parameter_yaml', type=str)
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-d', '--dryrun', action='store_true')

    args = parser.parse_args()

    parameter_yaml = args.parameter_yaml
    verbose = args.verbose
    dryrun = args.dryrun

    if verbose:
        print(f'parameter_yaml = {parameter_yaml}')

    _run_cleanup_nsbs_pre_align(parameter_yaml, verbose=verbose, dryrun=dryrun)


import os


def _run_nsbs_pre_align(parameter_yaml, verbose=False):

    from .lib import load_parameter_yaml
    parameter_dict = load_parameter_yaml(parameter_yaml)
    if verbose:
        print(f'parameter_dict = {parameter_dict}')

    output_dirpath = parameter_dict['general']['output_dirpath']

    if not os.path.exists(output_dirpath):
        os.makedirs(output_dirpath, exist_ok=True)

    if 'stack_to_ome_zarr' in parameter_dict and parameter_dict['stack_to_ome_zarr']['active']:
        from .lib import run_stack_to_ome_zarr
        run_stack_to_ome_zarr(parameter_dict, parameter_key='stack_to_ome_zarr', verbose=verbose)
        input_dirpath = os.path.join(output_dirpath, 'stack_to_ome_zarr/input-raw.ome.zarr')
        if not os.path.exists(os.path.join(output_dirpath, 'stack_to_ome_zarr', 'stack_to_ome_zarr.done')):
            return
        parameter_dict['general']['stack_key'] = 's0'
    else:
        input_dirpath = None

    from .lib import run_nsbs_alignment
    parameter_dict['sbs_alignment']['apply_final'] = True
    run_nsbs_alignment(
        parameter_dict, parameter_key='sbs_alignment',
        input_dirpath=input_dirpath,
        verbose=verbose
    )

    if not os.path.exists(os.path.join(output_dirpath, 'sbs_alignment', 'sbs_alignment.done')):
        return

    parameter_dict['nsbs_alignment']['stack_key'] = 's0'
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
    if parameter_dict['general']['auto_pad']:
        from squirrel.workflows.transformation import apply_auto_pad_workflow
        apply_auto_pad_workflow(
            os.path.join(output_dirpath, 'nsbs_alignment/nsbs.meta/combined.json'),
            os.path.join(output_dirpath, 'nsbs-pre-align.json'),
            verbose=verbose
        )
    else:
        import shutil
        shutil.copy(
            os.path.join(output_dirpath, 'nsbs_alignment/nsbs.meta/combined.json'),
            os.path.join(output_dirpath, 'nsbs-pre-align.json')
        )
    if not os.path.exists(os.path.join(output_dirpath, 'nsbs-pre-align.json')):
        return

    from .lib import run_apply_transformation
    run_apply_transformation(
        parameter_dict,
        parameter_key='apply_pre_align',
        transforms_filepath=os.path.join(output_dirpath, 'nsbs-pre-align.json'),
        input_dirpath=input_dirpath,
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


def get_default_parameter_file():

    # Argument parsing -----------------------------------

    import argparse

    parser = argparse.ArgumentParser(
        description='Creates a default parameter file',
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('-fp', '--output_filepath', type=str, default=None,
                        help='Filepath of the parameter file. By default it is created in the current directory')
    parser.add_argument('-pi', '--param_input_dirpath', type=str, default=None,
                        help='The parameter "general:input_dirpath" will be pre-set to this value')
    parser.add_argument('-psk', '--param_stack_key', type=str, default=None,
                        help='Stack key of the input dataset (relevant for h5, n5 or ome.zarr); Defaults to "s0"')
    parser.add_argument('-psp', '--param_stack_pattern', type=str, default=None,
                        help='Stack pattern of the input dataset (relevant for tif stack); Defaults to "*.tif"')
    parser.add_argument('-po', '--param_output_dirpath', type=str, default=None,
                        help='The parameter "general:output_dirpath" will be pre-set to this value')
    parser.add_argument('-p', '--params', nargs='+', default=None,
                        help='Specify any parameter and its value in this format:\n'
                             '  --params group:parameter_name:value\n'
                             '  Examples:\n'
                             '    --params stack_to_ome_zarr:active:true general:cores:32  -> enables input ome zarr conversion and sets the number of compute cores to 32\n'
                             '    --params general:resolution:[0.02,0.02,0.02]  -> Sets resolution to 0.02 isotropic\n'
                             '  Note: One Item consisting of group:parameter_name:value must not contain any spaces!\n'
                             '  Also note: You can generate any entry like this - also ones that have no effect!\n'
                             '  The idea of having this input argument is to make a workflow fully scriptable without the requirement for manually adjusting the parameter file.')
    parser.add_argument('--slurm', action='store_true',
                        help='Creates the parameter file for a slurm cluster; Note that the compute settings may '
                             'require adjustment')
    parser.add_argument('--estimate_crop_xy', action='store_true',
                        help='Estimates how to optimally crop the data to reduce slice size. \n'
                             'Note that enabling this requires conversion of the input to ome.zarr')
    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()

    output_filepath = args.output_filepath
    param_input_dirpath = args.param_input_dirpath
    param_stack_key = args.param_stack_key
    param_stack_pattern = args.param_stack_pattern
    param_output_dirpath = args.param_output_dirpath
    params = args.params
    estimate_crop_xy = args.estimate_crop_xy
    slurm = args.slurm

    verbose = args.verbose

    from amst2.workflows.lib import get_default_parameter_file_from_repo
    get_default_parameter_file_from_repo(
        'pre_align',
        output_filepath=output_filepath,
        param_input_dirpath=param_input_dirpath,
        param_stack_key=param_stack_key,
        param_stack_pattern=param_stack_pattern,
        param_output_dirpath=param_output_dirpath,
        params=params,
        slurm=slurm,
        estimate_crop_xy=estimate_crop_xy,
        verbose=verbose
    )

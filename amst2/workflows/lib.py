
import os


def run_snakemake_workflow(
        run_script, wf_name,
):
    import subprocess
    import sys

    try:
        process = subprocess.run(
            run_script,
            shell=True,
            executable='/bin/bash',
            check=True,
            text=True,
            capture_output=True  # Capture both stdout and stderr
        )
    except subprocess.CalledProcessError as e:
        print(f'\n\n{e.stderr}\n')
        print(f'\n{e.stdout}\n')
        print('_____________________________________________\n')
        print(f'{wf_name} failed. Check out error message above!')
        print('_____________________________________________\n')

        return 1

    if process.stderr.strip().split("\n")[-2].endswith("(100%) done"):
        print('_____________________________________________\n')
        print(f'{wf_name} was successful!')
        print('_____________________________________________\n')

        return 0

    else:
        print(process.stderr)
        print('_____________________________________________\n')
        print(f'{wf_name} failed. ')
        print('Check out {} for details or check out the output above!'.format(process.stderr.strip().split("\n")[-1]))
        print('_____________________________________________\n')

        return 1


def get_parameters_for_snakemake_workflow(parameter_dict, parameter_key, verbose=False):

    this_dict = parameter_dict['general'].copy()
    if verbose:
        print(f'parameter_key = {parameter_key}')
        print(f'this_dict = {this_dict}')

    for k, v in parameter_dict[parameter_key].items():
        if verbose:
            print(f'k = {k}; v = {v}')
        this_dict[k] = v

    print(f'this_dict = {this_dict}')
    return this_dict


def load_parameter_yaml(parameter_yaml):

    import yaml

    with open(parameter_yaml, 'r') as file:
        parameter_dict = yaml.safe_load(file)

    return parameter_dict


def get_resource_str(parameter_dict, resource_name='mem'):

    if resource_name in parameter_dict:
        mem_str = ''
        for mem_k, mem_v in parameter_dict[resource_name].items():
            mem_str += f'{mem_k}:{mem_v} '
        return mem_str

    return ''


def run_stack_to_ome_zarr(
        parameter_dict,
        parameter_key = 'stack_to_ome_zarr',
        verbose=False
):

    this_param_dict = get_parameters_for_snakemake_workflow(parameter_dict, parameter_key)
    output_dirpath = os.path.join(this_param_dict['output_dirpath'], parameter_key)
    if not os.path.exists(output_dirpath):
        os.mkdir(output_dirpath)

    done_fp = os.path.join(output_dirpath, f'{parameter_key}.done')
    if os.path.exists(done_fp):
        print('_____________________________________________\n')
        print(f'{parameter_key} already computed!')
        print('_____________________________________________\n')
        return

    mem_str = get_resource_str(this_param_dict, 'mem')
    runtime_str = get_resource_str(this_param_dict, 'runtime')

    # run_script = (
    #     "snk_stack_to_ome_zarr "
    #     f"{this_param_dict['input_dirpath']} "
    #     f"{output_dirpath} "
    #     f"-out-oz-fn input-raw.ome.zarr' "
    #     f"--resolution {' '.join([str(x) for x in this_param_dict['resolution']])} "
    #     f"--unit {this_param_dict['unit']} "
    #     f"--downsample_type {this_param_dict['downsample_type'] if 'downsample_type' in this_param_dict else 'Sample'} "
    #     f"--downsample_factors {' '.join(this_param_dict['downsample_factors']) if 'downsample_factors' in this_param_dict else '2 2 2 2 2'} "
    #     f"--cores {this_param_dict['cores'] if 'cores' in this_param_dict else 2048} "
    #     f"--batch_size {this_param_dict['batch_size'] if 'batch_size' in this_param_dict else 32} "
    #     f"--max_cores_per_task {this_param_dict['max_cores_per_task'] if 'max_cores_per_task' in this_param_dict else this_param_dict['batch_size'] if 'batch_size' in this_param_dict else min(this_param_dict['cores'], 32)} "
    #     f"{'--mem {}'.format(mem_str) if 'mem' in this_param_dict else ''} "
    #     f"{'--runtime {}'.format(runtime_str) if 'runtime' in this_param_dict else ''} "
    #     f"{'--cluster slurm' if 'cluster' in this_param_dict and this_param_dict['cluster'] == 'slurm' else ''} "
    #     f"{'-v' if verbose else ''} "
    #     "--continue_run"
    # )
    run_script = [
        "snk_stack_to_ome_zarr ",
        f"{this_param_dict['input_dirpath']} ",
        f"{output_dirpath} ",
        f"-out-oz-fn input-raw.ome.zarr' ",
        f"--resolution {' '.join([str(x) for x in this_param_dict['resolution']])} ",
        f"--unit {this_param_dict['unit']} ",
        f"--downsample_type {this_param_dict['downsample_type'] if 'downsample_type' in this_param_dict else 'Sample'} ",
        f"--downsample_factors {' '.join(this_param_dict['downsample_factors']) if 'downsample_factors' in this_param_dict else '2 2 2 2 2'} ",
        f"--cores {this_param_dict['cores'] if 'cores' in this_param_dict else 2048} ",
        f"--batch_size {this_param_dict['batch_size'] if 'batch_size' in this_param_dict else 32} ",
        f"--max_cores_per_task {this_param_dict['max_cores_per_task'] if 'max_cores_per_task' in this_param_dict else this_param_dict['batch_size'] if 'batch_size' in this_param_dict else min(this_param_dict['cores'], 32)} ",
        f"{'--mem {}'.format(mem_str) if 'mem' in this_param_dict else ''} ",
        f"{'--runtime {}'.format(runtime_str) if 'runtime' in this_param_dict else ''} ",
        f"{'--cluster slurm' if 'cluster' in this_param_dict and this_param_dict['cluster'] == 'slurm' else ''} ",
        f"{'-v' if verbose else ''} ",
        "--continue_run"
    ]

    print(run_script)

    error = run_snakemake_workflow(run_script, parameter_key)

    if not error:
        open(done_fp, 'w').close()


def run_nsbs_alignment(
        parameter_dict,
        parameter_key,
        input_dirpath=None,
        verbose=False
):

    this_param_dict = get_parameters_for_snakemake_workflow(parameter_dict, parameter_key, verbose=verbose)
    if verbose:
        print(f'parameter_key = {parameter_key}')
        print(f'parameter_dict = {parameter_dict}')
        print(f'this_param_dict = {this_param_dict}')
    output_dirpath = os.path.join(this_param_dict['output_dirpath'], parameter_key)
    if not os.path.exists(output_dirpath):
        os.mkdir(output_dirpath)

    done_fp = os.path.join(output_dirpath, f'{parameter_key}.done')
    if os.path.exists(done_fp):
        print('_____________________________________________\n')
        print(f'{parameter_key} already computed!')
        print('_____________________________________________\n')
        return

    mem_str = get_resource_str(this_param_dict, 'mem')
    runtime_str = get_resource_str(this_param_dict, 'runtime')

    if input_dirpath is None:
        input_dirpath = f"{this_param_dict['input_dirpath']} "

    run_script = (
        "snk_elastix_stack_alignment "
        f"{input_dirpath} "
        f"{output_dirpath} "
        f"-out-oz-fn nsbs-step-n.ome.zarr "
        f"--z_step {this_param_dict['z_step'] if 'z_step' in this_param_dict else 1} "
        f"{'--apply_final' if 'apply_final' in this_param_dict and this_param_dict['apply_final'] else ''} "
        f"{'--auto_mask' if 'auto_mask' in this_param_dict and this_param_dict['auto_mask'] else ''} "
        f"--downsample_type {this_param_dict['downsample_type'] if 'downsample_type' in this_param_dict else 'Sample'} "
        f"--downsample_factors {' '.join(this_param_dict['downsample_factors']) if 'downsample_factors' in this_param_dict else '2 2 2 2 2'} "
        f"{'--auto_pad' if 'auto_pad' in this_param_dict else ''} "
        f"--preview_downsample_level {this_param_dict['preview_downsample_level'] if 'preview_downsample_level' in this_param_dict else 2} "
        f"{'--cluster slurm' if 'cluster' in this_param_dict and this_param_dict['cluster'] == 'slurm' else ''} "
        f"--cores {this_param_dict['cores'] if 'cores' in this_param_dict else 2048} "
        f"--batch_size {this_param_dict['batch_size'] if 'batch_size' in this_param_dict else 32} "
        f"--max_cores_per_task {this_param_dict['max_cores_per_task'] if 'max_cores_per_task' in this_param_dict else this_param_dict['batch_size'] if 'batch_size' in this_param_dict else min(this_param_dict['cores'], 32)} "
        f"{'--mem {}'.format(mem_str) if 'mem' in this_param_dict else ''} "
        f"{'--runtime {}'.format(runtime_str) if 'runtime' in this_param_dict else ''} "
        f"{'-v' if verbose else ''} "
        "--continue_run"
    )

    print(run_script)

    error = run_snakemake_workflow(run_script, parameter_key)

    if not error:
        open(done_fp, 'w').close()


def run_apply_transformation(
        parameter_dict,
        parameter_key,
        transforms_filepath,
        input_dirpath=None,
        output_filename=None,
        verbose=False
):
    this_param_dict = get_parameters_for_snakemake_workflow(parameter_dict, parameter_key, verbose=verbose)
    output_dirpath = os.path.join(this_param_dict['output_dirpath'], parameter_key)
    if not os.path.exists(output_dirpath):
        os.mkdir(output_dirpath)

    done_fp = os.path.join(output_dirpath, f'{parameter_key}.done')
    if os.path.exists(done_fp):
        print('_____________________________________________\n')
        print(f'{parameter_key} already computed!')
        print('_____________________________________________\n')
        return

    mem_str = get_resource_str(this_param_dict, 'mem')
    runtime_str = get_resource_str(this_param_dict, 'runtime')

    if input_dirpath is None:
        input_dirpath = f"{this_param_dict['input_dirpath']} "

    run_script = (
        "snk_apply_transformation "
        f"{input_dirpath} "
        f"{transforms_filepath} "
        f"{output_dirpath} "
        f"{'-out-oz-fn {}'.format(output_filename) if output_filename is not None else ''} "
        f"{'--cluster slurm' if 'cluster' in this_param_dict and this_param_dict['cluster'] == 'slurm' else ''} "
        f"--cores {this_param_dict['cores'] if 'cores' in this_param_dict else 2048} "
        f"--batch_size {this_param_dict['batch_size'] if 'batch_size' in this_param_dict else 32} "
        f"--max_cores_per_task {this_param_dict['max_cores_per_task'] if 'max_cores_per_task' in this_param_dict else this_param_dict['batch_size'] if 'batch_size' in this_param_dict else min(this_param_dict['cores'], 32)} "
        f"{'--mem {}'.format(mem_str) if 'mem' in this_param_dict else ''} "
        f"{'--runtime {}'.format(runtime_str) if 'runtime' in this_param_dict else ''} "
        f"{'-v' if verbose else ''} "
        "--continue_run"
    )

    print(run_script)

    error = run_snakemake_workflow(run_script, parameter_key)

    if not error:
        open(done_fp, 'w').close()


def run_ome_zarr_to_stack(
        parameter_dict,
        parameter_key,
        input_dirpath=None,
        output_dirname=None,
        verbose=False
):
    this_param_dict = get_parameters_for_snakemake_workflow(parameter_dict, parameter_key, verbose=verbose)
    output_dirpath = os.path.join(this_param_dict['output_dirpath'], parameter_key)
    if not os.path.exists(output_dirpath):
        os.mkdir(output_dirpath)

    done_fp = os.path.join(output_dirpath, f'{parameter_key}.done')
    if os.path.exists(done_fp):
        print('_____________________________________________\n')
        print(f'{parameter_key} already computed!')
        print('_____________________________________________\n')
        return

    mem_str = get_resource_str(this_param_dict, 'mem')
    runtime_str = get_resource_str(this_param_dict, 'runtime')

    if input_dirpath is None:
        input_dirpath = f"{this_param_dict['input_dirpath']} "

    run_script = (
        "snk_ome_zarr_to_stack " 
        f"{input_dirpath} "
        f"{output_dirpath} "
        f"{output_dirname} " 
        f"{'--cluster slurm' if 'cluster' in this_param_dict and this_param_dict['cluster'] == 'slurm' else ''} "
        f"--cores {this_param_dict['cores'] if 'cores' in this_param_dict else 2048} "
        f"--batch_size {this_param_dict['batch_size'] if 'batch_size' in this_param_dict else 32} "
        f"--max_cores_per_task {this_param_dict['max_cores_per_task'] if 'max_cores_per_task' in this_param_dict else this_param_dict['batch_size'] if 'batch_size' in this_param_dict else min(this_param_dict['cores'], 32)} "
        f"{'--mem {}'.format(mem_str) if 'mem' in this_param_dict else ''} "
        f"{'--runtime {}'.format(runtime_str) if 'runtime' in this_param_dict else ''} "
        f"{'-v' if verbose else ''} "
        "--continue_run"
    )

    print(run_script)

    error = run_snakemake_workflow(run_script, parameter_key)

    if not error:
        open(done_fp, 'w').close()


if __name__ == '__main__':

    batch_idx = int(snakemake.wildcards['idx'])
    output = snakemake.output[0]
    run_info = snakemake.params['run_info']
    n_threads = snakemake.threads
    verbose = run_info['verbose']

    print(f'batch_idx = {batch_idx}')
    print(f'output = {output}')
    print(f'run_info = {run_info}')
    print(f'n_threads = {n_threads}')

    z_range = [batch_idx, batch_idx + run_info['batch_size']]

    print(f'z_range = {z_range}')

    from squirrel.workflows.template_matching import template_matching_stack_alignment_workflow

    template_matching_stack_alignment_workflow(
        run_info['input_ome_zarr_filepath'],
        output,
        **snakemake.params.template_matching_stack_alignment_workflow_params,
        z_range=z_range,
        determine_bounds=snakemake.params['determine_bounds'],
        n_threads=n_threads,
        verbose=verbose
    )

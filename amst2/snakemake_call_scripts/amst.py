
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

    from squirrel.workflows.amst import amst_workflow

    amst_workflow(
        run_info['input_ome_zarr_filepath'],
        output,
        pre_align_key=run_info['input_key'],
        transform=run_info['transform'],
        auto_mask_off=run_info['auto_mask_off'],
        median_radius=run_info['median_radius'],
        z_smooth_method=run_info['z_smooth_method'],
        gaussian_sigma=run_info['gaussian_sigma'],
        elastix_parameters=run_info['elastix_parameter_file'],
        z_range=z_range,
        quiet=False,
        n_workers=run_info['max_cores_per_task'],
        verbose=verbose
    )


if __name__ == '__main__':

    batch_idx = int(snakemake.wildcards['idx'])
    output = snakemake.output[0]
    run_info = snakemake.params['run_info']
    n_threads = snakemake.threads

    print(f'batch_idx = {batch_idx}')
    print(f'output = {output}')
    print(f'run_info = {run_info}')
    print(f'n_threads = {n_threads}')

    z_range = [batch_idx, batch_idx + run_info['batch_size']]

    from squirrel.workflows.convert import stack_to_ome_zarr_workflow
    stack_to_ome_zarr_workflow(
        run_info['stack_path'],
        run_info['output_ome_zarr_filepath'],
        stack_pattern=run_info['stack_pattern'],
        stack_key=run_info['stack_key'],
        resolution=run_info['resolution'],
        unit=run_info['unit'],
        downsample_type=run_info['downsample_type'],
        downsample_factors=run_info['downsample_factors'],
        name=run_info['name'],
        chunk_size=run_info['chunk_size'],
        z_range=z_range,
        xy_range=run_info['crop_xy'],
        save_bounds=run_info['save_bounds'],
        append=True,
        n_threads=n_threads,
        verbose=run_info['verbose'],
    )

    # The dummy output required by snakemake
    open(output, mode='w').close()

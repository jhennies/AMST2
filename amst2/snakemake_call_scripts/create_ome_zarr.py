
if __name__ == '__main__':

    output = snakemake.output[0]
    run_info = snakemake.params['run_info']
    n_threads = snakemake.threads

    print(f'output = {output}')
    print(f'run_info = {run_info}')
    print(f'n_threads = {n_threads}')

    from squirrel.library.ome_zarr import create_ome_zarr
    from squirrel.library.io import load_data_handle, load_data_from_handle_stack

    dtype = load_data_from_handle_stack(
        load_data_handle(run_info['stack_path'], run_info['stack_key'], run_info['stack_pattern'])[0],
        0
    )[0].dtype

    create_ome_zarr(
        run_info['output_ome_zarr_filepath'],
        shape=run_info['stack_shape'],
        resolution=run_info['resolution'],
        unit=run_info['unit'],
        downsample_type=run_info['downsample_type'],
        downsample_factors=run_info['downsample_factors'],
        chunk_size=run_info['chunk_size'],
        dtype=dtype,
        name=run_info['name'],
    )

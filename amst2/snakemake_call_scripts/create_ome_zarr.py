
if __name__ == '__main__':

    output = snakemake.output[0]
    run_info = snakemake.params['run_info']
    n_threads = snakemake.threads

    print(f'output = {output}')
    print(f'run_info = {run_info}')
    print(f'n_threads = {n_threads}')

    from squirrel.library.ome_zarr import create_ome_zarr

    create_ome_zarr(
        run_info['ome_zarr_filepath'],
        shape=run_info['stack_shape'],
        resolution=run_info['resolution'],
        unit=run_info['unit'],
        downsample_type=run_info['downsample_type'],
        chunk_size=run_info['chunk_size'],
        name=run_info['name'],
    )

    # # The dummy output required by snakemake
    # open(output, mode='w').close()

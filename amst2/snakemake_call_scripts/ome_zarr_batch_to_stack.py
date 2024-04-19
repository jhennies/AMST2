
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

    from squirrel.workflows.convert import ome_zarr_to_stack_workflow
    ome_zarr_to_stack_workflow(
        run_info['ome_zarr_filepath'],
        run_info['tif_stack_dirpath'],
        ome_zarr_key=run_info['ome_zarr_key'],
        z_range=z_range,
        n_threads=n_threads,
        verbose=run_info['verbose']
    )

    # The dummy output required by snakemake
    open(output, mode='w').close()

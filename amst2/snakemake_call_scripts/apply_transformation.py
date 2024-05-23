import numpy as np


if __name__ == '__main__':

    batch_idx = int(snakemake.wildcards['idx'])
    output_ome_zarr_filepath = snakemake.input[0]
    transformations_filepath = snakemake.input[1]

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

    from squirrel.library.ome_zarr import get_ome_zarr_handle
    input_ome_zarr_filepath = run_info['input_ome_zarr_filepath']
    input_ome_zarr_dataseth = get_ome_zarr_handle(input_ome_zarr_filepath, key='s0', mode='r')

    from squirrel.library.affine_matrices import AffineStack
    transforms = AffineStack(filepath=transformations_filepath)

    stack_shape = input_ome_zarr_dataseth.shape
    if transforms.exists_meta('stack_shape'):
        stack_shape = np.array(transforms.get_meta('stack_shape'))
    print(f'stack_shape = {stack_shape}')
    print(f'output_shape = {get_ome_zarr_handle(output_ome_zarr_filepath, key="s0", mode="r").shape}')

    # Serialize and apply the transformations
    from squirrel.library.transformation import apply_stack_alignment
    from squirrel.library.affine_matrices import AffineStack
    print(f'z_range = {z_range}')
    print(f'transforms = {len(transforms)}')
    print(f'transforms[z_range[0]: z_range[1]] = {transforms[z_range[0]: z_range[1]]}')
    result_stack = apply_stack_alignment(
        input_ome_zarr_dataseth,
        stack_shape,
        AffineStack(stack=transforms[z_range[0]: z_range[1]], is_sequenced=transforms.is_sequenced),
        z_range=z_range,
        n_workers=n_threads,
        verbose=verbose
    )

    from squirrel.library.ome_zarr import chunk_to_ome_zarr
    chunk_to_ome_zarr(
        result_stack,
        [batch_idx, 0, 0],
        get_ome_zarr_handle(output_ome_zarr_filepath, mode='a'),
        key='s0',
        populate_downsample_layers=True,
        verbose=verbose
    )

    open(output, 'w').close()

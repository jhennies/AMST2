
if __name__ == '__main__':

    input = snakemake.input
    output = snakemake.output[0]
    run_info = snakemake.params['run_info']
    n_threads = snakemake.threads
    preview_downsample_level = run_info['preview_downsample_level']
    input_ome_zarr_filepath = run_info['input_ome_zarr_filepath']
    verbose = run_info['verbose']

    if verbose:
        print(f'input = {input}')
        print(f'output = {output}')

    # Collect the transformations
    import numpy as np

    from squirrel.library.affine_matrices import load_affine_stack_from_multiple_files
    transforms = load_affine_stack_from_multiple_files(input)

    # Downsample the transformations
    from squirrel.library.ome_zarr import get_ome_zarr_handle, get_scale_of_downsample_level

    input_ome_zarr_fileh = get_ome_zarr_handle(input_ome_zarr_filepath, mode='r')
    scale_full = get_scale_of_downsample_level(input_ome_zarr_fileh, 0)
    scale_ds = get_scale_of_downsample_level(input_ome_zarr_fileh, preview_downsample_level)
    scale = (np.array(scale_ds) / np.array(scale_full)).astype(int)
    assert scale[0] == scale[1] == scale[2], 'Implemented only for isotropic scaling!'
    scale = 1 / scale[0]
    if verbose:
        print(f'scale = {scale}')
    if not transforms.is_sequenced:
        transforms = transforms.get_sequenced_stack()
    transforms = transforms.get_scaled(scale)

    # Serialize and apply the transformations
    from squirrel.library.transformation import apply_stack_alignment
    input_ome_zarr_dataseth = input_ome_zarr_fileh[f's{preview_downsample_level}']
    result_stack = apply_stack_alignment(
        input_ome_zarr_dataseth,
        input_ome_zarr_dataseth.shape,
        transforms,
        no_adding_of_transforms=True,
        verbose=verbose
    )

    # Write the result stack to disk
    from squirrel.library.io import write_h5_container
    write_h5_container(output, result_stack)

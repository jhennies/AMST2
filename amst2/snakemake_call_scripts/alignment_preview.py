
if __name__ == '__main__':

    input = snakemake.input
    output = snakemake.output[0]
    transforms_filepath = None
    if len(output) > 1:
        transforms_filepath = snakemake.output[1]
    run_info = snakemake.params['run_info']
    n_threads = snakemake.threads
    preview_downsample_level = run_info['preview_downsample_level']
    input_ome_zarr_filepath = run_info['input_ome_zarr_filepath']
    save_joined_transforms = snakemake.params['save_joined_transforms']
    compute_auto_pad = snakemake.params['compute_auto_pad']
    verbose = run_info['verbose']
    if save_joined_transforms:
        assert transforms_filepath is not None

    if verbose:
        print(f'input = {input}')
        print(f'output = {output}')

    # Collect the transformations
    import numpy as np

    from squirrel.library.affine_matrices import load_affine_stack_from_multiple_files
    transforms = load_affine_stack_from_multiple_files(input, sequence_stack=False)
    print(f'meta = {transforms.get_meta()}')
    print(f'len(transforms) = {len(transforms)}')

    if save_joined_transforms:
        transforms.to_file(transforms_filepath)

    transforms = transforms.get_sequenced_stack()

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
    # if not transforms.is_sequenced:
    #     transforms = transforms.get_sequenced_stack()
    assert transforms.is_sequenced

    # Perform auto-pad
    stack_shape = None
    if compute_auto_pad:
        from squirrel.library.image import apply_auto_pad
        transforms, stack_shape = apply_auto_pad(
            transforms,
            [len(transforms), 0., 0.],
            transforms.get_meta('bounds'),
            extra_padding=16
        )
        transforms.set_meta('stack_shape', stack_shape)

    if stack_shape is None:
        if transforms.exists_meta('stack_shape'):
            stack_shape = transforms.get_meta('stack_shape')
        else:
            stack_shape = input_ome_zarr_fileh['s0'].shape
    if not transforms.exists_meta('stack_shape'):
        transforms.set_meta('stack_shape', stack_shape)

    stack_shape = (np.array(stack_shape) * scale).astype(int)

    print(f'len(transforms) = {len(transforms)}')
    if transforms.exists_meta('z_step'):
        transforms = transforms.apply_z_step()
    print(f'len(transforms) = {len(transforms)}')
    print(f'meta = {transforms.get_meta()}')

    # Scale the transformations
    transforms = transforms.get_scaled(scale)

    print(f'stack_shape = {stack_shape}')
    print(f'scale = {scale}')

    # Apply the transformations
    from squirrel.library.transformation import apply_stack_alignment
    input_ome_zarr_dataseth = input_ome_zarr_fileh[f's{preview_downsample_level}']

    result_stack = apply_stack_alignment(
        input_ome_zarr_dataseth,
        stack_shape,
        transforms,
        no_adding_of_transforms=True,
        verbose=verbose
    )

    # Write the result stack to disk
    from squirrel.library.io import write_h5_container
    write_h5_container(output, result_stack)


if __name__ == '__main__':

    input = snakemake.input
    output = snakemake.output[0]
    transforms_filepath = None
    if len(snakemake.output) > 1:
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
    if not transforms.is_sequenced:
        transforms = transforms.get_sequenced_stack()
    print(f'meta = {transforms.get_meta()}')
    print(f'len(transforms) = {len(transforms)}')

    from squirrel.library.io import get_filetype
    input_filetype = get_filetype(input_ome_zarr_filepath)
    if verbose:
        print(f'input_ome_zarr_filepath = {input_ome_zarr_filepath}')
        print(f'input_filetype = {input_filetype}')
    if input_filetype == 'ome_zarr':
        from squirrel.library.ome_zarr import get_ome_zarr_handle
        input_fileh = get_ome_zarr_handle(input_ome_zarr_filepath, mode='r')
    else:
        from squirrel.library.io import load_data_handle
        input_fileh, _ = load_data_handle(input_ome_zarr_filepath, key=run_info['stack_key'], pattern=run_info['stack_pattern'])
    # import zarr
    # try:
    #     input_ome_zarr_fileh = get_ome_zarr_handle(input_ome_zarr_filepath, mode='r')
    # except zarr.errors.PathNotFoundError:
    #     input_ome_zarr_fileh = None

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
            if input_filetype == 'ome_zarr':
                stack_shape = input_fileh[run_info['stack_key']].shape
            else:
                stack_shape = input_fileh.shape

    if not transforms.exists_meta('stack_shape'):
        transforms.set_meta('stack_shape', stack_shape)

    print(f'len(transforms) = {len(transforms)}')
    if transforms.exists_meta('z_step'):
        transforms = transforms.apply_z_step()
    print(f'len(transforms) = {len(transforms)}')
    print(f'meta = {transforms.get_meta()}')

    if save_joined_transforms:
        transforms.to_file(transforms_filepath)

    if input_filetype == 'ome_zarr':

        from squirrel.library.ome_zarr import get_scale_of_downsample_level

        scale_full = get_scale_of_downsample_level(input_fileh, 0)
        scale_ds = get_scale_of_downsample_level(input_fileh, preview_downsample_level)
        scale = (np.array(scale_ds) / np.array(scale_full)).astype(int)
        assert scale[0] == scale[1] == scale[2], 'Implemented only for isotropic scaling!'
        scale = 1 / scale[0]
        if verbose:
            print(f'scale = {scale}')

        # Scale the transformations
        stack_shape = (np.array(stack_shape) * scale).astype(int)
        transforms = transforms.get_scaled(scale)
        print(f'stack_shape = {stack_shape}')
        print(f'scale = {scale}')

        # Apply the transformations
        from squirrel.library.transformation import apply_stack_alignment
        input_ome_zarr_dataseth = input_fileh[f's{preview_downsample_level}']

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

    else:
        open(output, mode='w').close()

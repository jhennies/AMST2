
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
    import json
    import numpy as np
    from squirrel.library.transformation import load_transform_matrices
    transforms = []
    for transform_filepath in input:
        transforms.extend(
            np.array(load_transform_matrices(transform_filepath, validate=True, ndim=2))
        )
        # with open(transform_filepath, mode='r') as f:
        #     transforms.extend(json.load(f))

    # Downsample the transformations
    from squirrel.library.transformation import scale_sequential_affines
    from squirrel.library.ome_zarr import get_ome_zarr_handle, get_scale_of_downsample_level

    input_ome_zarr_fileh = get_ome_zarr_handle(input_ome_zarr_filepath, mode='r')
    scale = get_scale_of_downsample_level(input_ome_zarr_fileh, preview_downsample_level)
    assert scale[0] == scale[1] == scale[2], 'Implemented only for isotropic scaling!'
    scale = 1 / scale[0]
    # if verbose:
    print(f'transforms[4] = {transforms[0]}')
    print(f'scale = {scale}')
    transforms = scale_sequential_affines(transforms, scale, xy_pivot=(0., 0.))
    # if verbose:
    print(f'transforms[1] = {transforms[0]}')

    # Serialize and apply the transformations
    from squirrel.library.transformation import apply_stack_alignment
    input_ome_zarr_dataseth = input_ome_zarr_fileh[f's{preview_downsample_level}']
    result_stack = apply_stack_alignment(
        input_ome_zarr_dataseth,
        input_ome_zarr_dataseth.shape,
        transforms,
        no_adding_of_transforms=False,
        xy_pivot=(0., 0.),
        param_order='M',
        verbose=verbose
    )

    # Write the result stack to disk
    from squirrel.library.io import write_h5_container
    write_h5_container(output, result_stack)

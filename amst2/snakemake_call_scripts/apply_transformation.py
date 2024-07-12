import numpy as np


if __name__ == '__main__':

    batch_idx = int(snakemake.wildcards['idx'])
    output_ome_zarr_filepath = snakemake.input[0]
    transformations_filepath = snakemake.input[1:]

    output = snakemake.output[0]
    run_info = snakemake.params['run_info']
    n_threads = snakemake.threads
    verbose = run_info['verbose']

    input_ome_zarr_filepath = run_info['input_ome_zarr_filepath']

    print(f'batch_idx = {batch_idx}')
    print(f'output = {output}')
    print(f'run_info = {run_info}')
    print(f'n_threads = {n_threads}')

    z_range = [batch_idx, batch_idx + run_info['batch_size']]

    print(f'z_range = {z_range}')

    # from squirrel.library.ome_zarr import get_ome_zarr_handle
    # input_ome_zarr_dataseth = get_ome_zarr_handle(input_ome_zarr_filepath, key='s0', mode='r')
    # stack_shape = input_ome_zarr_dataseth.shape
    # image_shape = stack_shape[1:]
    #
    # import os
    # if True:  # os.path.isdir(transformations_filepath[0]):
    #     print(transformations_filepath)
    #     from squirrel.library.elastix import ElastixStack, ElastixMultiStepStack
    #     from squirrel.library.affine_matrices import AffineStack
    #     stacks = []
    #     for transform_path in transformations_filepath:
    #         if os.path.isdir(transform_path):
    #             stacks.append(ElastixStack(dirpath=transform_path))  # , image_shape=target_image_shape))
    #         else:
    #             stack = AffineStack(filepath=transform_path)
    #             if stack.exists_meta('stack_shape'):
    #                 image_shape = stack.get_meta('stack_shape')[1:]
    #             else:
    #                 image_shape = stack_shape[1:]
    #             if not stack.is_sequenced:
    #                 stack = stack.get_sequenced_stack()
    #             stacks.append(ElastixStack(stack=stack, image_shape=image_shape))
    #
    #     emss = ElastixMultiStepStack(stacks=stacks, image_shape=image_shape)
    #
    #     result_stack = emss.apply_on_image_stack(
    #         input_ome_zarr_dataseth,
    #         target_image_shape=image_shape,
    #         z_range=z_range,
    #         n_workers=n_threads,
    #         verbose=verbose
    #     )
    #
    # else:
    #     assert len(transformations_filepath) == 0
    #     transformations_filepath = transformations_filepath[0]
    #     from squirrel.library.affine_matrices import AffineStack
    #     transforms = AffineStack(filepath=transformations_filepath)
    #
    #     if transforms.exists_meta('stack_shape'):
    #         stack_shape = np.array(transforms.get_meta('stack_shape'))
    #     print(f'stack_shape = {stack_shape}')
    #     print(f'output_shape = {get_ome_zarr_handle(output_ome_zarr_filepath, key="s0", mode="r").shape}')
    #
    #     # Serialize and apply the transformations
    #     from squirrel.library.transformation import apply_stack_alignment
    #     from squirrel.library.affine_matrices import AffineStack
    #     print(f'z_range = {z_range}')
    #     print(f'transforms = {len(transforms)}')
    #     print(f'transforms[z_range[0]: z_range[1]] = {transforms[z_range[0]: z_range[1]]}')
    #     result_stack = apply_stack_alignment(
    #         input_ome_zarr_dataseth,
    #         stack_shape,
    #         AffineStack(stack=transforms[z_range[0]: z_range[1]], is_sequenced=transforms.is_sequenced),
    #         z_range=z_range,
    #         n_workers=n_threads,
    #         verbose=verbose
    #     )

    from squirrel.workflows.elastix import apply_multi_step_stack_alignment_workflow

    result_stack = apply_multi_step_stack_alignment_workflow(
        input_ome_zarr_filepath,
        transformations_filepath,
        key='s0',
        auto_pad=False,
        z_range=z_range,
        start_transform_id=z_range[0],
        n_workers=1,  # This is using elastix which is parallelized internally!
        quiet=False,
        verbose=verbose
    )

    from squirrel.library.ome_zarr import get_ome_zarr_handle

    print(f'result_stack.shape = {result_stack.shape}')
    print(f'output dataset shape = {get_ome_zarr_handle(output_ome_zarr_filepath, mode="r")["s0"].shape}')
    print(f'input dataset shape = {get_ome_zarr_handle(input_ome_zarr_filepath, mode="r")["s0"].shape}')
    from squirrel.library.volume import pad_volume
    result_stack = pad_volume(result_stack, get_ome_zarr_handle(output_ome_zarr_filepath, mode="r")["s0"].shape, axes=[1, 2])
    print(f'output dataset shape (adjusted) = {get_ome_zarr_handle(output_ome_zarr_filepath, mode="r")["s0"].shape}')
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

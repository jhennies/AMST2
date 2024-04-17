
if __name__ == '__main__':

    sn_input = snakemake.input
    sn_output = snakemake.output[0]

    run_info = snakemake.params['run_info']
    use_tm = snakemake.params['use_tm']
    use_local = snakemake.params['use_local']

    from squirrel.library.affine_matrices import load_affine_stack_from_multiple_files

    if use_tm and use_local:

        local_filepaths = sn_input[:int(len(sn_input) / 2)]
        tm_filepaths = sn_input[int(len(sn_input) / 2):]

        local_transforms = load_affine_stack_from_multiple_files(local_filepaths)
        tm_transforms = load_affine_stack_from_multiple_files(tm_filepaths)

        assert tm_transforms.is_sequenced, 'Template matching sequences are always sequenced, this one must be too!'
        if not local_transforms.is_sequenced:
            local_transforms = local_transforms.get_sequenced_stack()

        tm_local = tm_transforms * (-local_transforms)
        final_transforms = local_transforms * tm_local.get_smoothed_stack(run_info['combine_sigma'])

        final_transforms.set_meta('bounds', tm_transforms.get_meta('bounds'))

    elif use_tm and not use_local:

        tm_filepaths = sn_input
        tm_transforms = load_affine_stack_from_multiple_files(tm_filepaths)

        assert tm_transforms.is_sequenced, 'Template matching sequences are always sequenced, this one must be too!'

        final_transforms = tm_transforms
        final_transforms.set_meta('bounds', tm_transforms.get_meta('bounds'))

    elif not use_local and use_local:

        local_filepaths = sn_input
        local_transforms = load_affine_stack_from_multiple_files(local_filepaths)

        if not local_transforms.is_sequenced:
            local_transforms = local_transforms.get_sequenced_stack()

        final_transforms = local_transforms
        final_transforms.set_meta('bounds', local_transforms.get_meta('bounds'))

    else:
        raise ValueError('Either template matching or local alignment must be active!')

    print(f'bounds = {final_transforms.get_meta("bounds")}')
    from squirrel.library.image import apply_auto_pad
    final_transforms, stack_shape = apply_auto_pad(
        final_transforms,
        [len(final_transforms), 0., 0.],
        final_transforms.get_meta('bounds'),
        extra_padding=16
    )

    final_transforms.set_meta('stack_shape', stack_shape)
    final_transforms.to_file(sn_output)


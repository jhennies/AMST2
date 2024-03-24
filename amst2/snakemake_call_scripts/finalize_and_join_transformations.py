
import numpy as np


if __name__ == '__main__':

    sn_input = snakemake.input
    sn_output = snakemake.output[0]

    run_info = snakemake.params['run_info']

    local_filepaths = sn_input[:int(len(sn_input) / 2)]
    tm_filepaths = sn_input[int(len(sn_input) / 2):]

    from squirrel.library.transformation import (
        load_transform_matrices_from_multiple_files,
        smooth_2d_affine_sequence,
        serialize_affine_sequence,
        save_transformation_matrices
    )
    from squirrel.library.linalg import dot_product_of_sequences
    from squirrel.library.elastix import save_transforms

    local_transforms, local_sequenced = load_transform_matrices_from_multiple_files(local_filepaths, validate=True, ndim=2)
    tm_transforms, tm_sequenced = load_transform_matrices_from_multiple_files(tm_filepaths, validate=True, ndim=2)

    assert tm_sequenced, 'Template matching sequences are always sequenced, this one must be too!'
    if not local_sequenced:
        local_transforms = serialize_affine_sequence(local_transforms, param_order='M', out_param_order='M')

    tm_local = dot_product_of_sequences(tm_transforms, local_transforms, inverse=(0, 1))
    tm_local_smoothed = smooth_2d_affine_sequence(tm_local, run_info['combine_sigma'])
    final_transforms = dot_product_of_sequences(local_transforms, tm_local_smoothed)

    final_transforms = np.array(save_transforms(final_transforms, None, param_order='M', save_order='C', ndim=2))[:, :6].tolist()

    save_transformation_matrices(sn_output, final_transforms, sequenced=True)

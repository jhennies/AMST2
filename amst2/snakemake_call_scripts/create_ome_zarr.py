import json
import numpy as np

if __name__ == '__main__':

    output = snakemake.output[0]
    params = snakemake.params
    output_ome_zarr_filepath = params['output_ome_zarr_filepath']
    stack_shape = params['stack_shape']
    resolution = params['resolution']
    unit = params['unit']
    downsample_type = params['downsample_type']
    downsample_factors = params['downsample_factors']
    chunk_size = params['chunk_size']
    dtype = np.dtype(params['dtype']).name
    name = params['name']
    n_threads = snakemake.threads

    print(f'output = {output}')
    print(f'n_threads = {n_threads}')

    if isinstance(stack_shape, str):
        stack_shape_split = str.split(stack_shape, ':')
        print(f'stack_shape_split = {stack_shape_split}')
        with open(stack_shape_split[0], mode='r') as f:
            stack_shape = json.load(f)
        for item in stack_shape_split[1:]:
            stack_shape = stack_shape[item]
        print(f'stack_shape = {stack_shape}')

    from squirrel.library.ome_zarr import OMEZarrStore
    OMEZarrStore.create(
        path=output_ome_zarr_filepath,
        shape=stack_shape,
        dtype=dtype,
        chunks=chunk_size,
        shards=None,
        downsample_factors=downsample_factors,
        resolution=resolution,
        unit=unit,
        downsample_method=downsample_type,
        ome_version='0.4',
        zarr_format=2,
        overwrite=False
    )

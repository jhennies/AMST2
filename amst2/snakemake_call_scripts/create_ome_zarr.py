import json

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
    dtype = params['dtype']
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

    from squirrel.library.ome_zarr import create_ome_zarr

    create_ome_zarr(
        output_ome_zarr_filepath,
        shape=stack_shape,
        resolution=resolution,
        unit=unit,
        downsample_type=downsample_type,
        downsample_factors=downsample_factors,
        chunk_size=chunk_size,
        dtype=dtype,
        name=name,
    )


if __name__ == '__main__':

    sn_input = snakemake.input
    sn_output = snakemake.output[0]

    run_info = snakemake.params['run_info']

    from squirrel.library.elastix import load_transform_stack_from_multiple_files

    transforms = load_transform_stack_from_multiple_files(sn_input)

    transforms.to_file(sn_output)

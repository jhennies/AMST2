
def make_cache_folder_structure(target_dirpath, cache_name, continue_run=False):

    import os

    target_dirpath = os.path.abspath(target_dirpath)

    if not os.path.exists(target_dirpath):
        os.mkdir(target_dirpath)
    cache_dirpath = os.path.join(target_dirpath, 'snk_cache')
    if not os.path.exists(cache_dirpath):
        os.mkdir(cache_dirpath)
    this_cache_dirpath = os.path.join(cache_dirpath, cache_name)
    if not os.path.exists(this_cache_dirpath):
        os.mkdir(this_cache_dirpath)
    elif os.path.exists(this_cache_dirpath) and not continue_run:
        raise RuntimeError('Cache directory exists. If you want to continue use --continue_run')

    return cache_dirpath, this_cache_dirpath


def solve_output_path_and_name(output_path, output_name, target_dirpath, fallback_name):

    import os

    if output_path is None:
        if output_name is not None:
            output_path = os.path.join(target_dirpath, output_name)
        else:
            output_path = os.path.join(target_dirpath, fallback_name)
    if output_name is None:
        output_name = os.path.split(output_path)[1]

    output_path = os.path.abspath(output_path)

    return output_path, output_name


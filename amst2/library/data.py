
import numpy as np


def estimate_mem_mb(handle):

    from squirrel.library.io import load_data_from_handle_stack
    data, _ = load_data_from_handle_stack(handle, 0)
    size = np.prod(data.shape)

    if data.dtype == 'uint8':
        return size / 1e6
    if data.dtype in ['uint16', 'float16']:
        return size * 2 / 1e6
    if data.dtype in ['uint32', 'float32']:
        return size * 4 / 1e6
    if data.dtype in ['uint64', 'float64']:
        return size * 8 / 1e6
    raise NotImplementedError(f'dtype == {data.dtype} is not implemented!')

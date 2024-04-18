import runpy
from setuptools import setup

version = runpy.run_path("amst2/__version__.py")["__version__"]

setup(
    name='amst2',
    version=version,
    author='Julian Hennies',
    author_email='hennies@embl.de',
    url='https://github.com/jhennies/amst2',
    license="GPLv3",
    packages=['amst2'],
    entry_points={
        'console_scripts': [
            'snk_stack_to_ome_zarr = amst2.conversion:snk_stack_to_ome_zarr'
            'snk_default_amst_pre_alignment = amst2.stack_alignment:snk_default_amst_pre_alignment'
        ]
    },
    install_requires=[
        'numpy',
        'h5py',
        'tifffile'
    ]
)

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
            'snk_stack_to_ome_zarr = amst2.data:snk_stack_to_ome_zarr'
        ]
    },
    install_requires=[
        'numpy',
        'h5py',
        'tifffile'
    ]
)

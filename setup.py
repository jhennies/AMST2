import runpy
from setuptools import setup

version = runpy.run_path("squirrel/__version__.py")["__version__"]

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
            'apply_stack_alignment = squirrel.apply_transformation:apply_stack_alignment'
        ]
    },
    install_requires=[
        'numpy',
        'h5py',
        'tifffile'
    ]
)

import runpy
from setuptools import setup, find_packages

version = runpy.run_path("amst2/__version__.py")["__version__"]

setup(
    name='amst2',
    version=version,
    author='Julian Hennies',
    author_email='hennies@embl.de',
    url='https://github.com/jhennies/amst2',
    license="GPLv3",
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'snk_stack_to_ome_zarr = amst2.conversion:snk_stack_to_ome_zarr',
            'snk_ome_zarr_to_stack = amst2.conversion:snk_ome_zarr_to_stack',
            'snk_default_amst_pre_alignment = amst2.stack_alignment:snk_default_amst_pre_alignment',
            'snk_amst = amst2.amst:snk_amst',
            'snk_elastix_stack_alignment = amst2.stack_alignment:snk_elastix_stack_alignment',
            'snk_template_matching_stack_alignment = amst2.stack_alignment:snk_template_matching_stack_alignment',
            'snk_apply_transformation = amst2.stack_alignment:snk_apply_transformation',
            'snk_normalize_stack = amst2.stack_operations:snk_normalize_stack',
            'amst-nsbs-pre-align = amst2.workflows.pre_align:nsbs_pre_align',
            'amst-cleanup-nsbs-pre-align = amst2.workflows.pre_align:cleanup_nsbs_pre_align',
            'amst-run = amst2.workflows.amst:amst',
            'amst-cleanup-run = amst2.workflows.amst:cleanup_amst'
        ]
    },
    install_requires=[
        'numpy',
        'h5py',
        'tifffile',
        'squirrel>=0.1.2'
    ]
)

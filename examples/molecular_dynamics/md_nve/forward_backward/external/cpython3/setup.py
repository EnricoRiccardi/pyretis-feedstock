# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""A script for building the C extension."""
from setuptools import (
    setup,
    Extension,
)
import numpy as np


LJMODULE = Extension(
    'ljc',
    sources=['ljc.c'],
    extra_compile_args=['-Ofast', '-march=native'],
)
setup(
    name='PyRETIS Lennard-Jones c extension',
    description='C extension for the Lennard-Jones potential.',
    ext_modules=[LJMODULE],
    include_dirs=[np.get_include()],
)

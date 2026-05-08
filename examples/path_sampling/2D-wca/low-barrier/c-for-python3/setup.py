# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""A script for building the WCA C extension."""
from setuptools import (
    setup,
    Extension,
)
import numpy as np


WCAFORCE = Extension(
    'wcaforces',
    sources=['wcaforces.c'],
    extra_compile_args=['-Ofast', '-march=native']
)
setup(
    name='PyRETIS WCA C extension',
    description='C extension for WCA potential.',
    ext_modules=[WCAFORCE],
    include_dirs=[np.get_include()],
)

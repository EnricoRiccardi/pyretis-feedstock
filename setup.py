# -*- coding: utf-8 -*-
# Copyright (c) 2023, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""PyRETIS - A simulation package for rare event simulations.

PyRETIS - A simulation package for rare event simulations.
Copyright (c) 2023, PyRETIS Development Team

This file is part of PyRETIS.

PyRETIS is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

PyRETIS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with PyRETIS. If not, see <http://www.gnu.org/licenses/>
"""
import ast
from codecs import open as openc
import pathlib
from setuptools import setup, find_namespace_packages

FULL_VERSION = '3.0.1.dev0'  # Automatically set by setup_version.py


def get_long_description():
    """Return the contents of README.rst as a string."""
    here = pathlib.Path(__file__).absolute().parent
    long_description = ''
    with openc(here.joinpath('README.rst'), encoding='utf-8') as fileh:
        long_description = fileh.read()
    return long_description


def get_version():
    """Return the version from version.py as a string."""
    here = pathlib.Path(__file__).absolute().parent
    filename = here.joinpath('pyretis', 'version.py')
    with openc(filename, encoding='utf-8') as fileh:
        for lines in fileh:
            if lines.startswith('FULL_VERSION ='):
                version = ast.literal_eval(lines.split('=')[1].strip())
                return version
    return FULL_VERSION


def get_requirements(docs=False):
    """Read requirements.txt and return a list of requirements."""
    here = pathlib.Path(__file__).absolute().parent
    requirements = []
    doc_package = ('sphinx', 'sphinx_bootstrap_theme')
    filename = here.joinpath('requirements.txt')
    with openc(filename, encoding='utf-8') as fileh:
        for lines in fileh:
            if '>=' in lines:
                package = lines.split('>=')[1].strip()
            elif '<' in lines:
                package = lines.split('<')[1].strip()
            if not docs and package in doc_package:
                continue
            requirements.append(lines.strip())
    return requirements


setup(
    name='pyretis',
    version=get_version(),
    description='A simulation package for rare events',
    long_description=get_long_description(),
    url='http://www.pyretis.org',
    author='PyRETIS Development Team',
    author_email='pyretis@pyretis.org',
    license='LGPLv2.1+',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        ('License :: OSI Approved :: '
         'GNU Lesser General Public License v2 or later (LGPLv2+)'),
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    keywords='rare-events md mc tps simulation tis retis',
    packages=find_namespace_packages(exclude=['docs']),
    package_data={
        'pyretis': [
            'pyretis.mplstyle',
            'pyretis/inout/report/templates/*',
            'pyretis/pyvisa/pyvisaGUI.ui'
        ]
    },
    include_package_data=True,
    install_requires=get_requirements(docs=False),
    entry_points={
        'console_scripts': [
            'pyretisrun = pyretis.bin.pyretisrun:entry_point',
            'pyretisanalyse = pyretis.bin.pyretisanalyse:entry_point',
            'pyvisa = pyretis.bin.pyvisa:entry_point'
        ]
    },
)

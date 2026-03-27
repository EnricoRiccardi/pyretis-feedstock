# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the common methods in pyretis.inout.common."""
import os
import logging
import pytest
from pyretis.inout.common import (
    make_dirs,
    simplify_ensemble_name,
    add_dirname,
    name_file,
    generate_file_name,
    create_empty_ensembles
)


logging.disable(logging.CRITICAL)


HERE = os.path.abspath(os.path.dirname(__file__))


def remove_dir(dirname):
    """Remove a directory."""
    try:
        os.removedirs(dirname)
    except OSError:
        pass


class TestMethods:
    """Test some of the methods from pyretis.inout.common."""
    def test_create_empty_ensembles(self):
        """Test that we can properly create empty ensembles."""
        settings = {'simulation': {'interfaces': [1, 2, 3, 4, 5],
                                   'flux': True,
                                   'zero_ensemble': True}}
        create_empty_ensembles(settings)
        settings = {'simulation': {'interfaces': [1, 2, 3, 4, 5],
                                   'flux': True,
                                   'zero_ensemble': True}}
        settings['ensemble'] = [settings]
        settings['ensemble'][0]['interface'] = 1
        settings['ensemble'][0]['tis'] = {'ensemble_number': 123321}
        create_empty_ensembles(settings)
        assert settings['ensemble'][0]['interface'] == 1
        assert settings['ensemble'][0]['tis']['ensemble_number'] == 123321

        settings = {'simulation': {'interfaces': [1, 2, 3, 4, 5],
                                   'flux': True,
                                   'zero_ensemble': True},
                    'tis': {'ensemble_number': 123321}}
        create_empty_ensembles(settings)
        assert settings['ensemble'][0]['tis']['ensemble_number'] == 123321
        assert len(settings['ensemble']) == 1

        settings = {'simulation': {'interfaces': [1, 2, 3, 4, 5],
                                   'flux': True,
                                   'zero_ensemble': True},
                    'ensemble': [{'interface': 3,
                                  'gnappo': 'lappo',
                                  'tis': {'ensemble_number': 3}}]}
        create_empty_ensembles(settings)
        assert settings['ensemble'][3]['interface'] == 3
        assert settings['ensemble'][3]['gnappo'] == 'lappo'
        assert settings['ensemble'][3]['tis']['ensemble_number'] == 3
        # todo. This should be part of the ensamble numering issue.
        # with pytest.raises(ValueError) as err:
        #    create_empty_ensembles(settings)
        # assert 'Invalid entry for setting ensemble',
        #    str(err.value))
        # settings['ensemble'][0] = {'tis': {'ensemble_number': 'flux'}}
        # create_empty_ensembles(settings)
        # assert
        #     settings['simulation']['ensemble'][0]['tis']['ensemble_number'],
        #     'flux')

    def test_make_dirs(self):
        """Test that we can create directories."""
        dirname = os.path.join(HERE, 'testdir')
        make_dirs(dirname)
        assert os.path.isdir(dirname)
        remove_dir(dirname)
        dirname = os.path.join(HERE, 'already_exists')
        with pytest.raises(OSError):
            make_dirs(dirname)
        dirname = os.path.join(HERE, 'dir_exists')
        remove_dir(dirname)
        make_dirs(dirname)
        msg = make_dirs(dirname)
        assert msg.endswith('already exist.')
        remove_dir(dirname)

    def test_simplify_ensemble_name(self):
        """Test that we can simplify ensemble names."""
        cases = [('[0^-]', '000'), ('[0^+]', '001'), ('[1^+]', '002')]
        for case in cases:
            txt = simplify_ensemble_name(case[0], fmt='{:03d}')
            assert txt == case[1]
        txt = simplify_ensemble_name('001', fmt='{:03d}')
        assert txt == '001'
        txt = simplify_ensemble_name('[1]', fmt='{:03d}')
        assert txt == '002'

    def test_add_dirname(self):
        """Test that we can add a directory to a filename."""
        path = add_dirname('filename.txt', 'path')
        assert path == os.path.join('path', 'filename.txt')
        path = add_dirname('filename.txt', None)
        assert path == 'filename.txt'

    def test_name_file(self):
        """Test that we can name a file."""
        name = name_file('test', 'txt', path='path')
        filename = ''.join(['test', os.extsep, 'txt'])
        filepath = os.path.join('path', filename)
        assert name == filepath

    def test_generate_filename(self):
        """Test the generation of file names."""
        settings = {'output': {}}
        name = generate_file_name('base.txt', 'dir', settings)
        assert name == os.path.join('dir', 'base.txt')
        settings = {'output': {'prefix': 'abc-'}}
        name = generate_file_name('base.txt', 'dir', settings)
        assert name == os.path.join('dir', 'abc-base.txt')

    def test_check_python_version(self):
        # DO custom import and go from there to prevent overriding anything
        # important
        # Set the min python version to next mayor version to remember updating
        import pyretis.inout.common as common
        common.MIN_PYTHON_VERSION = (4, 0)
        with pytest.raises(SystemExit):
            common.check_python_version()

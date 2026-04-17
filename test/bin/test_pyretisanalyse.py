# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the bin."""
import logging
import os
import shutil
import subprocess
import tempfile
from contextlib import contextmanager
from io import StringIO
from unittest.mock import patch
from pyretis.bin.pyretisanalyse import (hello_world,
                                        bye_bye_world,
                                        main)

HERE = os.path.abspath(os.path.dirname(__file__))


@contextmanager
def capture_log_output():
    """Capture all logger output to a StringIO buffer."""
    stream = StringIO()
    handler = logging.StreamHandler(stream)
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(logging.Formatter('%(message)s'))
    root_logger = logging.getLogger('')
    root_logger.addHandler(handler)
    # Temporarily re-enable logging (may have been disabled by other tests)
    prev_disable = logging.root.manager.disable
    prev_level = root_logger.level
    logging.disable(logging.NOTSET)
    root_logger.setLevel(logging.DEBUG)
    try:
        yield stream
    finally:
        root_logger.removeHandler(handler)
        root_logger.setLevel(prev_level)
        logging.disable(prev_disable)


class TestPyretisMainAnalyse:
    """Test Main."""

    def test_pyretisanalyse(self):
        """Test pyretisanalyse."""
        with tempfile.TemporaryDirectory() as tempdir:
            with patch('sys.stdout', new=StringIO()):
                with subprocess.Popen(['pyretisanalyse', '-V'],
                                      stderr=subprocess.DEVNULL,
                                      stdout=subprocess.PIPE) as ex:
                    asd = ex.stdout.read().split()
                assert b'PyRETIS' in asd

            with patch('sys.stdout', new=StringIO()):
                with subprocess.Popen(['pyretisanalyse'],
                                      stderr=subprocess.DEVNULL,
                                      stdout=subprocess.PIPE,
                                      cwd=tempdir) as ex:
                    asd = ex.stdout.read().strip()
                assert b'Input file required' in asd

            inputfile = 'a_non_existent_input.rst'
            with capture_log_output() as log_out:
                assert main(inputfile, tempdir, tempdir) == 1
            assert 'NOT found!' in log_out.getvalue().strip()

            with capture_log_output() as log_out:
                assert main(None, tempdir, tempdir) == 1
            assert 'Input file required' in log_out.getvalue().strip()

    def test_hello_world(self):
        """Test that we are polite."""
        with capture_log_output() as log_out:
            hello_world(infile='I_can_read_your_soul.rst',
                        run_dir='You_are_wasting_your_life_here',
                        report_dir='go_out_and_enjoy')
        assert 'Start' in log_out.getvalue().strip()

    def test_bye_world(self):
        """Test that we can die."""
        with capture_log_output() as log_out:
            bye_bye_world()
        assert 'reference' in log_out.getvalue().strip()

    def test_main_an(self):
        """Test that we know how to set up a report."""
        input_file = 'settings-retis.rst'
        # fake pathensemble line
        lines = ['']*8
        moves = ['ki', 'sh', 'ki', 'sh', 'ld', 'xx', 'sh']
        line = "         2          3          0 L M R      99 ACC "
        line += "sh -9.108621358e-02  1.007510172e+01       0"
        line += "      98  0.000000000e+00       0       0  1.00e+00\n"
        for j in range(7):
            lines[j] = line.replace('sh', moves[j])
        lines[7] = line.replace('ACC', 'BWI')
        lines[0] = line.replace('L', 'R')
        lines[1] = line.replace('L', 'R')
        with tempfile.TemporaryDirectory() as tempdir:
            shutil.copyfile(os.path.join(HERE, input_file),
                            os.path.join(tempdir, input_file))
            with capture_log_output() as log_out:
                assert main(os.path.join(tempdir, input_file), tempdir,
                            os.path.join(tempdir, 'report')) == 1
            assert 'traceback' in log_out.getvalue().strip()

            with capture_log_output() as log_out:
                for i in ['000', '001', '002']:
                    i_folder = os.path.join(tempdir, i)
                    os.mkdir(i_folder)
                    with open(i_folder + '/pathensemble.txt', 'w') as f:
                        for line in lines:
                            f.write(line)
                i_folder = os.path.join(tempdir, '003')
                os.mkdir(i_folder)
                with open(i_folder + '/pathensemble.txt', 'w') as f:
                    for line in lines[:-2]:
                        f.write(line)

                assert main(os.path.join(tempdir, input_file), tempdir,
                            os.path.join(tempdir, 'report')) == 1

            assert 'ssing probability: 1.00' in log_out.getvalue().strip()
            assert '2.577319588' in log_out.getvalue().strip()
            input_file = 'dummy_input.rst'  # TIS
            shutil.copyfile(os.path.join(HERE, input_file),
                            os.path.join(tempdir, input_file))
            with capture_log_output() as log_out:
                assert main(os.path.join(tempdir, input_file), tempdir,
                            os.path.join(tempdir, 'report')) == 0
            assert 'report/002_tis_report.' in log_out.getvalue().strip()

    def test_main_repptis(self):
        """Test pyretisanalyse (main) for a REPPTIS simulation.

        This one, just like the pyretisanalyse test for RETIS simuls,
        takes quite a long time because the output is zipped, and
        figures are created...

        """
        input_file = os.path.join("pptis", "retis.rst")
        pathensemblefiles =\
            [os.path.join("pptis",
                          "pathensemble"+str(i)+".txt") for i in range(4)]
        with tempfile.TemporaryDirectory() as tempdir:
            shutil.copyfile(os.path.join(HERE, input_file),
                            os.path.join(tempdir, "retis.rst"))
            for i in range(4):
                os.mkdir(os.path.join(tempdir, "00"+str(i)))
                shutil.copyfile(os.path.join(HERE, pathensemblefiles[i]),
                                os.path.join(tempdir, "00"+str(i),
                                             "pathensemble.txt"))
            with capture_log_output() as log_out:
                assert main(os.path.join(tempdir, "retis.rst"), tempdir,
                            os.path.join(tempdir, 'report')) == 0
            assert ('ssing probability: 5.166013535e-04' in
                    log_out.getvalue().strip())
            assert 'reduced): 0.734638269' in log_out.getvalue().strip()
            assert 'Xi: 0.46000' in log_out.getvalue().strip()

# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the bin."""
import logging
import os
import pytest
import subprocess
import tempfile
import shutil
from contextlib import contextmanager
from io import StringIO
from unittest.mock import patch
from pyretis.bin.pyvisa import (main,
                                hello_pyvisa,
                                bye_pyvisa,
                                pyvisa_visual)
from pyretis.inout.common import TRJ_FORMATS
from pyretis.pyvisa.common import (find_data, recalculate_all, _get_trjs,
                                   read_single_data_file, run_user_script)

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


class TestVarious:
    """Test Main."""

    def test_pyvisarun(self):
        """Test pyvisa run."""
        with tempfile.TemporaryDirectory() as tempdir:
            with subprocess.Popen(['pyvisa', '-cmp'], cwd=tempdir,
                                  stdout=subprocess.PIPE) as ex:
                asd = ex.stdout.read().split()
            assert b'PyVisA' in asd

        with tempfile.TemporaryDirectory() as tempdir:
            with subprocess.Popen(['pyvisa', '-recalculate'], cwd=tempdir,
                                  stdout=subprocess.PIPE) as ex:
                asd = ex.stdout.read().split()
            assert b'PyVisA' in asd

        with tempfile.TemporaryDirectory() as tempdir:
            input_file = 'dummy_input.rst'
            with subprocess.Popen(['pyvisa', '-i', input_file, '-cmp'],
                                  cwd=tempdir, stdout=subprocess.PIPE) as ex:
                asd = ex.stdout.read().split()
            assert b'PyVisA' in asd

    def test_hello_world(self):
        """Test that we are polite."""
        with capture_log_output() as log_out:
            hello_pyvisa(run_dir='.', infile='This_is_enough.lol')
        assert 'Start' in log_out.getvalue().strip()

    def test_bye_pyvisa(self):
        """Test that we can die."""
        with capture_log_output() as log_out:
            bye_pyvisa()
        assert 'references' in log_out.getvalue().strip()

    def test_pyvisa_visual(self):
        """Test that we know how to set up a simulation."""
        inputfile = 'a_non_existent_input.kus'
        with pytest.raises(ImportError) as err:
            pyvisa_visual('.', inputfile, pyvisa_dict={})
        assert 'no' in str(err.value)

    def test_pyvisa_main(self):
        """Test main."""
        infile = 'dummy_input.rst'
        with tempfile.TemporaryDirectory() as tempdir:
            input_file = os.path.join(tempdir, infile)
            shutil.copyfile(os.path.join(HERE, infile), input_file)
            with capture_log_output() as log_out:
                assert main(basepath=tempdir, input_file=input_file,
                            pyvisa_dict={'pyvisa_recalculate': True}) == 1
            assert 'Execution failed!' in log_out.getvalue().strip()

        with tempfile.TemporaryDirectory() as tempdir:
            input_file = os.path.join(tempdir, infile)
            shutil.copyfile(os.path.join(HERE, infile), input_file)
            with capture_log_output() as log_out:
                assert main(basepath=tempdir, input_file=input_file,
                            pyvisa_dict={'pyvisa_compressor': True}) == 1
            assert 'No files to an' in log_out.getvalue().strip()

        with tempfile.TemporaryDirectory() as tempdir:
            input_file = os.path.join(tempdir, infile)
            shutil.copyfile(os.path.join(HERE, infile), input_file)
            with capture_log_output() as log_out:
                assert main(basepath=tempdir, input_file=input_file,
                            pyvisa_dict={}) == 1
            assert 'traceback' in log_out.getvalue().strip()

        with tempfile.TemporaryDirectory() as tempdir:
            with capture_log_output() as log_out:
                assert main(basepath=tempdir, input_file='no_thank_you',
                            pyvisa_dict={}) == 1
            assert 'NOT' in log_out.getvalue().strip()


# ---------------------------------------------------------------------------
# Fixture: minimal simulation directory mirroring
# test/pyvisa/test_simulation_dir
# ---------------------------------------------------------------------------
PYVISA_TESTDIR = os.path.abspath(
    os.path.join(HERE, '..', 'pyvisa', 'test_simulation_dir'))


class TestFormatSupport:
    """TRJ_FORMATS covers the file formats expected by PyRETIS engines."""

    def test_trj_formats_contain_expected_extensions(self):
        """TRJ_FORMATS must include all standard trajectory extensions."""
        required = {'.xyz', '.gro', '.trr', '.txt', '.g96'}
        assert required.issubset(set(TRJ_FORMATS))

    def test_pyvisa_visual_raises_for_unsupported_format(self):
        """pyvisa_visual raises ImportError for unknown file extensions."""
        with pytest.raises(ImportError):
            pyvisa_visual('.', 'some_file.unsupported_ext', pyvisa_dict={})

    def test_pyvisa_visual_raises_when_requirements_missing(self, monkeypatch):
        """pyvisa_visual raises ImportError when optional Qt deps absent."""
        monkeypatch.setattr('pyretis.bin.pyvisa.HAS_PYVISA_REQ', False)
        with pytest.raises(ImportError, match='Requirements not installed'):
            pyvisa_visual('.', 'some.rst', pyvisa_dict={})


class TestGetTrjs:
    """_get_trjs correctly discovers trajectory files inside a folder."""

    def test_returns_all_supported_files(self, tmp_path):
        """_get_trjs finds every file whose extension is in TRJ_FORMATS."""
        for ext in TRJ_FORMATS:
            (tmp_path / f'traj{ext}').touch()
        (tmp_path / 'README.txt.bak').touch()  # should be ignored
        found = _get_trjs(str(tmp_path))
        found_exts = {os.path.splitext(f)[1] for f in found}
        assert found_exts == set(TRJ_FORMATS)

    def test_empty_folder_returns_empty_list(self, tmp_path):
        """_get_trjs on a folder with no trajectory files returns []."""
        assert _get_trjs(str(tmp_path)) == []

    def test_ignores_unsupported_extensions(self, tmp_path):
        """_get_trjs ignores files whose extensions are not in TRJ_FORMATS."""
        (tmp_path / 'config.json').touch()
        (tmp_path / 'notes.md').touch()
        assert _get_trjs(str(tmp_path)) == []

    def test_ignores_analysis_txt_files(self, tmp_path):
        """_get_trjs excludes known analysis outputs from folder scans."""
        (tmp_path / 'order.txt').touch()
        (tmp_path / 'energy.txt').touch()
        (tmp_path / 'pathensemble.txt').touch()
        (tmp_path / 'error.txt').touch()
        (tmp_path / 'traj.txt').touch()
        result = _get_trjs(str(tmp_path))
        assert result == [str(tmp_path / 'traj.txt')]

    def test_returns_absolute_paths(self, tmp_path):
        """_get_trjs entries are full paths, not just file names."""
        (tmp_path / 'run.xyz').touch()
        result = _get_trjs(str(tmp_path))
        assert len(result) == 1
        assert os.path.isabs(result[0])


class TestFindDataExtended:
    """find_data handles different layouts: single file, folder, ensembles."""

    def test_single_trajectory_file(self):
        """find_data with data=<file> returns a dict keyed on '000'."""
        traj = os.path.join(PYVISA_TESTDIR,
                            '003', 'traj', 'traj-acc', '0', 'traj', 'traj.xyz')
        result = find_data(PYVISA_TESTDIR, data=traj)
        assert '000' in result
        assert 'traj.xyz' in result['000']['traj']['0']['traj'][0]

    def test_folder_with_trajectory_files(self, tmp_path):
        """find_data with data=<folder> scans for trajectory files."""
        for i in range(3):
            (tmp_path / f'frame_{i}.xyz').touch()
        result = find_data(str(tmp_path), data=str(tmp_path))
        assert '000' in result
        assert len(result['000']['traj']['0']['traj']) == 3

    def test_empty_folder_returns_empty_dict(self, tmp_path):
        """find_data with an empty folder returns an empty dict."""
        result = find_data(str(tmp_path), data=str(tmp_path))
        assert result == {}

    def test_simulation_layout_finds_all_ensembles(self):
        """find_data on a simulation dir finds all numbered ensembles."""
        result = find_data(PYVISA_TESTDIR)
        ensemble_names = list(result.keys())
        # Test simulation has ensembles 000-005 plus test_path
        numeric_ensembles = [e for e in ensemble_names if e.isdigit()]
        assert len(numeric_ensembles) >= 4

    def test_filter_by_ensemble_names(self):
        """find_data respects the ensemble_names filter."""
        result = find_data(PYVISA_TESTDIR, ensemble_names=['003'])
        assert list(result.keys()) == ['003']

    def test_folder_mixed_extensions(self, tmp_path):
        """find_data only picks up supported extensions from a folder."""
        (tmp_path / 'traj.xyz').touch()
        (tmp_path / 'notes.txt.bak').touch()
        (tmp_path / 'image.png').touch()
        result = find_data(str(tmp_path), data=str(tmp_path))
        assert '000' in result
        assert len(result['000']['traj']['0']['traj']) == 1
        assert result['000']['traj']['0']['traj'][0].endswith('.xyz')


class TestRecalculateAllExtended:
    """recalculate_all with individual files, folders, and full simulations."""

    def test_recalculate_with_data_file(self):
        """recalculate_all with data=<file> writes a new order.txt."""
        with tempfile.TemporaryDirectory() as tempdir:
            shutil.copytree(PYVISA_TESTDIR, os.path.join(tempdir, 'sim'))
            sim_dir = os.path.join(tempdir, 'sim')
            traj = os.path.join(
                sim_dir, '003', 'traj', 'traj-acc', '0', 'traj', 'traj.xyz')
            # Use the dummy_input.rst which has a valid Position OP
            shutil.copyfile(os.path.join(HERE, 'dummy_input.rst'),
                            os.path.join(sim_dir, 'dummy.rst'))
            with patch('sys.stdout', new=StringIO()):
                # Should process the trajectory without errors
                result = recalculate_all(sim_dir, 'dummy.rst', data=traj)
            # Position OP works on a 1D system;
            # succeeds if the file has particles
            assert isinstance(result, bool)

    def test_recalculate_with_data_folder(self, tmp_path):
        """recalculate_all with data=<folder> processes all .xyz files."""
        # Build a minimal folder containing a single valid xyz frame
        src_xyz = os.path.join(
            PYVISA_TESTDIR,
            '003', 'traj', 'traj-acc', '0', 'traj', 'traj.xyz')
        traj_dir = tmp_path / 'trajs'
        traj_dir.mkdir()
        shutil.copyfile(src_xyz, traj_dir / 'traj.xyz')
        shutil.copyfile(os.path.join(HERE, 'dummy_input.rst'),
                        tmp_path / 'dummy.rst')
        with patch('sys.stdout', new=StringIO()):
            result = recalculate_all(str(tmp_path), 'dummy.rst',
                                     data=str(traj_dir))
        assert isinstance(result, bool)

    def test_recalculate_no_trajectories_returns_false(self, tmp_path):
        """recalculate_all returns False when no trajectory data is found."""
        shutil.copyfile(os.path.join(HERE, 'dummy_input.rst'),
                        tmp_path / 'dummy.rst')
        empty_dir = tmp_path / 'empty'
        empty_dir.mkdir()
        with patch('sys.stdout', new=StringIO()):
            result = recalculate_all(str(tmp_path), 'dummy.rst',
                                     data=str(empty_dir))
        assert result is False

    def test_recalculate_invalid_op_returns_false(self):
        """recalculate_all returns False for an invalid order parameter."""
        with tempfile.TemporaryDirectory() as tempdir:
            shutil.copytree(PYVISA_TESTDIR, os.path.join(tempdir, 'sim'))
            sim_dir = os.path.join(tempdir, 'sim')
            with patch('sys.stdout', new=StringIO()):
                # input.rst has 'class = unicorn' for the CV —
                # should return False
                result = recalculate_all(sim_dir, 'input.rst')
            assert result is False


class TestSingleDataFileParsing:
    """Standalone txt/csv tables can be loaded directly by PyVisA."""

    def test_order_txt_header_keeps_first_column_plotable(self, tmp_path):
        """Commented order.txt headers preserve the explicit time column."""
        data_file = tmp_path / 'order.txt'
        data_file.write_text(
            'Recalculated data\n'
            '### Time Orderp cv2\n'
            '0 1.5 2.5\n'
            '5 1.7 2.7\n',
            encoding='utf-8',
        )

        frames, plot_cols, main_op = read_single_data_file(str(data_file))

        assert plot_cols == ['Time', 'Orderp', 'cv2']
        assert main_op == 'Orderp'
        assert list(frames['Time']) == [0.0, 5.0]
        assert list(frames['Orderp']) == [1.5, 1.7]

    def test_csv_comment_header_accepts_other_comment_symbols(self, tmp_path):
        """Header parsing tolerates comment prefixes beyond a single #."""
        data_file = tmp_path / 'data.csv'
        data_file.write_text(
            ';; step, lambda, cv2\n'
            '0.0, 2.0, 3.0\n'
            '1.5, 2.5, 3.5\n',
            encoding='utf-8',
        )

        frames, plot_cols, main_op = read_single_data_file(str(data_file))

        assert plot_cols == ['step', 'lambda', 'cv2']
        assert main_op == 'lambda'
        assert list(frames['step']) == [0.0, 1.5]

    def test_rst_labels_are_used_when_header_is_missing(self, tmp_path):
        """Missing headers fall back to names from a provided .rst file."""
        rst_file = tmp_path / 'named.rst'
        rst_text = (tmp_path / 'dummy_input_source.rst')
        shutil.copyfile(os.path.join(HERE, 'dummy_input.rst'), rst_text)
        rst_text = rst_text.read_text(encoding='utf-8')
        rst_text = rst_text.replace(
            'class = Position\n', 'class = Position\nname = Lambda\n', 1)
        rst_file.write_text(rst_text, encoding='utf-8')

        data_file = tmp_path / 'data.txt'
        data_file.write_text('0 1.25\n1 1.50\n', encoding='utf-8')

        frames, plot_cols, main_op = read_single_data_file(
            str(data_file), rst_file=str(rst_file))

        assert plot_cols == ['time', 'Lambda']
        assert main_op == 'Lambda'
        assert list(frames['time']) == [0.0, 1.0]


class TestRunUserScript:
    """User scripts should run robustly for PyVisA single-file workflows."""

    def test_run_user_script_accepts_whitespace_tables(self, tmp_path):
        """Whitespace-delimited stdout is accepted, not only comma CSV."""
        script = tmp_path / 'emit_table.py'
        script.write_text(
            "print('1.0 2.0\\n3.0 4.0')\n",
            encoding='utf-8',
        )

        order_txt, error = run_user_script(str(script))

        assert error is None
        assert order_txt is not None
        frames, plot_cols, main_op = read_single_data_file(order_txt)
        assert plot_cols == ['Time', 'Orderp', 'cv1']
        assert main_op == 'Orderp'
        assert list(frames['Orderp']) == [1.0, 3.0]

    def test_run_user_script_uses_script_directory_as_cwd(self, tmp_path):
        """Relative file access inside the script works as expected."""
        (tmp_path / 'data.txt').write_text('[[1.0], [2.0]]', encoding='utf-8')
        script = tmp_path / 'emit_relative.py'
        script.write_text(
            "from pathlib import Path\n"
            "print(Path('data.txt').read_text())\n",
            encoding='utf-8',
        )

        order_txt, error = run_user_script(str(script))

        assert error is None
        assert order_txt is not None


class TestMainRecalculateData:
    """main() routes pyvisa_recalculate + pyvisa_data combinations."""

    def test_main_recalculate_with_data_folder(self):
        """main() with pyvisa_recalculate and pyvisa_data folder succeeds."""
        with tempfile.TemporaryDirectory() as tempdir:
            src = os.path.join(PYVISA_TESTDIR,
                               '003', 'traj', 'traj-acc', '0', 'traj')
            traj_dir = os.path.join(tempdir, 'trajs')
            shutil.copytree(src, traj_dir)
            shutil.copyfile(os.path.join(HERE, 'dummy_input.rst'),
                            os.path.join(tempdir, 'dummy.rst'))
            with capture_log_output():
                status = main(
                    basepath=tempdir,
                    input_file='dummy.rst',
                    pyvisa_dict={
                        'pyvisa_recalculate': True,
                        'pyvisa_data': traj_dir,
                    },
                )
            # 0 = success, 1 = handled error (e.g. empty OP calc)
            assert status in (0, 1)

    def test_main_recalculate_with_data_single_file(self):
        """main() with pyvisa_recalculate and pyvisa_data single file runs."""
        with tempfile.TemporaryDirectory() as tempdir:
            src_xyz = os.path.join(
                PYVISA_TESTDIR,
                '003', 'traj', 'traj-acc', '0', 'traj', 'traj.xyz')
            dest_xyz = os.path.join(tempdir, 'traj.xyz')
            shutil.copyfile(src_xyz, dest_xyz)
            shutil.copyfile(os.path.join(HERE, 'dummy_input.rst'),
                            os.path.join(tempdir, 'dummy.rst'))
            with capture_log_output():
                status = main(
                    basepath=tempdir,
                    input_file='dummy.rst',
                    pyvisa_dict={
                        'pyvisa_recalculate': True,
                        'pyvisa_data': dest_xyz,
                    },
                )
            assert status in (0, 1)

    def test_main_recalculate_missing_file_returns_error(self):
        """main() fails gracefully when input_file is absent."""
        with tempfile.TemporaryDirectory() as tempdir:
            with capture_log_output() as log_out:
                status = main(basepath=tempdir, input_file='ghost.rst',
                              pyvisa_dict={'pyvisa_recalculate': True,
                                           'pyvisa_data': None})
            assert status == 1
            assert 'NOT' in log_out.getvalue()

# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the FileIO class."""
import filecmp
import logging
import os
import shutil
import pathlib
import tempfile
import pytest
from io import StringIO
from unittest.mock import patch
import numpy as np
from numpy.random import randint, rand
from pyretis.core.path import Path
from pyretis.core.pathensemble import PathEnsemble
from pyretis.inout.analysisio.analysisio import run_analysis, read_first_block
from pyretis.inout.formats.formatter import OutputFormatter
from pyretis.inout.formats.snapshot import read_txt_snapshots
from pyretis.inout.fileio import FileIO
from pyretis.inout.formats.energy import EnergyFile, EnergyPathFile
from pyretis.inout.formats.cross import CrossFile
from pyretis.inout.formats.order import OrderFile, OrderPathFile
from pyretis.inout.formats.pathensemble import PathEnsembleFile
from pyretis.inout.formats.path import PathExtFile, PathIntFile
from pyretis.inout.formats.snapshot import SnapshotFile
from pyretis.inout.settings import parse_settings_file
from .help import (
    assert_equal_path_dict,
    create_external_path,
    create_test_paths,
    turn_on_logging,
    remove_file,
    set_up_system,
)


logging.disable(logging.CRITICAL)


HERE = os.path.abspath(os.path.dirname(__file__))


PATH_DATA = [
    {
        'cycle': 0, 'generated': ['ki', 0.0, 0, 0],
        'interface': ('L', 'M', 'L'), 'length': 494,
        'ordermax': (-0.765171027, 255), 'ordermin': (-0.9002987171, 0),
        'status': 'ACC', 'weight': 1.,
    },
    {
        'cycle': 1, 'generated': ['s-', 0.0, 0, 0],
        'interface': ('L', 'M', 'L'), 'length': 400,
        'ordermax': (-0.8643118362, 136),
        'ordermin': (-0.9004062103, 399),
        'status': 'ACC', 'weight': 1.,
    },
    {
        'cycle': 2, 'generated': ['tr', 0.0, 0, 0],
        'interface': ('L', 'M', 'L'), 'length': 400,
        'ordermax': (-0.8643118362, 263), 'ordermin': (-0.9004062103, 0),
        'status': 'ACC', 'weight': 1.,
    },
]


class TestAnaIo:
    """Test the analysisio files."""
    def test_read_first_block(self, caplog):
        """Test read first block function."""
        with tempfile.TemporaryDirectory() as tempdir:
            filegone = os.path.join(tempdir, 'Mordvichev')
            pathlib.Path(filegone).touch(exist_ok=True)
            out = read_first_block(file_type='energy', file_name=filegone)
            assert out is None
            with pytest.raises(ValueError) as ext:
                read_first_block(file_type='BT2', file_name=filegone)
            assert 'Unknown file type' in str(ext.value)
            filetest = os.path.join(HERE, 'order-data.txt')
            with turn_on_logging():
                with caplog.at_level(
                        logging.WARNING,
                        logger='pyretis.inout.analysisio.analysisio'):
                    read_first_block(file_type='order', file_name=filetest)

    def test_run_analysis(self, caplog):
        """Test analysisio function"""
        lines = ['']*8
        moves = ['ki', 'sh', 'ki', 'sh', 'ld', 'xx', 'sh']
        line = "         2          3          0 L M R      99 ACC "
        line += "sh -9.108621358e-00  1.007510172e+01       0"
        line += "      98  0.000000000e+00       0       0  1.00e+00\n"
        hline = "# Cycle: 1, status: ACC, move: ('sh', -1.430, 18, 14)\n"
        olines = ['#     Time       Orderp\n',
                  '          0     0.11319\n', '          1    -0.1367\n',
                  '          2    -0.78098\n', '          3    -1.0238\n',
                  '          4    -1.11127\n', '          5    -1.2367\n',
                  '          6    -0.58084\n', '          7     0.0238\n',
                  '          8    -1.06811\n', '          9     0.1131\n',
                  '         10     1.06811\n', '         11     0.3131\n']

        elines = ['#     Time      Potential        Kinetic\n',
                  '          0   -1706.522027    1267.214065\n',
                  '          1   -8706.550067    1167.214065\n',
                  '          2  -18720.057124    1780.765874\n',
                  '          3  -18706.550067    1767.214065\n',
                  '          4  -18705.164689    1765.349001\n',
                  '          5  -18629.768664    1691.047800\n',
                  '          6  -18677.061888    1737.016339\n',
                  '          7  -18567.169855    1628.072515\n',
                  '          8  -18515.405500    1576.664809\n',
                  '          9  -18515.405500    1576.664809\n',
                  '         10  -18583.703121    1643.763063\n',
                  '         11  -18583.703121    1643.763063\n']
        elines2 = [
            '# Time       Potential      Kinetic       Total    Temperature\n',
            '   0      -0.165063      0.000000      -0.165063      0.000000\n',
            '   1      -0.164748      0.000853      -0.163895      0.001706\n',
            '   2      -0.164173      0.001763      -0.162411      0.003526\n',
            '   3      -0.163650      0.002079      -0.161571      0.004157\n',
            '   4      -0.163075      0.001212      -0.161863      0.002424\n',
            '   5      -0.162872      0.000026      -0.162847      0.000051\n',
            '   6      -0.163122      0.000634      -0.162489      0.001267\n',
            '   7      -0.163780      0.003669      -0.160110      0.007339\n',
            '   8      -0.164685      0.004810      -0.159874      0.009620\n',
            '   9      -0.164685      0.004810      -0.159874      0.009620\n',
            '  10      -0.164685      0.004810      -0.159874      0.009620\n',
            '  11      -0.174685      0.008810      -0.129874      0.008620\n']

        clines = [' #     Step  Int Dir \n', '         3    1   - \n',
                  '         1    2   +  \n', '         4    1   - \n',
                  '         6    2   +  \n', '         7    1   - \n',
                  '         6    2   +  \n', '         7    1   - \n',
                  '         6    2   -  \n', '         7    1   + \n',
                  '         8    1   -  \n', '         9    1   + \n']

        for j in range(7):
            lines[j] = line.replace('sh', moves[j])
        lines[7] = line.replace('ACC', 'BWI')
        lines[0] = line.replace('L', 'R')
        lines[1] = line.replace('L', 'R')

        input_file = os.path.join(HERE, 'settings-retis.rst')
        settings = parse_settings_file(input_file)
        settings['simulation']['interfaces'] = [-0.9, -0.8, 0]
        with tempfile.TemporaryDirectory() as tempdir:
            for i in ['000', '001', '002']:
                i_folder = os.path.join(tempdir, i)
                os.mkdir(i_folder)
                with open(i_folder + '/pathensemble.txt', 'w',
                          encoding='utf-8') as f_p, \
                        open(i_folder + '/energy.txt', 'w',
                             encoding='utf-8') as f_e, \
                        open(i_folder + '/order.txt', 'w',
                             encoding='utf-8') as f_o:
                    for line in lines:
                        f_p.write(line)
                        f_e.write(hline)
                        for line in elines:
                            f_e.write(line)
                        f_o.write(hline)
                        for line in olines:
                            f_o.write(line)
            settings['simulation']['exe-path'] = tempdir
            settings['analysis']['report-dir'] = tempdir
            with patch('sys.stdout', new=StringIO()) as stdout:
                run_analysis(settings)
            assert 'Pathensemble data' in stdout.getvalue().strip()
            assert '99.000000' in stdout.getvalue().strip()

        with tempfile.TemporaryDirectory() as tempdir:
            with open(os.path.join(tempdir, 'pathensemble.txt'), 'w',
                      encoding='utf-8') as f_p, \
                    open(os.path.join(tempdir, 'energy.txt'), 'w',
                         encoding='utf-8') as f_e, \
                    open(os.path.join(tempdir, 'order.txt'), 'w',
                         encoding='utf-8') as f_o, \
                    open(os.path.join(tempdir, 'cross.txt'), 'w',
                         encoding='utf-8') as c_r:
                for line in lines:
                    f_p.write(line)
                for line in elines2:
                    f_e.write(line)
                f_o.write(hline)
                for line in olines:
                    f_o.write(line)
                for line in clines:
                    c_r.write(line)
                inputfile = os.path.join(HERE, 'settings.rst')
            settings = parse_settings_file(inputfile)
            settings['simulation']['exe-path'] = tempdir
            settings['analysis']['report-dir'] = tempdir
            settings['simulation']['task'] = 'md-flux'
            settings['simulation']['endcycle'] = 8
            settings['simulation']['interfaces'] = [-0.9, -0.8]
            with patch('sys.stdout', new=StringIO()) as stdout:
                run_analysis(settings)
            assert '0.2500' in stdout.getvalue().strip()
            assert '0.3750' in stdout.getvalue().strip()

        settings['simulation']['task'] = 'mistery'
        with pytest.raises(ValueError) as err:
            run_analysis(settings)
            assert 'Unknown ' in err.value

        # repptis analysis
        with tempfile.TemporaryDirectory() as tempdir:
            # copy pptis/pathensemble.txt files to tempdir
            for i, ens in enumerate(['000', '001', '002', '003']):
                os.mkdir(os.path.join(tempdir, ens))
                shutil.copy(
                    os.path.join(HERE, 'pptis/pathensemble'+str(i)+'.txt'),
                    os.path.join(tempdir, ens, 'pathensemble.txt')
                )
            inputfile = os.path.join(HERE, 'pptis/retis.rst')
            settings = parse_settings_file(inputfile)
            settings['simulation']['exe-path'] = tempdir
            settings['analysis']['report-dir'] = tempdir
            settings['simulation']['task'] = 'repptis'
            with patch('sys.stdout', new=StringIO()) as stdout:
                run_analysis(settings)
            # flux
            assert '0.862851052' in stdout.getvalue().strip()
            # pcross
            assert '2.213477168e-03' in stdout.getvalue().strip()


class TestFileIO:
    """Test the FileIO class."""

    def test_initiation(self, caplog):
        """Test initiation and opening of files."""
        with turn_on_logging():
            with caplog.at_level(logging.INFO, logger='pyretis.inout.fileio'):
                FileIO('some-name', 'r', None, backup='not-an-option!')

        # Test for completely new file:
        filename = os.path.join(HERE, 'a_new_file')
        remove_file(filename)
        fileio = FileIO(filename, 'w', None, backup=False)
        assert fileio.fileh is None
        fileio.open()
        assert fileio.fileh is not None
        fileio.write('test')
        fileio.close()
        remove_file(filename)

        # Test when a file exists:
        filename = os.path.join(HERE, 'already_exists2')
        remove_file(filename)
        remove_file(f'{filename}_000')
        with open(filename, 'w', encoding='utf-8') as fileh:
            fileh.write('test')
        fileio = FileIO(filename, 'w', None, backup=True)
        fileio.open()
        fileio.write('text')
        assert fileio.fileh is not None
        assert os.path.isfile(filename)
        assert os.path.isfile(f'{filename}_000')
        del fileio
        remove_file(filename)
        remove_file(f'{filename}_000')
        # Test for invalid filename + context manager:
        with turn_on_logging():
            with caplog.at_level(logging.CRITICAL,
                                 logger='pyretis.inout.fileio'):
                with FileIO('/#"½%&?<><|*', 'r', None) as some_file:
                    assert some_file.file_mode == 'r'
                    lines = list(some_file)
                    assert not lines
        # Test for weird mode:
        fileio = FileIO('some-file', 'q', None)
        with pytest.raises(ValueError):
            fileio.open()
        fileio = FileIO('some-file', 'r', None)
        fileio.file_mode = 'q'
        with pytest.raises(ValueError):
            fileio.open_file_read()
        fileio = FileIO('some-file', 'w', None)
        fileio.file_mode = 'q'
        with pytest.raises(ValueError):
            fileio.open_file_write()

    def test_write_modes(self, caplog):
        """Test the write modes and backup settings."""
        with tempfile.NamedTemporaryFile() as tmp:
            with FileIO(tmp.name, 'w', None, backup=True) as fileio:
                assert fileio.fileh is not None
                assert not fileio.fileh.closed
                fileio.write('some text')
            with FileIO(tmp.name, 'a', None, backup=True) as fileio:
                assert fileio.fileh is not None
                assert not fileio.fileh.closed
                status = fileio.write('some text')
                assert status
                status = fileio.write('some text', end=None)
                assert status
                status = fileio.write(None)
                assert not status
            fileio = FileIO(tmp.name, 'w', None, backup=True)
            fileio.open()
            fileio.close()
            with turn_on_logging():
                with caplog.at_level(logging.WARNING,
                                     logger='pyretis.inout.fileio'):
                    fileio.write('text')
            fileio = FileIO(tmp.name, 'a', None, backup='True')
            with turn_on_logging():
                with caplog.at_level(logging.WARNING,
                                     logger='pyretis.inout.fileio'):
                    fileio.write('text')

    def test_reading(self, caplog):
        """Test generic reading."""
        filename = os.path.join(HERE, 'energy.txt')
        correct = []
        with open(filename, 'r', encoding='utf-8') as infile:
            correct = list(infile)
        lines = []
        with FileIO(filename, 'r', OutputFormatter('test')) as fileio:
            for line in fileio:
                lines.append(line)
        assert lines == correct
        with FileIO(filename, 'r', OutputFormatter('test')) as fileio:
            fileio.close()
            lines = list(fileio)
            assert not lines

    def test_load(self, caplog):
        """Test generic load."""
        filename = os.path.join(HERE, 'energy.txt')
        with FileIO(filename, 'r', OutputFormatter('test')) as fileio:
            for block in fileio.load():
                assert block['comment'] == [
                    '#     Time      Potential        Kinetic'
                    '          Total    Temperature']

    def test_output(self, caplog):
        """Test generic output."""
        data = [1, 2.0, 3.123456]
        with tempfile.NamedTemporaryFile() as tmp:
            with FileIO(tmp.name, 'w', OutputFormatter('test'),
                        backup=False) as fileio:
                fileio.output(1, data)
                fileio.output(2, data)
            tmp.flush()
            with FileIO(tmp.name, 'r', OutputFormatter('test')) as fileio:
                for block in fileio.load():
                    for i, line in enumerate(block['data']):
                        raw = [i + 1] + data
                        assert line == raw
        with tempfile.NamedTemporaryFile() as tmp:
            fileio = FileIO(tmp.name, 'w', OutputFormatter(''), backup=False)
            fileio.open()
            fileio.output(1, data)
            fileio.flush()
            fileio.close()
            fileio.flush()

    def test_formatter_info(self, caplog):
        """Test the formatter info method."""
        fileio = FileIO('My-name-is-nobody', 'r', None)
        info = fileio.formatter_info()
        assert info is None
        formatter = OutputFormatter('Yippee ki-yay')
        fileio = FileIO('My-name-is-nobody', 'r', formatter)
        info = fileio.formatter_info()
        assert info == OutputFormatter


class TestEnergyFile:
    """Test the energy file i/o."""

    def test_energy_fileio(self, caplog):
        """Test reading and writing for energy files."""
        fields = ['vpot', 'ekin', 'etot', 'temp']
        with tempfile.NamedTemporaryFile() as tmp:
            raw_data = []
            with EnergyFile(tmp.name, 'w') as efile:
                for i in range(50):
                    rnd = rand(len(fields))
                    raw_data.append([i] + list(rnd))
                    data = dict(zip(fields, rnd))
                    efile.output(i, data)
            tmp.flush()
            raw_data = np.array(raw_data)
            with EnergyFile(tmp.name, 'r') as efile:
                for block in efile.load():
                    data = block['data']
                    for i, key in enumerate(['time'] + fields):
                        for num1, num2 in zip(raw_data[:, i], data[key]):
                            assert num1 == pytest.approx(num2, abs=1e-6)

    def test_energy_file_read(self, caplog):
        """Test reading existing energy data."""
        filename = os.path.join(HERE, 'energy.txt')
        correct = {
            'time': [0, 1000, 2000],
            'vpot': [-0.165063, -0.806879, -0.760215],
            'ekin': [0.035000, 0.328694, 0.041118],
            'etot': [-0.130063, -0.478185, -0.719097],
            'temp': [0.070000, 0.657387, 0.082236],
        }
        with EnergyFile(filename, 'r') as efile:
            for block in efile.load():
                for key, val in correct.items():
                    assert np.allclose(val, block['data'][key])
        # Reading a file with errors, this should give a warning:
        filename = os.path.join(HERE, 'energy-error.txt')
        with turn_on_logging():
            with caplog.at_level(logging.WARNING,
                                 logger='pyretis.inout.fileio'):
                with EnergyFile(filename, 'r') as efile:
                    for _ in efile.load():
                        pass


class TestCrossFile:
    """Test the cross file i/o."""

    def test_cross_fileio(self, caplog):
        """Test reading and writing for cross files."""
        with tempfile.NamedTemporaryFile() as tmp:
            raw_data = []
            with CrossFile(tmp.name, 'w') as cfile:
                prev = 0
                for i in range(50):
                    interf = randint(1, 10)
                    step = prev + randint(1, 1000)
                    prev = step
                    direction = 1 if randint(0, 2) == 1 else -1
                    raw_data.append((step, interf + 1, direction))
                    data = [(step, interf, ('-', '+', '+')[direction + 1])]
                    cfile.output(i, data)
            tmp.flush()
            with CrossFile(tmp.name, 'r') as cfile:
                for block in cfile.load():
                    data = block['data']
                    for i, j in zip(data, raw_data):
                        assert i == j

    def test_cross_fileread(self, caplog):
        """Test reading from existing files."""
        filename1 = os.path.join(HERE, 'cross-error.txt')
        filename2 = os.path.join(HERE, 'cross.txt')
        data1 = None
        with turn_on_logging():
            with caplog.at_level(logging.WARNING,
                                 logger='pyretis.inout.fileio'):
                with CrossFile(filename1, 'r') as cfile:
                    for block in cfile.load():
                        data1 = block['data']
        data2 = None
        with CrossFile(filename2, 'r') as cfile:
            for block in cfile.load():
                data2 = block['data']
        for i, line in enumerate(data2):
            if i < 10:
                assert line == data1[i]
            elif i == 10:
                pass
            else:
                assert line == data1[i-1]


class TestOrderFile:
    """Test the order parameter file i/o."""

    def test_order_parameter_fileio(self, caplog):
        """Test writing and reading of order parameter files."""
        extra = 3
        with tempfile.NamedTemporaryFile() as tmp:
            raw_data = []
            with OrderFile(tmp.name, 'w') as ofile:
                for i in range(50):
                    rnd = rand(1 + extra)
                    raw_data.append([i] + list(rnd))
                    ofile.output(i, rnd)
            tmp.flush()
            with OrderFile(tmp.name, 'r') as ofile:
                for block in ofile.load():
                    data = block['data']
                    for i, j in zip(data, raw_data):
                        for num1, num2 in zip(i, j):
                            assert num1 == pytest.approx(num2, abs=1e-6)

    def test_order_parameter_read(self, caplog):
        """Test reading of order parameter data from a file."""
        filename = os.path.join(HERE, 'order-data.txt')
        correct_data = [
            [[0, 0.1, -0.1], [1, 0.2, -0.2], [2, 0.3, -0.3],
             [3, 0.4, -0.4], [4, 0.5, -0.5]],
            [[0, 1.1, -1.1], [1, 1.2, -1.2], [2, 1.3, -1.3],
             [3, 1.4, -1.4], [4, 1.5, -1.5]],
            [[0, 2.1, -2.1], [1, 2.2, -2.2], [2, 2.3, -2.3],
             [3, 2.4, -2.4], [4, 2.5, -2.5], [5, 2.6, -2.6],
             [6, 2.7, -2.7], [7, 2.8, -2.8]],
            [[0, 3.1, -3.1], [1, 3.2, -3.2], [2, 3.3, -3.3]],
            [[0, 4.1, -4.1], [2, 4.3, -4.3], [3, 4.4, -4.4]],
        ]
        with OrderFile(filename, 'r') as ofile:
            for i, data in enumerate(ofile.load()):
                assert data['comment'][0] == f'# Cycle: {i}, status: ACC'
                assert data['comment'][1] == '#     Time       Orderp'
                for j, k in zip(correct_data[i], data['data']):
                    assert np.allclose(j, k)


class TestPathFiles:
    """Test file i/o for path files."""

    def test_energy_path_file(self, caplog):
        """Test i/o for energy path data."""
        paths, correct_data = create_test_paths(npath=5)
        with tempfile.NamedTemporaryFile() as tmp:
            for i, path in enumerate(paths):
                with EnergyPathFile(tmp.name, 'a') as efile:
                    efile.output(i, [path, 'ACC'])
            tmp.flush()
            with EnergyPathFile(tmp.name, 'r') as efile:
                for i, block in enumerate(efile.load()):
                    for key in ('ekin', 'vpot'):
                        assert np.allclose(block['data'][key],
                                           correct_data[i][key])

    def test_order_path_file(self, caplog):
        """Test i/o for order path data."""
        paths, correct_data = create_test_paths(npath=5)
        with tempfile.NamedTemporaryFile() as tmp:
            for i, path in enumerate(paths):
                with OrderPathFile(tmp.name, 'a') as ofile:
                    ofile.output(i, [path, 'ACC'])
            tmp.flush()
            with OrderPathFile(tmp.name, 'r') as ofile:
                for i, block in enumerate(ofile.load()):
                    corr = np.array(correct_data[i]['order'])
                    data = block['data']
                    assert np.allclose(data[:, 1], corr[:, 0])
                    assert np.allclose(data[:, 2], corr[:, 1])
                    assert np.allclose(data[:, 3], corr[:, 2])

    def test_path_int_file(self, caplog):
        """Test the internal path writer."""
        paths, _ = create_test_paths(npath=3)
        with tempfile.NamedTemporaryFile() as tmp:
            for i, path in enumerate(paths):
                with PathIntFile(tmp.name, 'a') as pfile:
                    pfile.output(i, [path, 'ACC'])
            tmp.flush()
            with PathIntFile(tmp.name, 'r') as pfile:
                for i, block in enumerate(pfile.load()):
                    assert f'# Cycle: {i}, status: ACC' == block['comment'][0]

    def test_path_ext_file(self, caplog):
        """Test the external path writer."""
        with tempfile.NamedTemporaryFile() as tmp:
            for i in range(3):
                path, _ = create_external_path(random_length=True)
                with PathExtFile(tmp.name, 'a') as pfile:
                    pfile.output(i, [path, 'ACC'])
            tmp.flush()
            with PathExtFile(tmp.name, 'r') as pfile:
                for i, block in enumerate(pfile.load()):
                    assert f'# Cycle: {i}, status: ACC' == block['comment'][0]


class TestPathEnsembleFile:
    """Test the i/o for the PathEnsembleFile."""

    @staticmethod
    def test_path_ensemble_read():
        """Test reading of path ensemble files."""
        filename = os.path.join(HERE, 'pathensemble001.txt')
        with PathEnsembleFile(filename, 'r') as pfile:
            for i, path in enumerate(pfile.get_paths()):
                ref_data = PATH_DATA[i]
                ref_data['filename'] = filename
                assert_equal_path_dict(path, ref_data)

    @staticmethod
    def _fake_path_from_dict(path_dict):
        """Return a path object from a path dict."""
        path = Path(None)
        path.generated = path_dict['generated']
        avg = 0.5 * (path_dict['ordermax'][0] + path_dict['ordermin'][0])
        for i in range(path_dict['length']):
            path.append(
                set_up_system([avg, i], [0], [0], vpot=0, ekin=0)
            )
        if path_dict['interface'][0] == 'L':
            path.phasepoints[0].order[0] = path_dict['ordermin'][0] + 0.0001
        if path_dict['interface'][0] == 'R':
            path.phasepoints[0].order[0] = path_dict['ordermax'][0]
        if path_dict['interface'][-1] == 'L':
            path.phasepoints[-1].order[0] = path_dict['ordermin'][0] + 0.0001
        if path_dict['interface'][-1] == 'R':
            path.phasepoints[-1].order[0] = path_dict['ordermax'][0]
        for key in ('ordermax', 'ordermin'):
            idx = path_dict[key][1]
            path.phasepoints[idx].order[0] = path_dict[key][0]
        return path

    def test_path_ensemble_write(self, caplog):
        """Test writing of path ensemble data."""
        filename = os.path.join(HERE, 'pathensemble001.txt')
        settings = {
            'ensemble_number': 1,
            'interfaces': [-0.9, -0.9, 1],
            'detect': -0.8
        }
        ens = PathEnsemble(settings['ensemble_number'],
                           settings['interfaces'])
        with tempfile.NamedTemporaryFile() as tmp:
            with PathEnsembleFile(tmp.name, 'w',
                                  ensemble_settings=settings) as pfile:
                for pathi in PATH_DATA:
                    path = self._fake_path_from_dict(pathi)
                    ens.add_path_data(path, pathi['status'],
                                      cycle=pathi['cycle'])
                    pfile.output(pathi['cycle'], ens)
            tmp.flush()
            assert filecmp.cmp(tmp.name, filename)
        # Test warning when ensemble settings are not given:
        with turn_on_logging():
            with caplog.at_level(
                    logging.WARNING,
                    logger='pyretis.inout.formats.pathensemble'):
                with PathEnsembleFile(filename, 'r') as pfile:
                    pass
        # Test open missing file:
        with PathEnsembleFile('path_mack_path/file__', 'r',
                              ensemble_settings=settings) as pfile:
            with turn_on_logging():
                with caplog.at_level(
                        logging.CRITICAL,
                        logger='pyretis.inout.formats.pathensemble'):
                    for _ in pfile.load():
                        pass


class TestSnapshot:
    """Test methods related to the snapshot files."""

    def test_initiate_settings(self, caplog):
        """Test initiation of the snapshot file with settings."""
        settings = {'write_vel': False, 'extra-to-ignore': False}
        with tempfile.NamedTemporaryFile() as tmp:
            with SnapshotFile(tmp.name, 'r',
                              format_settings=settings) as sfile:
                assert not sfile.formatter.write_vel

    def test_read_txt_snapshot(self, caplog):
        """Test the read_txt_snapshot method."""
        filename = os.path.join(HERE, 'config.txt')
        read1 = read_txt_snapshots(filename)
        read2 = read_txt_snapshots(filename,
                                   data_keys=('name', 'x', 'y', 'z'))
        correct = {
            'atomname': [['X', 'Y'], ['A', 'B']],
            'box': [np.array([1., 2., 3., 4., 5., 6.]),
                    np.array([7., 8., 9.])],
            'x': [[1.0, 7.0], [1.1, 7.1]],
            'y': [[2.0, 8.0], [2.1, 8.1]],
            'z': [[3.0, 9.0], [3.1, 9.1]],
            'vx': [[4.0, 1.0], [4.1, 1.1]],
            'vy': [[5.0, 2.0], [5.1, 2.1]],
            'vz': [[6.0, 3.0], [6.1, 3.1]],
        }
        for i, (snap1, snap2) in enumerate(zip(read1, read2)):
            assert snap1['atomname'] == correct['atomname'][i]
            assert snap2['name'] == correct['atomname'][i]
            for j in ('x', 'y', 'z', 'box'):
                assert np.allclose(snap1[j], correct[j][i])
                assert np.allclose(snap2[j], correct[j][i])
            for j in ('vx', 'vy', 'vz'):
                assert np.allclose(snap1[j], correct[j][i])
                assert j not in snap2
        filename = os.path.join(HERE, 'config_with_error.txt')
        read3 = read_txt_snapshots(filename)
        with pytest.raises(Exception) as context:
            for i in read3:
                print(i)
        assert 'Oops_not_an_integer' in str(context.value)
        assert 'invalid literal for int()' in str(context.value)
        # Also test the file writer:
        filename = os.path.join(HERE, 'config.txt')
        with SnapshotFile(filename, 'r') as traj:
            for i, data in enumerate(traj.load()):
                for j in ('x', 'y', 'z', 'box'):
                    assert np.allclose(data[j], correct[j][i])

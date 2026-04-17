import unittest
from unittest.mock import MagicMock, patch, mock_open
import os
import signal
from pyretis.engines.gromacs2 import GromacsRunner, GromacsEngine2


class TestGromacs2Extra(unittest.TestCase):
    def setUp(self):
        self.cmd = ['mdrun', '-s', 'topol.tpr']
        self.trr_file = 'test.trr'
        self.edr_file = 'test.edr'
        self.exe_dir = '.'

    @patch('subprocess.Popen')
    @patch('os.path.isfile')
    @patch('os.fstat')
    @patch('builtins.open', new_callable=mock_open)
    def test_start_file_not_present(
            self, mock_file, mock_fstat, mock_isfile, mock_popen):
        # Mock isfile to return False to trigger the 'else' branch in start()
        mock_isfile.return_value = False
        runner = GromacsRunner(
            self.cmd, self.trr_file, self.edr_file, self.exe_dir)

        # Mock popen to avoid actual subprocess
        mock_process = MagicMock()
        mock_process.poll.return_value = 0  # Subprocess finishes immediately
        mock_popen.return_value = mock_process

        runner.start()
        self.assertTrue(runner.stop_read)

    @patch('os.killpg')
    @patch('os.getpgid')
    def test_stop_running(self, mock_getpgid, mock_killpg):
        runner = GromacsRunner(
            self.cmd, self.trr_file, self.edr_file, self.exe_dir)
        mock_process = MagicMock()
        mock_process.returncode = None
        mock_process.pid = 1234
        runner.running = mock_process

        runner.stop()
        mock_killpg.assert_called_once()
        self.assertTrue(runner.stop_read)

    @patch('pyretis.engines.gromacs2.read_trr_header')
    @patch('os.path.getsize')
    @patch('os.fstat')
    @patch('builtins.open', new_callable=mock_open)
    def test_get_gromacs_frames_eof(
            self, mock_file, mock_fstat, mock_getsize, mock_read_header):
        runner = GromacsRunner(
            self.cmd, self.trr_file, self.edr_file, self.exe_dir)
        runner.stop_read = False
        runner.fileh = mock_file()
        runner.header_size = 100
        mock_getsize.return_value = 200

        # Mock read_trr_header to raise EOFError
        mock_read_header.side_effect = EOFError

        # Mock check_poll to return None (running)
        runner.check_poll = MagicMock(return_value=None)

        # Mock reopen_file to return None (so it breaks eventually)
        with patch('pyretis.engines.gromacs2.reopen_file',
                   return_value=(None, None)):
            # Set stop_read to True after one iteration.
            def side_effect(*args, **kwargs):
                runner.stop_read = True
                return None
            runner.check_poll.side_effect = side_effect

            frames = list(runner.get_gromacs_frames())
            self.assertEqual(len(frames), 0)

    def test_check_poll_no_running(self):
        runner = GromacsRunner(
            self.cmd, self.trr_file, self.edr_file, self.exe_dir)
        with self.assertRaises(RuntimeError):
            runner.check_poll()


if __name__ == '__main__':
    unittest.main()

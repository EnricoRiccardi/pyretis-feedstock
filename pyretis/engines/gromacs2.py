# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""A GROMACS external MD integrator interface.

This module defines a class for using GROMACS as an external engine.

Important classes defined here
------------------------------

GromacsEngine2
    A class responsible for interfacing GROMACS.
"""
import logging
import os
import shlex
import signal
import subprocess
from time import sleep
from pyretis.core.box import box_matrix_to_list
from pyretis.engines.external import ExternalMDEngine
from pyretis.engines.gromacs import GromacsEngine
from pyretis.inout.formats.gromacs import (
    read_trr_header,
    read_trr_data,
    TRR_DATA_ITEMS
)
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())

TRR_HEAD_SIZE = 1000
# Actually we don't know the header size, this is just a "large"
# number so that we don't try to read the header before at least
# this number of bytes have been written.


class GromacsEngine2(GromacsEngine):
    """
    A class for interfacing GROMACS.

    This class defines an interface to GROMACS. Attributes are similar
    to :py:class:`.GromacsEngine`. In this particular interface,
    GROMACS is executed without starting and stopping and we rely on
    reading the output TRR file from GROMACS while a simulation is
    running.
    """

    def __init__(self, gmx, mdrun, input_path, timestep, subcycles,
                 exe_path=os.path.abspath('.'),
                 maxwarn=0, gmx_format='gro', write_vel=True,
                 write_force=False, domain_decomp=False):
        """Set up the GROMACS engine.

        Parameters
        ----------
        gmx : string
            The GROMACS executable.
        mdrun : string
            The GROMACS mdrun executable.
        input_path : string
            The absolute path to where the input files are stored.
        timestep : float
            The time step used in the GROMACS MD simulation.
        subcycles : integer
            The number of steps each GROMACS MD run is composed of.
        exe_path : string, optional
            The absolute path at which the main PyRETIS simulation will be run.
        maxwarn : integer, optional
            Setting for the GROMACS ``grompp -maxwarn`` option.
        gmx_format : string, optional
            The format used for GROMACS configurations.
        write_vel : boolean, optional
            Determines if GROMACS should write velocities or not.
        write_force : boolean, optional
            Determines if GROMACS should write forces or not.
        domain_decomp: boolean, optional
            Wheter GROMACS uses domain decomposition or not. This will not
            set the domain decomposition behavior, the user should know
            know this.
        """
        super().__init__(gmx, mdrun, input_path, timestep, subcycles,
                         exe_path=exe_path,
                         maxwarn=maxwarn, gmx_format=gmx_format,
                         write_vel=write_vel, write_force=write_force,
                         domain_decomp=domain_decomp)

    def _propagate_from(self, name, path, ensemble, msg_file, reverse=False):
        """
        Propagate with GROMACS from the current system configuration.

        Here, we assume that this method is called after the propagate()
        has been called in the parent. The parent is then responsible
        for reversing the velocities and also for setting the initial
        state of the system.

        Parameters
        ----------
        name : string
            A name to use for the trajectory we are generating.
        path : object like :py:class:`pyretis.core.Path.PathBase`
            This is the path we use to fill in phase-space points.
        ensemble: dict
            it contains:

            * `system` : object like :py:class:`.System`
              The system object gives the initial state for the
              integration. The initial state is stored and the system is
              reset to the initial state when the integration is done.
            * `order_function` : object like :py:class:`.OrderParameter`
              The object used for calculating the order parameter.
            * `interfaces` : list of floats
              These interfaces define the stopping criterion.

        msg_file : object like :py:class:`.FileIO`
            An object we use for writing out messages that are useful
            for inspecting the status of the current propagation.
        reverse : boolean, optional
            If True, the system will be propagated backward in time.

        Returns
        -------
        success : boolean
            This is True if we generated an acceptable path.
        status : string
            A text description of the current status of the propagation.

        """
        status = f'propagating with GROMACS (reverse = {reverse})'
        system = ensemble['system']
        interfaces = ensemble['interfaces']
        order_function = ensemble['order_function']
        logger.debug(status)
        success = False
        left, _, right = interfaces
        # Dumping of the initial config were done by the parent, here
        # we will just use it:
        initial_conf = system.particles.get_pos()[0]
        # Get the current order parameter:
        order = self.calculate_order(ensemble)
        msg_file.write(
            f'# Initial order parameter: {" ".join([str(i) for i in order])}'
        )
        # So, here we will just blast off GROMACS and check the .trr
        # output when we can.
        # 1) Create mdp_file with updated number of steps:
        settings = {'gen_vel': 'no',
                    'nsteps': path.maxlen * self.subcycles,
                    'continuation': 'no'}
        mdp_file = os.path.join(self.exe_dir, f'{name}.mdp')
        self._modify_input(self.input_files['input'], mdp_file, settings,
                           delim='=')
        # 2) Run GROMACS preprocessor:
        out_files = self._execute_grompp(mdp_file, initial_conf, name)
        # Generate some names that will be created by mdrun:
        confout = f'{name}.{self.ext}'
        out_files['conf'] = confout
        out_files['cpt_prev'] = f'{name}_prev.cpt'
        for key in ('cpt', 'edr', 'log', 'trr'):
            out_files[key] = f'{name}.{key}'
        # Remove some of these files if present (e.g. left over from a
        # crashed simulation). This is so that GromacsRunner will not
        # start reading a .trr left from a previous simulation.
        remove = [val for key, val in out_files.items() if key != 'tpr']
        self._remove_files(self.exe_dir, remove)
        tpr_file = out_files['tpr']
        trr_file = os.path.join(self.exe_dir, out_files['trr'])
        edr_file = os.path.join(self.exe_dir, out_files['edr'])
        cmd = shlex.split(self.mdrun.format(tpr_file, name, confout))
        # 3) Fire off GROMACS mdrun:
        logger.debug('Executing GROMACS.')
        msg_file.write(f'# Trajectory file is: {trr_file}')
        msg_file.write('# Starting GROMACS.')
        msg_file.write('# Step order parameter cv1 cv2 ...')
        with GromacsRunner(cmd, trr_file, edr_file, self.exe_dir) as gro:
            for i, data in enumerate(gro.get_gromacs_frames()):
                # Update the configuration file:
                system.particles.set_pos((trr_file, i))
                # Also provide the loaded positions since they are
                # available:
                system.particles.pos = data['x']
                system.particles.vel = data.get('v', None)
                if system.particles.vel is not None and reverse:
                    system.particles.vel *= -1
                length = box_matrix_to_list(data['box'])
                system.update_box(length)
                order = order_function.calculate(system)
                msg_file.write(f'{i} {" ".join([str(j) for j in order])}')
                snapshot = {'order': order,
                            'config': (trr_file, i),
                            'vel_rev': reverse}
                phase_point = self.snapshot_to_system(system, snapshot)
                status, success, stop, _ = self.add_to_path(path, phase_point,
                                                            left, right)
                if stop:
                    logger.debug('Ending propagate at %i. Reason: %s',
                                 i, status)
                    break
        logger.debug('GROMACS propagation done, obtaining energies!')
        msg_file.write('# Propagation done.')
        msg_file.write(f'# Reading energies from: {out_files["edr"]}')
        energy = self.get_energies(out_files['edr'])
        path.update_energies(energy['kinetic en.'], energy['potential'])
        logger.debug('Removing GROMACS output after propagate.')
        remove = [val for key, val in out_files.items()
                  if key not in ('trr', 'gro', 'g96')]
        self._remove_files(self.exe_dir, remove)
        self._remove_gromacs_backup_files(self.exe_dir)
        msg_file.flush()
        return success, status

    def integrate(self, ensemble, steps, thermo='full'):
        """
        Perform several integration steps.

        This method will perform several integration steps using
        GROMACS. It will also calculate order parameter(s) and energy
        terms if requested.

        Parameters
        ----------
        ensemble: dict
            it contains:

            * `system` : object like :py:class:`.System`
              The system object gives the initial state for the
              integration. The initial state is stored and the system is
              reset to the initial state when the integration is done.
            * `order_function` : object like :py:class:`.OrderParameter`
              The object used for calculating the order parameter.

        steps : integer
            The number of steps we are going to perform. Note that we
            do not integrate on the first step (e.g. step 0) but we do
            obtain the other properties. This is to output the starting
            configuration.
        thermo : string, optional
            Select the thermodynamic properties we are to obtain.

        Yields
        ------
        results : dict
            The results from a MD step. This contains the state of the system
            and order parameter(s) and energies (if calculated).

        """
        logger.debug('Integrating with GROMACS')
        # Dump the initial config:
        system = ensemble['system']
        order_function = ensemble.get('order_function')
        initial_file = self.dump_frame(system)
        self.energy_terms = self.select_energy_terms(thermo)
        if order_function:
            order = self.calculate_order(ensemble)
        else:
            order = None
        name = 'pyretis-gmx'
        # 1) Create mdp_file with updated number of steps:
        # Note the -1 here due do different numbering in GROMACS and PyRETIS.
        settings = {'nsteps': (steps - 1) * self.subcycles,
                    'continuation': 'no'}
        mdp_file = os.path.join(self.exe_dir, f'{name}.mdp')
        self._modify_input(self.input_files['input'], mdp_file, settings,
                           delim='=')
        # 2) Run GROMACS preprocessor:
        out_files = self._execute_grompp(mdp_file, initial_file, name)
        # Generate some names that will be created by mdrun:
        confout = f'{name}.{self.ext}'
        out_files['conf'] = confout
        out_files['cpt_prev'] = f'{name}_prev.cpt'
        for key in ('cpt', 'edr', 'log', 'trr'):
            out_files[key] = f'{name}.{key}'
        # Remove some of these files if present (e.g. left over from a
        # crashed simulation). This is so that GromacsRunner will not
        # start reading a .trr left from a previous simulation.
        remove = [val for key, val in out_files.items() if key != 'tpr']
        self._remove_files(self.exe_dir, remove)
        tpr_file = out_files['tpr']
        trr_file = os.path.join(self.exe_dir, out_files['trr'])
        edr_file = os.path.join(self.exe_dir, out_files['edr'])
        cmd = shlex.split(self.mdrun.format(tpr_file, name, confout))
        # 3) Fire off GROMACS mdrun:
        logger.debug('Executing GROMACS.')
        with GromacsRunner(cmd, trr_file, edr_file, self.exe_dir) as gro:
            for i, data in enumerate(gro.get_gromacs_frames()):
                system.particles.pos = data['x']
                system.particles.vel = data.get('v', None)
                length = box_matrix_to_list(data['box'])
                system.update_box(length)
                results = {}
                if order:
                    results['order'] = order
                if order_function:
                    order = order_function.calculate(system)
                time1 = (i * self.timestep * self.subcycles -
                         0.1 * self.timestep)
                time2 = ((i + 1) * self.timestep * self.subcycles +
                         0.1 * self.timestep)
                energy = self.get_energies(out_files['edr'], begin=time1,
                                           end=time2)
                results['thermo'] = self.rename_energies(energy)
                yield results
        logger.debug('GROMACS execution done.')


def get_data(fileh, header):
    """Read data from the TRR file.

    Parameters
    ----------
    fileh : file object
        The file we are reading.
    header : dict
        The previously read header. Contains sizes and what to read.

    Returns
    -------
    data : dict
        The data read from the file.
    data_size : integer
        The size of the data read.

    """
    data_size = sum(header[key] for key in TRR_DATA_ITEMS)
    data = read_trr_data(fileh, header)
    return data, data_size


def reopen_file(filename, fileh, inode, bytes_read):
    """Reopen a file if the inode has changed.

    Parameters
    ----------
    filename : string
        The name of the file we are working with.
    fileh : file object
        The current open file object.
    inode : integer
        The current inode we are using.
    bytes_read : integer
        The position we should start reading at.

    Returns
    -------
    out[0] : file object or None
        The new file object.
    out[1] : integer or None
        The new inode.

    """
    if os.stat(filename).st_ino != inode:
        new_fileh = open(filename, 'rb')  # pylint: disable=consider-using-with
        fileh.close()
        new_inode = os.fstat(new_fileh.fileno()).st_ino
        new_fileh.seek(bytes_read)
        return new_fileh, new_inode
    return None, None


def read_remaining_trr(filename, fileh, start):
    """Read remaining frames from the TRR file.

    Parameters
    ----------
    filename : string
        The file we are reading from.
    fileh : file object
        The file object we are reading from.
    start : integer
        The current position we are at.

    Yields
    ------
    out[0] : string
        The header read from the file
    out[1] : dict
        The data read from the file.
    out[2] : integer
        The size of the data read.

    """
    stop = False
    bytes_read = start
    bytes_total = os.path.getsize(filename)
    logger.debug('Reading remaing data from: %s', filename)
    while not stop:
        if bytes_read >= bytes_total:
            stop = True
            continue
        header = None
        new_bytes = bytes_read
        try:
            header, new_bytes = read_trr_header(fileh)
        except EOFError:  # pragma: no cover
            # Just assume that we have reached the end of the
            # file and we just stop here. It should not be reached,
            # kept for safety
            stop = True
            continue
        if header is not None:
            bytes_read += new_bytes
            try:
                data, new_bytes = get_data(fileh, header)
                if data is not None:
                    bytes_read += new_bytes
                    yield header, data, bytes_read
            except EOFError:  # pragma: no cover
                # Hopefully, this code should not be reached.
                # kept for safety
                stop = True
                continue


class GromacsRunner:
    """A helper class for running GROMACS.

    This class handles the reading of the TRR on the fly and
    it is used to decide when to end the GROMACS execution.

    Attributes
    ----------
    cmd : string
        The command for executing GROMACS.
    trr_file : string
        The GROMACS TRR file we are going to read.
    edr_file : string
        A .edr file we are going to read.
    exe_dir : string
        Path to where we are currently running GROMACS.
    fileh : file object
        The current open file object.
    running : None or object like :py:class:`subprocess.Popen`
        The process running GROMACS.
    bytes_read : integer
        The number of bytes read so far from the TRR file.
    ino : integer
        The current inode we are using for the file.
    stop_read : boolean
        If this is set to True, we will stop the reading.
    SLEEP : float
        How long we wait after an unsuccessful read before
        reading again.
    data_size : integer
        The size of the data (x, v, f, box, etc.) in the TRR file.
    header_size : integer
        The size of the header in the TRR file.

    """

    SLEEP = 0.1

    def __init__(self, cmd, trr_file, edr_file, exe_dir):
        """Set the GROMACS command and the files we need.

        Parameters
        ----------
        cmd : string
            The command for executing GROMACS.
        trr_file : string
            The GROMACS TRR file we are going to read.
        edr_file : string
            A .edr file we are going to read.
        exe_dir : string
            Path to where we are currently running GROMACS.

        """
        self.cmd = cmd
        self.trr_file = trr_file
        self.edr_file = edr_file
        self.exe_dir = exe_dir
        self.fileh = None
        self.running = None
        self.bytes_read = 0
        self.ino = 0
        self.stop_read = True
        self.data_size = 0
        self.header_size = 0
        self.stdout_name = None
        self.stderr_name = None
        self.stdout = None
        self.stderr = None

    def start(self):
        """Start execution of GROMACS and wait for output file creation."""
        logger.debug('Starting GROMACS execution in %s', self.exe_dir)

        self.stdout_name, self.stderr_name = \
            ExternalMDEngine._engine_output_files(self.exe_dir)
        # pylint: disable=consider-using-with
        self.stdout = open(self.stdout_name, 'ab')
        # pylint: disable=consider-using-with
        self.stderr = open(self.stderr_name, 'ab')

        self.running = subprocess.Popen(  # nosec B603 B607
            self.cmd,
            stdin=subprocess.PIPE,
            stdout=self.stdout,
            stderr=self.stderr,
            shell=False,
            cwd=self.exe_dir,
            start_new_session=True,
        )
        present = []
        # Wait for the TRR/EDR files to appear:
        for fname in (self.trr_file, self.edr_file):
            while not os.path.isfile(fname):
                logger.debug('Waiting for GROMACS file "%s"', fname)
                sleep(self.SLEEP)
                poll = self.check_poll()
                if poll is not None:
                    logger.debug('GROMACS execution stopped')
                    break
            if os.path.isfile(fname):
                present.append(fname)
        # Prepare and open the TRR file:
        self.bytes_read = 0
        # Ok, so GROMACS might have crashed in between writing the
        # files. Check that both files are indeed here:
        if self.trr_file in present and self.edr_file in present:
            # pylint: disable=consider-using-with
            self.fileh = open(self.trr_file, 'rb')
            self.ino = os.fstat(self.fileh.fileno()).st_ino
            self.stop_read = False
        else:
            self.stop_read = True

    def __enter__(self):
        """Start running GROMACS, for a context manager."""
        self.start()
        return self

    def get_gromacs_frames(self):
        """Read the GROMACS TRR file on-the-fly."""
        first_header = True
        header = None
        while not self.stop_read:
            poll = self.check_poll()
            if poll is not None:
                # GROMACS is done, read remaining data.
                self.stop_read = True
                if os.path.getsize(self.trr_file) - self.bytes_read > 0:
                    for _, data, _ in read_remaining_trr(self.trr_file,
                                                         self.fileh,
                                                         self.bytes_read):
                        yield data

            else:
                # First we try to get the header from the file:
                size = os.path.getsize(self.trr_file)
                if self.header_size == 0:
                    header_size = TRR_HEAD_SIZE
                else:
                    header_size = self.header_size
                if size >= self.bytes_read + header_size:
                    # Try to read next frame:
                    try:
                        header, new_bytes = read_trr_header(self.fileh)
                    except EOFError:
                        new_fileh, new_ino = reopen_file(self.trr_file,
                                                         self.fileh,
                                                         self.ino,
                                                         self.bytes_read)
                        if new_fileh is not None:
                            self.fileh = new_fileh
                            self.ino = new_ino
                    if header is not None:
                        self.bytes_read += new_bytes
                        self.header_size = new_bytes
                        if first_header:
                            logger.debug('TRR header was: %i', new_bytes)
                            first_header = False
                        # Calculate the size of the data:
                        self.data_size = sum(header[key] for key in
                                             TRR_DATA_ITEMS)
                        data = None
                        while data is None:
                            size = os.path.getsize(self.trr_file)
                            if size >= self.bytes_read + self.data_size:
                                try:
                                    data, new_bytes = get_data(self.fileh,
                                                               header)
                                except EOFError:
                                    new_fileh, new_ino = reopen_file(
                                        self.trr_file,
                                        self.fileh,
                                        self.ino,
                                        self.bytes_read)
                                    if new_fileh is not None:
                                        self.fileh = new_fileh
                                        self.ino = new_ino
                                if data is None:
                                    # Data is not ready, just wait:
                                    sleep(self.SLEEP)
                                else:
                                    self.bytes_read += new_bytes
                                    yield data
                            else:
                                # Data is not ready, just wait:
                                sleep(self.SLEEP)
                else:
                    # Header was not ready, just wait before trying again.
                    sleep(self.SLEEP)

    def close(self):
        """Close the file, in case that is explicitly needed."""
        if self.fileh is not None and not self.fileh.closed:
            logger.debug('Closing GROMACS file: "%s"', self.trr_file)
            self.fileh.close()
        for handle in (self.stdout, self.stderr):
            if handle is not None and not handle.closed:
                handle.close()
        if self.running is not None:
            ExternalMDEngine._finalize_engine_output(
                self.running.returncode, self.stderr_name
            )

    def stop(self):
        """Stop the current GROMACS execution."""
        if self.running:
            for handle in (self.running.stdin, self.running.stdout,
                           self.running.stderr):
                if handle:
                    try:
                        handle.close()
                    except AttributeError:
                        pass
            if self.running.returncode is None:
                logger.debug('Terminating GROMACS execution')
                os.killpg(os.getpgid(self.running.pid), signal.SIGTERM)

            self.running.wait(timeout=360)
        self.stop_read = True
        self.close()  # Close the TRR file.

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Just stop execution and close file for a context manager."""
        self.stop()

    def __del__(self):
        """Just stop execution and close file."""
        self.stop()

    def check_poll(self):
        """Check the current status of the running subprocess."""
        if self.running:
            poll = self.running.poll()
            if poll is not None:
                logger.debug('Execution of GROMACS stopped')
                logger.debug('Return code was: %i', poll)
                if poll != 0:
                    logger.error('STDOUT, see file: %s', self.stdout_name)
                    logger.error('STDERR, see file: %s', self.stderr_name)
                    raise RuntimeError('Error in GROMACS execution.')
            return poll
        raise RuntimeError('GROMACS is not running.')

# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""A template for a new external MD engine."""
import logging
import os
from pyretis.engines.external import ExternalMDEngine


logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())

INPUT_CONF = []
OUTPUT_FILES = []
REQUIRED_FILES = []


def get_energy():
    """Fictional function to read energy."""
    return 42, 42


def kinetic_energy_output_from_velocity_generation(_, __):
    """Fictional function to read kin eng."""
    return 42


def method_to_read_config_file(_):
    """Fictional function to read config file."""
    return 42, 42


def method_to_write_config_file(_, __, ___):
    """Fictional function to write a config file."""
    return


class NewEngine(ExternalMDEngine):
    """A new external engine class.

    Attributes
    ----------
    exe : string
        The command for executing the engine.
    input_path : string
        The directory where the input files are stored.
    timestep : float
        The time step used in the external MD simulation.
    subcycles : integer
        The number of steps each external step is composed of.
    """

    def __init__(self, exe, input_path, timestep, subcycles):
        """Initiate the script.

        Parameters
        ----------
        exe : string
            The engine executable.
        input_path : string
            The absolute path to where the input files are stored.
        timestep : float
            The time step used in the external engine.
        subcycles : integer
            The number of steps each external step is composed of.
        """
        super().__init__('New external engine', timestep, subcycles)
        # Store variables for the engine:
        self.exe = exe
        # store input path:
        self.input_path = os.path.abspath(input_path)
        # store input files:
        self.input_files = {}
        for key in REQUIRED_FILES:
            self.input_files[key] = os.path.join(self.input_path, key)
            if not os.path.isfile(self.input_files[key]):
                msg = 'Engine could not find file "{}"!'.format(key)
                raise ValueError(msg)

        # Create a pyretis version of the input file:
        settings = {'timestep': self.timestep,
                    'subcycles': self.subcycles}

        self.input_files['input'] = os.path.join(self.input_path,
                                                 'PYRETIS_INPUT')
        self._modify_input(self.input_files['ORIGINAL_INPUT'],
                           self.input_files['input'], settings, delim='=')

    def _run_program(self):
        """Method to execute the new engine.

        Returns
        -------
        out : dict
            The files created by the run.
        """
        self.execute_command(self.exe, cwd=self.exe_dir, inputs=None)
        out = {}
        for key in OUTPUT_FILES:
            filename = os.path.join(self.exe_dir, key)
            if os.path.isfile(filename):
                out[key] = filename
        return out

    def _extract_frame(self, traj_file, idx, out_file):
        """Extract a frame from a trajectory file.

        This method is used by `self.dump_config` when we are
        dumping from a trajectory file. It is not used if we are
        dumping from a single config file. Here you can write
        your own code for dumping or make use of existing libraries
        such as mdtraj, mdanalysis, ase, ...

        Parameters
        ----------
        traj_file : string
            The trajectory file to dump from.
        idx : integer
            The frame number we look for.
        out_file : string
            The file to dump to.
        """
        return

    def _propagate_from(self, name, path, ensemble, reverse=False):
        """Propagate from the current system configuration.

        Here, we assume that this method is called after the propagate()
        has been called in the parent. The parent is then responsible
        for reversing the velocities and also for setting the initial
        state of the system.

        Parameters
        ----------
        name : string
            A name to use for the trajectory we are generating.
        path : object like :py:class:`.PathBase`
            This is the path we use to fill in phase-space points.
            We are here not returning a new path - this since we want
            to delegate the creation of the path to the method
            that is running `propagate`.
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
        reverse : boolean
            If True, the system will be propagated backwards in time.

        Returns
        -------
        success : boolean
            This is True if we generated an acceptable path.
        status : string
            A text description of the current status of the propagation.
        """
        # Dumping of the initial config were done by the parent, here
        # we will just the dumped file:
        initial_conf = ensemble['system'].particles.get_pos()[0]
        # Copy this file, in case the external program requires
        # specific names for the input configurations:
        input_conf = os.path.join(self.exe_dir, 'input_conf')
        output_conf = os.path.join(self.exe_dir, 'output_conf')
        self._copyfile(initial_conf, output_conf)
        # Set the positions to be just the file we copied:
        ensemble['system'].particles.set_pos((input_conf, None))
        # Create the output trajectory file:
        traj_file = os.path.join(self.exe_dir, name)
        # Add other files here as well, input settings, force field etc...
        success = False
        left, _, right = ensemble['interfaces']
        phase_point = ensemble['system'].particles.get_particle_state()
        # Just to be keep the same order as in the loop
        vpot = phase_point['vpot']
        ekin = phase_point['ekin']
        order = self.calculate_order(ensemble)
        for i in range(path.maxlen):
            # We first add the current phase point, and then we propagate.
            phase_point = {
                'order': order,
                'pos': (traj_file, i),
                'vel': reverse,
                'vpot': vpot,
                'ekin': ekin}
            status, success, stop = self.add_to_path(path, phase_point,
                                                     left, right)
            # Add the current frame to the output trajectory also:
            # e.g.: add_frame_to_traj('CURRENT_FRAME_FILE', traj_file, i)
            if stop:
                logger.debug('Ending propagate at %i. Reason: %s', i, status)
                break
            # Here execute the external program
            # for instance by using self.step() or self._run_program()
            # get energies after execution:
            ekin, vpot = get_energy()
            # Recalculate order parameter
            # If the output file is named different, update the
            # system object so that it points to this file, e.g.:
            # system.particles.set_pos((NAME_OF_FILE, None)), or
            # if the energies are needed as well, use
            # system.particles.set_particle_state(...)
            order = self.calculate_order(ensemble)
        # We are done with the execution:
        # clean up and remove files:
        list_of_files_to_remove = []
        self._remove_files(self.exe_dir, list_of_files_to_remove)
        return success, status

    def step(self, system, name):
        """Perform a single step with the external engine.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The system we are integrating.
        name : string
            To name the output files from the step.

        Returns
        -------
        out_files : dict
            The output files created by this step.
        """
        # Get initial config:
        input_file = self.dump_frame(system)
        # copy input files, maybe rename input_file etc.
        # Then, run it
        out_files = self._run_program()
        # Store the frame after the step, make use of the name variable:
        config = os.path.join(self.exe_dir, 'OUTPUT.{}'.format(name))
        # Get energies
        ekin, vpot = get_energy()
        # remove files that we dont need after the step:
        files_to_remove = []
        self._remove_files(self.exe_dir, files_to_remove)
        # Update the state of the system to match the current config
        phase_point = {'pos': (config, None),
                       'vel': False, 'vpot': vpot, 'ekin': ekin}
        system.particles.set_particle_state(phase_point)
        return config

    @staticmethod
    def _read_configuration(filename):
        """Method to read config files for the external engine.

        This method is used when we calculate the order parameter.

        Parameters
        ----------
        filename : string
            The file to read the configuration from.

        Returns
        -------
        xyz : numpy.array
            The positions.
        vel : numpy.array
            The velocities.
        """
        xyz, vel = method_to_read_config_file(filename)
        return xyz, vel

    @staticmethod
    def _reverse_velocities(filename, outfile):
        """Method to reverse velocity in a given snapshot.

        Parameters
        ----------
        filename : string
            The configuration to reverse velocities in.
        outfile : string
            The output file for storing the configuration with
            reversed velocities.
        """
        xyz, vel = method_to_read_config_file(filename)
        method_to_write_config_file(outfile, xyz, -vel)

    def modify_velocities(self, system, rgen, vel_settings):
        """Modify the velocities of the current state.

        This method will modify the velocities of a time slice.

        Parameters
        ----------
        system : object like :py:class:`.System`
            System is used here since we need access to the particle
            list.
        rgen : object like :py:class:`.RandomGenerator`
            This is the random generator that will be used.
        vel_settings: dict
            It contains all the info for the velocity:

            * `sigma_v` : numpy.array, optional
              These values can be used to set a standard deviation (one
              for each particle) for the generated velocities.
            * `aimless` : boolean, optional
              Determines if we should do aimless shooting or not.
            * `momentum` : boolean, optional
              If True, we reset the linear momentum to zero after
              generating.
            * `rescale` : float, optional
              In some NVE simulations, we may wish to re-scale the
              energy to a fixed value. If `rescale` is a float > 0,
              we will re-scale the energy (after modification of
              the velocities) to match the given float.

        Returns
        -------
        dek : float
            The change in the kinetic energy.
        kin_new : float
            The new kinetic energy.
        """
        dek = None
        kin_old = None
        kin_new = None
        rescale = vel_settings.get('rescale')

        if rescale is not None and rescale is not False and rescale > 0:
            # Code to support energy rescale, or just ignore this as
            # set:
            msgtxt = 'External engine does not support energy rescale!'
            logger.error(msgtxt)
            raise NotImplementedError(msgtxt)
        else:
            kin_old = system.particles.ekin
        if vel_settings.get('aimless', False):
            # pos = self.dump_frame(system)
            phase_point = system.particles.get_particle_state()
            # Generate new velocities, either by hand or by asking the
            # external program.
            genvel = os.path.join(self.exe_dir, 'g_POSCAR')
            # Also get the new kinetic energy
            kin_new = kinetic_energy_output_from_velocity_generation(genvel,
                                                                     rgen)
            # Update the phase point (note we keep the potential energy):
            phase_point['pos'] = (genvel, None)
            phase_point['vel'] = False
            phase_point['ekin'] = kin_new
            system.particles.set_particle_state(phase_point)
        else:  # soft velocity change, add from Gaussian dist
            # You might not need to support this:
            msgtxt = 'External engine only support aimless shooting!'
            logger.error(msgtxt)
            raise NotImplementedError(msgtxt)
        if kin_old is None or kin_new is None:
            dek = float('inf')
            logger.warning(('Kinetic energy not found for previous point.'
                            '\n(This happens when the initial configuration '
                            'does not contain energies.)'))
        else:
            dek = kin_new - kin_old
        return dek, kin_new

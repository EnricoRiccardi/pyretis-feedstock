# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""A GROMACS external MD integrator interface.

This module defines a class for using GROMACS as an external engine.

Important classes defined here
------------------------------

GromacsEngine (:py:class:`.GromacsEngine`)
    A class responsible for interfacing GROMACS.
"""
import logging
import os
import shlex
from pyretis.core.box import box_matrix_to_list
from pyretis.engines.external import ExternalMDEngine
from pyretis.inout.settings import look_for_input_files
from pyretis.inout.formats.gromacs import (
    read_gromos96_file,
    read_gromacs_gro_file,
    write_gromacs_gro_file,
    write_gromos96_file,
    read_xvg_file,
    read_trr_frame,
)
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


class GromacsEngine(ExternalMDEngine):
    """
    A class for interfacing GROMACS.

    This class defines the interface to GROMACS.

    Attributes
    ----------
    gmx : string
        The command for executing GROMACS. Note that we are assuming
        that we are using version 5 (or later) of GROMACS.
    mdrun : string
        The command for executing GROMACS mdrun. In some cases, this
        executable can be different from ``gmx mdrun``.
    mdrun_c : string
        The command for executing GROMACS mdrun when continuing a
        simulation. This is derived from the ``mdrun`` command.
    input_path : string
        The directory where the input files are stored.
    subcycles : int,
        The number of simulation steps of the external engine for each
        PyRETIS step (e.g. interaction between the softwares frequency)
    exe_path : string, optional
        The absolute path at which the main PyRETIS simulation will be run.
    maxwarn : integer
        Setting for the GROMACS ``grompp -maxwarn`` option.
    gmx_format : string
        This string selects the output format for GROMACS.
    write_vel : boolean, optional
        True if we want to output the velocities.
    write_force : boolean, optional
        True if we want to output the forces.
    domain_decomp: boolean, optional
        Whether gromacs uses domain decomposition (DD) or not. This
        does not set domain decomposition. The user should know whether
        DD is executed or not.
    """

    def __init__(self, gmx, mdrun, input_path, timestep, subcycles,
                 exe_path=os.path.abspath('.'),
                 maxwarn=0, gmx_format='gro',
                 write_vel=True,
                 write_force=False,
                 domain_decomp=False,
                 conf=None, template='grompp.mdp', top='topol.top'):
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
            Whether domain decomposition (DD) is executed by GROMACS or not.
            User should know if this happens or not.
        conf : string, optional
            The default configuration file (e.g. 'conf.gro').
        template : string, optional
            The GROMACS mdp template (default is 'grompp.mdp').
        top : string, optional
            The GROMACS topology file (default is 'topol.top').
        """
        super().__init__('GROMACS engine', timestep, subcycles)
        self.ext = gmx_format
        if self.ext not in ('g96', 'gro'):
            msg = 'Unknown GROMACS format: "%s"'
            logger.error(msg, self.ext)
            raise ValueError(msg % self.ext)
        # Define the GROMACS GMX command:
        self.gmx = gmx
        # Define GROMACS GMX MDRUN commands:
        self.mdrun = mdrun + ' -s {} -deffnm {} -c {}'
        # This is for continuation of a GROMACS simulation:
        self.mdrun_c = mdrun + ' -s {} -cpi {} -append -deffnm {} -c {}'
        # Whether DD is used or not:
        self.domain_decomp = domain_decomp
        self.ext_time = self.timestep * self.subcycles
        self.maxwarn = maxwarn
        # Define the energy terms, these are hard-coded, but
        # here we open up for changing that:
        self.energy_terms = self.select_energy_terms('path')
        self.input_path = os.path.join(exe_path, input_path)
        # Set the defaults input files:
        if conf is None:
            conf = f'conf.{self.ext}'
        default_files = {
            'conf': conf,
            'input_o': template,  # "o" = original input file.
            'topology': top}
        extra_files = {
            'index': 'index.ndx',
        }

        # An user doesn't need to have problems with g96 and mdtraj.
        file_g = os.path.join(self.input_path, os.path.splitext(conf)[0])
        if self.ext == 'gro':
            self.top, _, _, _ = read_gromacs_gro_file(
                os.path.join(self.input_path, conf)
            )
        elif self.ext == 'g96':
            if not os.path.isfile(file_g+'.gro'):
                cmd = [self.gmx, 'editconf',
                       '-f', os.path.join(self.input_path, conf),
                       '-o', file_g+'.gro']
                self.execute_command(cmd, cwd=None)

            self.top, _, _, _ = read_gromos96_file(
                os.path.join(self.input_path, conf)
            )
            self.top['VELOCITY'] = self.top['POSITION'].copy()

        # Check the presence of the defaults input files or, if absent,
        # try to find them by their expected extension.
        self.input_files = look_for_input_files(self.input_path,
                                                default_files,
                                                extra_files)
        logger.debug('Input files: %s', self.input_files)
        # Check the input file and create a PyRETIS version with
        # consistent settings:
        settings = {
            'dt': self.timestep,
            'nstxout-compressed': 0,
            'gen_vel': 'no'
        }
        for key in ('nsteps', 'nstxout', 'nstvout', 'nstfout', 'nstlog',
                    'nstcalcenergy', 'nstenergy'):
            settings[key] = self.subcycles
        if not write_vel:
            settings['nstvout'] = 0
        if not write_force:
            settings['nstfout'] = 0

        # PyRETIS construct its own mdp file
        self.input_files['input'] = os.path.join(self.input_path,
                                                 'pyretis.mdp')
        self._modify_input(self.input_files['input_o'],
                           self.input_files['input'], settings, delim='=')
        logger.info(('Created GROMACS mdp input from %s. You might '
                     'want to check the input file: %s'),
                    self.input_files['input_o'], self.input_files['input'])

        # Generate a tpr file using the input files:
        logger.info('Creating ".tpr" for GROMACS in %s', self.input_path)
        self.exe_dir = self.input_path
        out_files = self._execute_grompp(self.input_files['input'],
                                         self.input_files['conf'], 'topol')

        # This will generate some noise, let's remove files we don't need:
        mdout = os.path.join(self.input_path, out_files['mdout'])
        self._removefile(mdout)
        # We also remove GROMACS backup files after creating the tpr:
        self._remove_gromacs_backup_files(self.input_path)
        # Keep the tpr file.
        self.input_files['tpr'] = os.path.join(self.input_path,
                                               out_files['tpr'])
        logger.info('GROMACS ".tpr" created: %s', self.input_files['tpr'])

    @staticmethod
    def select_energy_terms(terms):
        """Select energy terms to extract from GROMACS.

        Parameters
        ----------
        terms : string
            This string will name the terms to extract. Currently
            we only allow for two types of output, but this can be
            customized in the future.

        """
        allowed_terms = {
            'full': ('\n'.join(('Potential', 'Kinetic-En.', 'Total-Energy',
                                'Temperature', 'Pressure'))).encode(),
            'path': b'Potential\nKinetic-En.',
        }
        if terms not in allowed_terms:
            return allowed_terms['path']
        return allowed_terms[terms]

    @staticmethod
    def rename_energies(gmx_energy):
        """Rename GROMACS energy terms to PyRETIS convention."""
        energy_map = {'potential': 'vpot',
                      'kinetic en.': 'ekin',
                      'temperature': 'temp',
                      'total energy': 'etot',
                      'pressure': 'press'}
        energy = {}
        for key, val in gmx_energy.items():
            name = energy_map.get(key, key)
            energy[name] = val[0]
        return energy

    def _execute_grompp(self, mdp_file, config, deffnm):
        """Execute the GROMACS preprocessor.

        Parameters
        ----------
        mdp_file : string
            The path to the mdp file.
        config : string
            The path to the GROMACS config file to use as input.
        deffnm : string
            A string used to name the GROMACS files.

        Returns
        -------
        out_files : dict
            This dict contains files that were created by the GROMACS
            preprocessor.

        """
        topol = self.input_files['topology']
        tpr = f'{deffnm}.tpr'
        cmd = [self.gmx, 'grompp', '-f', mdp_file, '-c', config,
               '-p', topol, '-o', tpr]
        if 'index' in self.input_files:
            cmd.extend(['-n', self.input_files['index']])
        if self.maxwarn > 0:
            cmd.extend(['-maxwarn', str(self.maxwarn)])
        self.execute_command(cmd, cwd=self.exe_dir)
        out_files = {'tpr': tpr, 'mdout': 'mdout.mdp'}
        return out_files

    def _execute_mdrun(self, tprfile, deffnm):
        """
        Execute GROMACS mdrun.

        This method is intended as the initial ``gmx mdrun`` executed.
        That is, we here assume that we do not continue a simulation.

        Parameters
        ----------
        tprfile : string
            The .tpr file to use for executing GROMACS.
        deffnm : string
            To give the GROMACS simulation a name.

        Returns
        -------
        out_files : dict
            This dict contains the output files created by ``mdrun``.
            Note that we here hard code the file names.

        """
        confout = f'{deffnm}.{self.ext}'
        cmd = shlex.split(self.mdrun.format(tprfile, deffnm, confout))
        self.execute_command(cmd, cwd=self.exe_dir)
        out_files = {'conf': confout,
                     'cpt_prev': f'{deffnm}_prev.cpt'}
        for key in ('cpt', 'edr', 'log', 'trr'):
            out_files[key] = f'{deffnm}.{key}'
        self._remove_gromacs_backup_files(self.exe_dir)
        return out_files

    def _execute_grompp_and_mdrun(self, config, deffnm):
        """
        Execute GROMACS ``grompp`` and ``mdrun``.

        Here we use the input file given in the input directory.

        Parameters
        ----------
        config : string
            The path to the GROMACS config file to use as input.
        deffnm : string
            A string used to name the GROMACS files.

        Returns
        -------
        out_files : dict of strings
            The files created by this command.

        """
        out_files = {}
        out_grompp = self._execute_grompp(self.input_files['input'],
                                          config, deffnm)
        tpr_file = out_grompp['tpr']
        for key, value in out_grompp.items():
            out_files[key] = value
        out_mdrun = self._execute_mdrun(tpr_file,
                                        deffnm)
        for key, value in out_mdrun.items():
            out_files[key] = value
        return out_files

    def _execute_mdrun_continue(self, tprfile, cptfile, deffnm):
        """
        Continue the execution of GROMACS.

        Here, we assume that we have already executed ``gmx mdrun`` and
        that we are to append and continue a simulation.

        Parameters
        ----------
        tprfile : string
            The .tpr file which defines the simulation.
        cptfile : string
            The last checkpoint file (.cpt) from the previous run.
        deffnm : string
            To give the GROMACS simulation a name.

        Returns
        -------
        out_files : dict
            The output files created/appended by GROMACS when we
            continue the simulation.

        """
        confout = f'{deffnm}.{self.ext}'.format(deffnm, self.ext)
        self._removefile(confout)
        cmd = shlex.split(self.mdrun_c.format(tprfile, cptfile,
                                              deffnm, confout))
        self.execute_command(cmd, cwd=self.exe_dir)
        out_files = {'conf': confout}
        for key in ('cpt', 'edr', 'log', 'trr'):
            out_files[key] = f'{deffnm}.{key}'
        self._remove_gromacs_backup_files(self.exe_dir)
        return out_files

    def _extend_gromacs(self, tprfile, time):
        """Extend a GROMACS simulation.

        Parameters
        ----------
        tprfile : string
            The file to read for extending.
        time : float
            The time (in ps) to extend the simulation by.

        Returns
        -------
        out_files : dict
            The files created by GROMACS when we extend.

        """
        tpxout = f'ext_{tprfile}'
        self._removefile(tpxout)
        cmd = [self.gmx, 'convert-tpr', '-s', tprfile,
               '-extend', str(time), '-o', tpxout]
        self.execute_command(cmd, cwd=self.exe_dir)
        out_files = {'tpr': tpxout}
        return out_files

    def _extend_and_execute_mdrun(self, tpr_file, cpt_file, deffnm):
        """Extend GROMACS and execute mdrun.

        Parameters
        ----------
        tpr_file : string
            The location of the "current" .tpr file.
        cpt_file : string
            The last checkpoint file (.cpt) from the previous run.
        deffnm : string
            To give the GROMACS simulation a name.

        Returns
        -------
        out_files : dict
            The files created by GROMACS when we extend.

        """
        out_files = {}
        out_grompp = self._extend_gromacs(tpr_file, self.ext_time)
        ext_tpr_file = out_grompp['tpr']
        for key, value in out_grompp.items():
            out_files[key] = value
        out_mdrun = self._execute_mdrun_continue(ext_tpr_file, cpt_file,
                                                 deffnm)
        for key, value in out_mdrun.items():
            out_files[key] = value
        # Move extended tpr so that we can continue extending:
        source = os.path.join(self.exe_dir, ext_tpr_file)
        dest = os.path.join(self.exe_dir, tpr_file)
        self._movefile(source, dest)
        out_files['tpr'] = tpr_file
        return out_files

    def _remove_gromacs_backup_files(self, dirname):
        """Remove files GROMACS has backed up.

        These are files starting with a '#'

        Parameters
        ----------
        dirname : string
            The directory where we are to remove files.

        """
        for entry in os.scandir(dirname):
            if entry.name.startswith('#') and entry.is_file():
                filename = os.path.join(dirname, entry.name)
                self._removefile(filename)

    def _extract_frame(self, traj_file, idx, out_file):
        """Extract a frame from a .trr, .xtc or .trj file.

        If the extension is different from .trr, .xtc or .trj, we will
        basically just copy the given input file.

        Parameters
        ----------
        traj_file : string
            The GROMACS file to open.
        idx : integer
            The frame number we look for.
        out_file : string
            The file to dump to.

        Note
        ----
        This will only properly work if the frames in the input
        trajectory are uniformly spaced in time.

        """
        trajexts = ['.trr', '.xtc', '.trj']

        logger.debug('Extracting frame, idx = %i', idx)
        logger.debug('Source file: %s, out file: %s', traj_file, out_file)
        if traj_file[-4:] in trajexts:
            _, data = read_trr_frame(traj_file, idx)
            xyz = data['x']
            vel = data.get('v')
            box = box_matrix_to_list(data['box'], full=True)
            if out_file[-4:] == '.gro':
                write_gromacs_gro_file(out_file, self.top, xyz, vel, box)
            elif out_file[-4:] == '.g96':
                write_gromos96_file(out_file, self.top, xyz, vel, box)

        else:
            cmd = [self.gmx, 'trjconv',
                   '-f', traj_file,
                   '-s', self.input_files['tpr'],
                   '-o', out_file]

            self.execute_command(cmd, inputs=b'0', cwd=None)

        # If domain decomposition is used, the molecules must be made
        # whole. This is a bug in GROMACS, where DD only works auto-
        # matically for initial/ending .gro files of a simulation
        # trajectory. I.e. only for the latter .gro 'frames', the bonds
        # are correctly represented over the periodic boundaries.
        # This is version specific and *could* happen only for >16cpus
        # Therefore, we don't test this in coverage...
        if self.domain_decomp:  # pragma: no cover
            cmd = [self.gmx, 'trjconv',
                   '-f', out_file,
                   '-s', self.input_files['tpr'],
                   '-pbc', 'whole',
                   '-o', out_file]
            logger.info(("We are making the molecules whole"))
            self.execute_command(cmd, inputs=b'0', cwd=None)

    def get_energies(self, energy_file, begin=None, end=None):
        """Return energies from a GROMACS run.

        Parameters
        ----------
        energy_file : string
            The file to read energies from.
        begin : float, optional
            Select the time for the first frame to read.
        end : float, optional
            Select the time for the last frame to read.

        Returns
        -------
        energy : dict fo numpy.arrays
            The energies read from the produced GROMACS xvg file.

        """
        cmd = [self.gmx, 'energy', '-f', energy_file]
        if begin is not None:
            begin = max(begin, 0)
            cmd.extend(['-b', str(begin)])
        if end is not None:
            cmd.extend(['-e', str(end)])
        self.execute_command(cmd, inputs=self.energy_terms,
                             cwd=self.exe_dir)
        xvg_file = os.path.join(self.exe_dir, 'energy.xvg')
        energy = read_xvg_file(xvg_file)
        self._removefile(xvg_file)
        return energy

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
        path : object like :py:class:`pyretis.core.path.PathBase`
            This is the path we use to fill in phase-space points.
        ensemble: dict
            It contains:

            * `system`: object like :py:class:`.System`
              The system object gives the initial state for the
              integration. The initial state is stored and the system is
              reset to the initial state when the integration is done.
            * `order_function`: object like :py:class:`.OrderParameter`
              The object used for calculating the order parameter.
            * `interfaces`: list of floats
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
        logger.debug(status)
        success = False
        interfaces = ensemble['interfaces']
        left, _, right = interfaces
        # Dumping of the initial config were done by the parent, here
        # we will just use it:
        initial_conf = ensemble['system'].particles.get_pos()[0]
        # Get the current order parameter:
        order = self.calculate_order(ensemble)
        msg_file.write(
            f'# Initial order parameter: {" ".join([str(i) for i in order])}'
        )
        # In some cases, we don't really have to perform a step as the
        # initial config might be left/right of the interface in
        # question. Here, we will perform a step anyway. This is to be
        # sure that we obtain energies and also a trajectory segment.
        # Note that all the energies are obtained after we are done
        # with the integration from the .edr file of the trajectory.
        msg_file.write('# Running grompp and mdrun (initial step).')
        out_files = self._execute_grompp_and_mdrun(initial_conf, name)
        # Define name of some files:
        tpr_file = out_files['tpr']
        cpt_file = out_files['cpt']
        traj_file = os.path.join(self.exe_dir, out_files['trr'])
        msg_file.write(f'# Trajectory file is: {traj_file}')
        conf_abs = os.path.join(self.exe_dir, out_files['conf'])
        # Note: The order parameter is calculated AT THE END of each iteration.
        msg_file.write('# Starting GROMACS.')
        msg_file.write('# Step order parameter cv1 cv2 ...')
        for i in range(path.maxlen):
            msg_file.write(
                f'{i} {" ".join([str(j) for j in order])}'
            )
            # We first add the previous phase point, and then we propagate.
            snapshot = {'order': order,
                        'config': (traj_file, i),
                        'vel_rev': reverse}
            phase_point = self.snapshot_to_system(ensemble['system'],
                                                  snapshot)
            status, success, stop, _ = self.add_to_path(path, phase_point,
                                                        left, right)
            if stop:
                logger.debug('GROMACS propagation ended at %i. Reason: %s',
                             i, status)
                break
            if i == 0:
                # This step was performed before entering the main loop.
                pass
            elif i > 0:
                out_extnd = self._extend_and_execute_mdrun(tpr_file, cpt_file,
                                                           name)
                out_files.update(out_extnd)
            # Calculate the order parameter using the current system:
            ensemble['system'].particles.set_vel(reverse)
            ensemble['system'].particles.set_pos((conf_abs, None))
            order = self.calculate_order(ensemble)
            # We now have the order parameter, for GROMACS just remove the
            # config file to avoid the GROMACS #conf_abs# backup clutter:
            self._removefile(conf_abs)
            msg_file.flush()
        logger.debug('GROMACS propagation done, obtaining energies')
        msg_file.write('# Propagation done.')
        msg_file.write(f'# Reading energies from: {out_files["edr"]}')
        msg_file.flush()
        energy = self.get_energies(out_files['edr'])
        path.update_energies(energy['kinetic en.'], energy['potential'])
        logger.debug('Removing GROMACS output after propagate.')
        remove = [val for key, val in out_files.items()
                  if key not in ('trr', 'gro', 'g96')]
        self._remove_files(self.exe_dir, remove)
        self._remove_gromacs_backup_files(self.exe_dir)
        return success, status

    def step(self, system, name):
        """Perform a single step with GROMACS.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The system we are integrating.
        name : string
            To name the output files from the GROMACS step.

        Returns
        -------
        out : string
            The name of the output configuration, obtained after
            completing the step.

        """
        initial_conf = self.dump_frame(system)
        # Save as a single snapshot file:
        system.particles.set_pos((initial_conf, None))
        system.particles.set_vel(False)
        out_grompp = self._execute_grompp(self.input_files['input'],
                                          initial_conf,
                                          name)
        out_mdrun = self._execute_mdrun(out_grompp['tpr'], name)
        conf_abs = os.path.join(self.exe_dir, out_mdrun['conf'])
        logger.debug('Obtaining GROMACS energies after single step.')
        energy = self.get_energies(out_mdrun['edr'])
        system.particles.set_pos((conf_abs, None))
        system.particles.set_vel(False)
        system.particles.vpot = energy['potential'][-1]
        system.particles.ekin = energy['kinetic en.'][-1]
        logger.debug('Removing GROMACS output after single step.')
        remove = [val for _, val in out_grompp.items()]
        remove += [val for key, val in out_mdrun.items() if key != 'conf']
        self._remove_files(self.exe_dir, remove)
        return out_mdrun['conf']

    def _prepare_shooting_point(self, input_file):
        """
        Create the initial configuration for a shooting move.

        This creates a new initial configuration with random velocities.
        Here, the random velocities are obtained by running a zero-step
        GROMACS simulation.

        Parameters
        ----------
        input_file : string
            The input configuration to generate velocities for.

        Returns
        -------
        output_file : string
            The name of the file created.
        energy : dict
            The energy terms read from the GROMACS .edr file.

        """
        gen_mdp = os.path.join(self.exe_dir, 'genvel.mdp')
        if os.path.isfile(gen_mdp):
            logger.debug('%s found. Re-using it!', gen_mdp)
        else:
            # Create output file to generate velocities:
            settings = {'gen_vel': 'yes', 'gen_seed': -1, 'nsteps': 0,
                        'continuation': 'no'}
            self._modify_input(self.input_files['input'], gen_mdp, settings,
                               delim='=')
        # Run GROMACS grompp for this input file:
        out_grompp = self._execute_grompp(gen_mdp, input_file, 'genvel')
        remove = [val for _, val in out_grompp.items()]
        # Run GROMACS mdrun for this tpr file:
        out_mdrun = self._execute_mdrun(out_grompp['tpr'], 'genvel')
        remove += [val for key, val in out_mdrun.items() if key != 'conf']
        confout = os.path.join(self.exe_dir, out_mdrun['conf'])
        energy = self.get_energies(out_mdrun['edr'])
        # Remove run-files:
        logger.debug('Removing GROMACS output after velocity generation.')
        self._remove_files(self.exe_dir, remove)
        return confout, energy

    def _read_configuration(self, filename):
        """Read output from GROMACS .g96/gro files.

        Parameters
        ----------
        filename : string
            The file to read the configuration from.

        Returns
        -------
        box : numpy.array
            The box dimensions.
        xyz : numpy.array
            The positions.
        vel : numpy.array
            The velocities.

        """
        box = None
        if self.ext == 'g96':
            _, xyz, vel, box = read_gromos96_file(filename)
        elif self.ext == 'gro':
            _, xyz, vel, box = read_gromacs_gro_file(filename)
        else:
            msg = 'GROMACS engine does not support reading "%s"'
            logger.error(msg, self.ext)
            raise ValueError(msg % self.ext)
        return box, xyz, vel

    def _reverse_velocities(self, filename, outfile):
        """Reverse velocity in a given snapshot.

        Parameters
        ----------
        filename : string
            The configuration to reverse velocities in.
        outfile : string
            The output file for storing the configuration with
            reversed velocities.

        """
        if self.ext == 'g96':
            txt, xyz, vel, _ = read_gromos96_file(filename)
            write_gromos96_file(outfile, txt, xyz, -1 * vel)
        elif self.ext == 'gro':
            txt, xyz, vel, _ = read_gromacs_gro_file(filename)
            write_gromacs_gro_file(outfile, txt, xyz, -1 * vel)
        else:
            msg = 'GROMACS engine does not support writing "%s"'
            logger.error(msg, self.ext)
            raise ValueError(msg % self.ext)

    def modify_velocities(self, ensemble, vel_settings):
        """Modify the velocities of the current state.

        This method will modify the velocities of a time slice.

        Parameters
        ----------
        ensemble : dict
            It contains:

            * `system`: object like :py:class:`.System`
              This is the system that contains the particles we are
              investigating.

        vel_settings: dict
            It contains:

            * `sigma_v`: numpy.array, optional
              These values can be used to set a standard deviation (one
              for each particle) for the generated velocities.
            * `aimless`: boolean, optional
              Determines if we should do aimless shooting or not.
            * `momentum`: boolean, optional
              If True, we reset the linear momentum to zero after
              generating.
            * `rescale or rescale_energy`: float, optional
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
        system = ensemble['system']
        rescale = vel_settings.get('rescale_energy',
                                   vel_settings.get('rescale'))
        if rescale is not None and rescale is not False and rescale > 0:
            msgtxt = 'GROMACS engine does not support energy re-scale.'
            logger.error(msgtxt)
            raise NotImplementedError(msgtxt)
        kin_old = system.particles.ekin
        if vel_settings.get('aimless', False):
            pos = self.dump_frame(system)
            posvel, energy = self._prepare_shooting_point(pos)
            kin_new = energy['kinetic en.'][-1]
            system.particles.set_pos((posvel, None))
            system.particles.set_vel(False)
            system.particles.ekin = kin_new
            system.particles.vpot = energy['potential'][-1]
        else:  # Soft velocity change, from a Gaussian distribution:
            msgtxt = 'GROMACS engine only support aimless shooting!'
            logger.error(msgtxt)
            raise NotImplementedError(msgtxt)
        if vel_settings.get('momentum', False):
            pass
        if kin_old is None or kin_new is None:
            dek = float('inf')
            logger.debug(('Kinetic energy not found for previous point.'
                          '\n(This happens when the initial configuration '
                          'does not contain energies.)'))
        else:
            dek = kin_new - kin_old
        return dek, kin_new

    def integrate(self, ensemble, steps, thermo='full'):
        """
        Perform several integration steps.

        This method will perform several integration steps using
        GROMACS. It will also calculate order parameter(s) and energy
        terms if requested.

        Parameters
        ----------
        ensemble: dict
            It contains:

            * `system`: object like :py:class:`.System`
              The system object gives the initial state for the
              integration. The initial state is stored and the system is
              reset to the initial state when the integration is done.
            * `order_function`: object like :py:class:`.OrderParameter`
              The object used for calculating the order parameter.

        steps : integer
            The number of steps we are going to perform. Note that we
            do not integrate on the first step (e.g. step 0) but we do
            obtain the other properties. This is to output the starting
            configuration.
        thermo : string, optional
            Select the thermodynamic properties we are to calculate.

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
        out_files = {}
        conf_abs = None
        self.energy_terms = self.select_energy_terms(thermo)
        # For step zero, obtain the order parameter:
        if order_function:
            order = self.calculate_order(ensemble)
        else:
            order = None

        for i in range(steps):
            if i == 0:
                out_files = self._execute_grompp_and_mdrun(
                    initial_file,
                    'pyretis-gmx'
                )
                conf_abs = os.path.join(self.exe_dir, out_files['conf'])
            elif 0 < i < steps - 1:
                out_extnd = self._extend_and_execute_mdrun(
                    out_files['tpr'],
                    out_files['cpt'],
                    'pyretis-gmx'
                )
                out_files.update(out_extnd)
            else:
                pass
            # Update with results from previous step:
            results = {}
            if order:
                results['order'] = order
            # Update for order parameter:
            if order_function:
                system.particles.set_pos((conf_abs, None, None))
                order = self.calculate_order(ensemble)
            # Obtain latest energies:
            time1 = i * self.timestep * self.subcycles
            time2 = (i + 1) * self.timestep * self.subcycles
            # time1 and time2 should be correct now, but we are victims
            # of floating points. Subtract/add something small so that
            # we round to correct time.
            time1 -= self.timestep * 0.1
            time2 += self.timestep * 0.1
            energy = self.get_energies(out_files['edr'], begin=time1,
                                       end=time2)
            # Rename energies into the PyRETIS convention:
            results['thermo'] = self.rename_energies(energy)
            yield results

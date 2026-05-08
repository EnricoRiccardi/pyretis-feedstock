# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Module defining the base classes for the PyRETIS output.

Important classes defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PathStorage (:py:class:`.PathStorage`)
    A class for handling storage of external trajectories to an archive.
    In this case, it is simply a folder based structure.

PathStorageTar (:py:class:`.PathStorageTar`)
    A class for handling storage of external trajectories to an archive
    format (tar in this case).

"""
import logging
import tempfile
import shutil
import tarfile
import os
from pyretis.inout.common import OutputBase, create_backup, make_dirs
from pyretis.inout.formats import (
    EnergyPathFormatter,
    OrderPathFormatter,
    PathExtFormatter,
)


logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


def add_to_tar_file(tar_file, files, file_mode='a'):
    """Open a tar file and add some files.

    Parameters
    ----------
    tar_file : string
        The path to the tar file to create.
    files : list of tuples
        The path(s) to the file(s) to add with names to
        create internaly in the tar file. That is:
        ``file[i] = (source[i], target[i])`` where ``target[i]`` will be
        used as the archname for ``source[i]``.
    file_mode : string
        Determines how we try to open the tar file.

    Returns
    -------
    out : boolean
        True if we managed to write to the file.

    """
    try:
        with tarfile.open(tar_file, file_mode) as tar:
            for src, trg in files:
                tar.add(src, arcname=trg)
        return True
    except tarfile.ReadError:
        logger.warning('Could not open tar file: "%s"', tar_file)
        return False


def generate_traj_names(path, target_dir):
    """Generate new file names for moving copying paths.

    Parameters
    ----------
    path : object like :py:class:`.PathBase`
        This is the path object we are going to store.
    target_dir : string
        The location where we are moving the path to.

    Returns
    -------
    files : list of tuples of strings
        The files to add and the file names to use internally
        for storing. That is: ``file[i] = (source[i], target[i])``
        where ``target[i]`` is the internal storage name.

    """
    source = set()
    files = []
    for phasepoint in path.phasepoints:
        pos_file, _ = phasepoint.particles.get_pos()
        if pos_file not in source:
            source.add(pos_file)
            dest = os.path.join(target_dir, os.path.basename(pos_file))
            files.append((pos_file, dest))
    return files


class PathStorage(OutputBase):
    """A class for handling storage of external trajectories.

    Attributes
    ----------
    target : string
        Determines the target for this output class. Here it will
        be a file archive (i.e. a directory based collection of
        files).
    archive_acc : string
        Basename for the archive with accepted trajectories.
    archive_rej : string
        Basename for the archive with rejected trajectories.
    archive_traj : string
        Basename for a sub-folder containing the actual files
        for a trajectory.
    formatters : dict
        This dict contains the formatters for writing path data,
        with default filenames used for them.
    out_dir_fmt : string
        A format to use for creating directories within the archive.
        This one is applied to the step number for the output.

    """

    target = 'file-archive'
    archive_acc = 'traj-acc'
    archive_rej = 'traj-rej'
    archive_traj = 'traj'
    formatters = {
        'order': {'fmt': OrderPathFormatter(), 'file': 'order.txt'},
        'energy': {'fmt': EnergyPathFormatter(), 'file': 'energy.txt'},
        'traj': {'fmt': PathExtFormatter(), 'file': 'traj.txt'},
    }
    out_dir_fmt = '{}'

    def __init__(self):
        """Set up the storage.

        Note that we here do not pass any formatters to the parent
        class. This is because this class is less flexible and
        only intended to do one thing - write path data for external
        trajectories.

        """
        super().__init__(None)

    def archive_name_from_status(self, status):
        """Return the name of the archive to use."""
        return self.archive_acc if status == 'ACC' else self.archive_rej

    def output_path_files(self, step, data, target_dir):
        """Write the output files for energy, path and order parameter.

        Parameters
        ----------
        step : integer
            The current simulation step.
        data : list
            Here, ``data[0]`` is assumed to be an object like
            :py:class:`.Path` and `data[1]`` a string containing the
            status of this path.
        target_dir : string
            The path to where we archive the files.

        Returns
        -------
        files : list of tuple of strings
            These are the files created. The tuple contains the files as
            a full path and a relative path (to the given
            target directory). The form is
            ``files[i] = (full_path[i], relative_path[i])``.
            The relative path is useful for organizing internally in
            archives.

        """
        path, status, = data[0], data[1]
        files = []
        for key, val in self.formatters.items():
            logger.debug('Storing: %s', key)
            fmt = val['fmt']
            full_path = os.path.join(target_dir, val['file'])
            relative_path = os.path.join(
                self.out_dir_fmt.format(step), val['file']
            )
            files.append((full_path, relative_path))
            with open(full_path, 'w', encoding="utf8") as output:
                for line in fmt.format(step, (path, status)):
                    output.write(f'{line}\n')
        return files

    def output(self, step, data):
        """Format the path data and store the path.

        Parameters
        ----------
        step : integer
            The current simulation step.
        data : list
            Here, ``data[0]`` is assumed to be an object like
            :py:class:`.Path`, ``data[1]`` a string containing the
            status of this path and ``data[2]`` the path ensemble for
            which the path was generated.

        Returns
        -------
        files : list of tuples of strings
            The files added to the archive.

        """
        path, status, path_ensemble = data
        archive = self.archive_name_from_status(status)
        # This is the path on form: /path/to/000/traj/11
        archive_path = os.path.join(
            path_ensemble.directory['traj'],
            archive,
            self.out_dir_fmt.format(step),
        )
        # To organize things we create a subfolder for storing the
        # files. This is on form: /path/to/000/traj/11/traj
        traj_dir = os.path.join(archive_path, self.archive_traj)
        # Create the needed directories:
        make_dirs(traj_dir)
        # Write order, energy and traj files to the archive:
        files = self.output_path_files(step, data, archive_path)
        # Get name of files we need to copy for the current path:
        traj_files = generate_traj_names(
            path,
            os.path.join(self.out_dir_fmt.format(step), self.archive_traj)
        )
        # Copy the actual trajectory files:
        for (file_src, file_dest) in traj_files:
            target = os.path.join(traj_dir, os.path.basename(file_dest))
            shutil.copyfile(file_src, target)
            files.append((target, file_dest))
        # Copy the ensemble.restart file if path was accepted:
        if data[1] == 'ACC':
            file_src = os.path.join(path_ensemble.directory['path_ensemble'],
                                    'ensemble.restart')
            if os.path.isfile(file_src):
                target = os.path.join(archive_path, 'ensemble.restart')
                shutil.copyfile(file_src, target)
                file_dest = os.path.join(self.out_dir_fmt.format(step),
                                         'ensemble.restart')
                files.append((target, file_dest))

        return files

    def write(self, towrite, end='\n'):
        """We do not need the write method for this object."""
        logger.critical(
            '%s does *not* support the "write" method!',
            self.__class__.__name__
        )

    def formatter_info(self):
        """Return info about the formatters."""
        return [val['fmt'].__class__ for val in self.formatters.values()]

    def __str__(self):
        """Return basic info."""
        return f'{self.__class__.__name__} - archive writer.'


class PathStorageTar(PathStorage):
    """A class for handling storage of external trajectories.

    Here, we collect the external trajectories into tar archives.

    """

    archive_acc = 'traj-acc.tar'
    archive_rej = 'traj-rej.tar'

    @staticmethod
    def _add_to_archive(archive_file, files):
        """Add some files to an archive file.

        Parameters
        ----------
        archive_file : string
            The path to the archive file to create.
        files : list of tuples
            The path(s) to the file(s) to add with names to
            create internally in the archive. That is, the form is:
            ``file[i] = (source[i], target[i])`` where ``target[i]``
            will be used as the archive name for file ``source[i]``.

        Returns
        -------
        add : boolean
            True if we added files to the archive.

        """
        add = add_to_tar_file(archive_file, files, file_mode='a')
        if not add:
            logger.info(
                ('Could not write the trajectory file.'
                 ' Will try recreating it.')
            )
            logger.info('Backing up old trajectory file first.')
            logger.info(create_backup(archive_file))
            add = add_to_tar_file(archive_file, files, file_mode='w')
        return add

    def output(self, step, data):
        """Format the path data and store the path.

        Parameters
        ----------
        step : integer
            The current simulation step.
        data : list
            Here, ``data[0]`` is assumed to be an object like
            :py:class:`.Path`, ``data[1]`` a string containing the
            status of this path and ``data[2]`` the path ensemble for
            which the path was generated.

        Returns
        -------
        files : list of tuples of strings
            The files added to the archive.

        """
        path, status, path_ensemble = data
        tar_file = self.archive_name_from_status(status)
        tar_path = os.path.join(path_ensemble.directory['traj'], tar_file)
        files = []
        # We create the files in a temporary directory and then add them
        # to the tar archive:
        with tempfile.TemporaryDirectory() as tempdir:
            # Write order, energy and traj files:
            files1 = self.output_path_files(step, data, tempdir)
            if self._add_to_archive(tar_path, files1):
                files.extend(files1)
        # Add the trajectory files to the archive:
        target_dir = os.path.join(self.out_dir_fmt.format(step),
                                  self.archive_traj)
        files2 = generate_traj_names(path, target_dir)
        if self._add_to_archive(tar_path, files2):
            files.extend(files2)
        # Include the ensemble.restart file if path was accepted:
        if data[1] == 'ACC':
            file_src = os.path.join(path_ensemble.directory['path_ensemble'],
                                    'ensemble.restart')
            if os.path.isfile(file_src):
                file_dest = os.path.join(self.out_dir_fmt.format(step),
                                         'ensemble.restart')
                files.append((file_src, file_dest))
        return files

    def __str__(self):
        """Return basic info."""
        return f'{self.__class__.__name__} - TAR archive writer.'

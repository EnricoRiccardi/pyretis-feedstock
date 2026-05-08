# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Definition of a class for handling output related to simulations.

Important classes defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Task (:py:class:`.Task`)
    Base class for tasks. This is used by :py:class:`.SimulationTask`
    and :py:class:`.OutputTask`.

OutputTask (:py:class:`.OutputTask`)
    A class representing a simulation output task.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

get_task_type (:py:meth:`.get_task_type`)
    Do additional handling for a path task.

get_file_mode (:py:meth:`.get_file_mode`)
    Determine if we should append or backup existing files.

task_from_settings (:py:meth:`.task_from_settings`)
    Create output task from simulation settings.

Important variables defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

OUTPUT_TASKS (:py:data:`.OUTPUT_TASKS`)
    A dictionary defining the different output tasks known to PyRETIS.

"""
import logging
from pyretis.core.common import initiate_instance
from pyretis.inout.common import generate_file_name
from pyretis.inout.formats import (
    EnergyFormatter,
    OrderFormatter,
    CrossFormatter,
    SnapshotFormatter,
    ThermoTableFormatter,
    PathTableFormatter,
    EnergyPathFormatter,
    OrderPathFormatter,
    PathIntFormatter,
    PathEnsembleFormatter,
    RETISResultFormatter,
)
from pyretis.inout.fileio import FileIO
from pyretis.inout.screen import ScreenOutput
from pyretis.inout.archive import PathStorage
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


OUTPUT_TASKS = {}
"""Define a set of known output tasks.

The output tasks are defined as dictionaries with the following keys:

* target : string
    "file" or "screen", defines where the task writes to.
* filename : string
    A default file name for an output file if writing to a file.
* result : tuple of strings
    Determines what item from the result dictionary we are outputting.
* when : string
    Determines what input setting from the "output" section is used to
    define the output frequency. Default values are defined by the
    output section, see: py:mod:`pyretis.inout.settings.settings`.
* formatter : object like :py:class:`.OutputFormatter`
    Selects the formatter for the output.
* writer : object like :py:class:`.OutputBase`
    Selects the writer for the output.
* settings : tuple of strings, optional
    A dict with additional settings which can be passed to the
    formatter if needed. These settings can, for instance, be
    defined in the output section of the input file.

The writer can be defined explicitly, or via the formatter. If a
formatter is given, then the generic :py:class:`.FileIO` will be used.
If no formatter is given, then the writer is assumed to be given.

"""
OUTPUT_TASKS['energy'] = {
    'target': 'file',
    'filename': 'energy.txt',
    'result': ('thermo',),
    'when': 'energy-file',
    'formatter': EnergyFormatter,
}
OUTPUT_TASKS['order'] = {
    'target': 'file',
    'filename': 'order.txt',
    'result': ('order',),
    'when': 'order-file',
    'formatter': OrderFormatter,
}
OUTPUT_TASKS['cross'] = {
    'target': 'file',
    'filename': 'cross.txt',
    'result': ('cross',),
    'when': 'cross-file',
    'formatter': CrossFormatter,
}
OUTPUT_TASKS['traj-txt'] = {
    'target': 'file',
    'filename': 'traj.txt',
    'result': ('system',),
    'when': 'trajectory-file',
    'formatter': SnapshotFormatter,
}
OUTPUT_TASKS['traj-xyz'] = {
    'target': 'file',
    'filename': 'traj.xyz',
    'result': ('system',),
    'when': 'trajectory-file',
    'formatter': SnapshotFormatter,
}
OUTPUT_TASKS['thermo-screen'] = {
    'target': 'screen',
    'result': ('thermo',),
    'when': 'screen',
    'formatter': ThermoTableFormatter,
}
OUTPUT_TASKS['thermo-file'] = {
    'target': 'file',
    'filename': 'thermo.txt',
    'result': ('thermo',),
    'when': 'energy-file',
    'formatter': ThermoTableFormatter,
}
OUTPUT_TASKS['pathensemble'] = {
    'target': 'file',
    'filename': 'pathensemble.txt',
    'result': ('pathensemble',),
    'when': 'pathensemble-file',
    'formatter': PathEnsembleFormatter,
}
OUTPUT_TASKS['pathensemble-screen'] = {
    'target': 'screen',
    'result': ('pathensemble',),
    'when': 'screen',
    'formatter': PathTableFormatter,
}
OUTPUT_TASKS['pathensemble-retis-screen'] = {
    'target': 'screen',
    'result': ('pathensemble',),
    'when': 'screen',
    'formatter': RETISResultFormatter,
}
OUTPUT_TASKS['path-order'] = {
    'target': 'file',
    'filename': 'order.txt',
    'result': ('path', 'status'),
    'when': 'order-file',
    'formatter': OrderPathFormatter,
}
OUTPUT_TASKS['path-energy'] = {
    'target': 'file',
    'filename': 'energy.txt',
    'result': ('path', 'status'),
    'when': 'energy-file',
    'formatter': EnergyPathFormatter,
}
OUTPUT_TASKS['path-traj-int'] = {
    'target': 'file',
    'filename': 'traj.txt',
    'result': ('path', 'status'),
    'when': 'trajectory-file',
    'formatter': PathIntFormatter,
}
OUTPUT_TASKS['path-traj-ext'] = {
    'target': 'file-archive',
    'filename': 'traj.txt',
    'result': ('path', 'status', 'pathensemble'),
    'when': 'trajectory-file',
    'writer': PathStorage,
}


class Task:
    """Base representation of a "task".

    A task is just something that is supposed to be executed at
    a certain point. This class will just set up functionality
    that is common for output tasks and for simulation tasks.

    Attributes
    ----------
    when : dict
        Determines when the task should be executed.

    """

    _ALLOWED_WHEN = {'every', 'at'}

    def __init__(self, when):
        """Initialise the task.

        Parameters
        ----------
        when : dict, optional
            Determines if the task should be executed.

        """
        self._when = None
        self.when = when

    @property
    def when(self):
        """Return the "when" property."""
        return self._when

    @when.setter
    def when(self, when):
        """Update self to new value(s).

        It will only update `self.when` for the keys given in the
        input `when`.

        Parameters
        ----------
        when : dict
            This dict contains the settings to update.

        Returns
        -------
        out : None
            Returns `None` but modifies `self.when`.

        """
        if when is None:
            self._when = when
        else:
            if self._when is None:
                self._when = {}
            for key, val in when.items():
                if key in self._ALLOWED_WHEN:
                    self._when[key] = val
                else:
                    logger.warning(
                        'Ignoring unknown "when" setting: "%s"', key
                    )

    def execute_now(self, step):
        """Determine if a task should be executed.

        Parameters
        ----------
        step : dict of ints
            Keys are 'step' (current cycle number), 'start' cycle number at
            start 'stepno' the number of cycles we have performed so far.

        Returns
        -------
        out : boolean
            True of the task should be executed.

        """
        if self.when is None:
            return True
        exe = False
        if 'every' in self.when:
            exe = step['stepno'] % self.when['every'] == 0
        if not exe and 'at' in self.when:
            try:
                exe = step['step'] in self.when['at']
            except TypeError:
                exe = step['step'] == self.when['at']
        return exe

    def task_dict(self):
        """Return basic info about the task."""
        return {'when': self.when}


class OutputTask(Task):
    """A base class for simulation output.

    This class will handle an output task for a simulation. The
    output task consists of one object which is responsible for
    formatting the output data and one object which is responsible
    for writing that data, for instance to the screen or to a file.

    Attributes
    ----------
    target : string
        This string identifies what kind of output we are dealing with.
        This will typically be either "screen" or "file".
    name : string
        This string identifies the task, it can, for instance, be used
        to reference the dictionary used to create the writer.
    result : tuple of strings
        This string defines the result we are going to output.
    writer : object like :py:class:`.OutputBase`
        This object will handle the actual outputting
        of the result.
    when : dict
        Determines if the task should be executed.

    """

    def __init__(self, name, result, writer, when):
        """Initialise the generic output task.

        Parameters
        ----------
        name : string
            This string identifies the task, it can, for instance, be used
            to reference the dictionary used to create the writer.
        result : list of strings
            These strings define the results we are going to output.
        writer : object like :py:class:`.IOBase`
            This object will handle formatting of the actual result and
            output to screen or to a file.
        when : dict
            Determines when the output should be written. Example:
            `{'every': 10}` will be executed at every 10th step.

        """
        super().__init__(when)
        self.target = writer.target
        self.name = name
        self.result = result
        self.writer = writer

    def output(self, simulation_result):
        """Output given results from simulation steps.

        This will output the task using the result found in the
        `simulation_result` which should be the dictionary returned
        from a simulation object (e.g. object like
        :py:class:`.Simulation`) after a step.

        Parameters
        ----------
        simulation_result : dict
            This is the result from a simulation step.

        Returns
        -------
        out : boolean
            True if the writer wrote something, False otherwise.

        """
        step = simulation_result['cycle']
        if not self.execute_now(step):
            return False
        result = []
        for res in self.result:
            if res not in simulation_result:
                return False  # Requested result was not ready at this step.
            result.append(simulation_result[res])
        if len(result) == 1:
            return self.writer.output(step['step'], result[0])
        return self.writer.output(step['step'], result)

    def __str__(self):
        """Print information about the output task."""
        msg = [f'Output task {self.name}']
        msg += [f'- Output frequency: {self.when}']
        msg += [f'- Acting on result(s): {self.result}']
        msg += [f'- Target is: {self.target}']
        msg += [f'- Writer is:\xa0{self.writer}']
        return '\n'.join(msg)

    def task_dict(self):
        """Return a dict with info about the task."""
        return {'name': self.name,
                'when': self.when,
                'result': self.result,
                'target': self.target,
                'writer': self.writer.__class__,
                'formatter': self.writer.formatter_info()}


def get_task_type(task, engine):
    """Do additional handling for a path task.

    The path task is special since we do very different things for
    external paths. The set-up required to do this is handled here.

    Parameters
    ----------
    task : dict
        Settings related to the specific task.
    engine : object like :py:class:`.EngineBase`
        This object is used to determine if we need to do something
        special for external engines. If no engine is given, we do
        not do anything special.

    Returns
    -------
    out : string
        The task type we are going to be creating for.

    """
    if task['type'] == 'path-traj-{}':
        if engine is None or engine.engine_type in ('internal', 'openmm'):
            fmt = 'int'
        else:
            fmt = 'ext'
        return task['type'].format(fmt)
    return task['type']


def get_file_mode(settings):
    """Determine if we should append or backup existing files.

    This method translates the backup settings into a file mode string.
    We assume here that the file is opened for writing.

    Parameters
    ----------
    settings : dict
        The simulation settings.

    Returns
    -------
    file_mode : string
        A string representing the file mode to use.

    """
    file_mode = 'w'
    try:
        old = settings['output']['backup'].lower()
        if old == 'append':
            logger.debug('Will append to existing files.')
            file_mode = 'a'
    except AttributeError:
        logger.warning('Could not understand setting for "backup"'
                       ' in "output" section.')
        old = 'backup'
        logger.warning('Handling of existing files is set to: "%s"', old)
        settings['output']['backup'] = old
    return file_mode


def task_from_settings(task, settings, directory, engine, progress=False):
    """Create output task from simulation settings.

    Parameters
    ----------
    task : dict
        Settings for creating a task. This dict contains the type and
        name of the task to create. It can also contain overrides to
        the default settings in :py:data:`OUTPUT_TASKS`.
    settings : dict
        Settings for the simulation.
    directory : string
        The directory to write output files to.
    engine : object like :py:class:`.EngineBase`
        This object is used to determine if we need to do something
        special for external engines. If no engine is given, we do
        not do anything special.
    progress : boolean, optional
        For some simulations, the user may select to display a
        progress bar. We will then just disable the other screen
        output.

    Returns
    -------
    out : object like :py:class:`.OutputTask`
        An output task we can use in the simulation.

    """
    task_type = get_task_type(task, engine)
    task_settings = OUTPUT_TASKS[task_type].copy()
    # Override defaults if any:
    for key in task_settings:
        if key in task:
            task_settings[key] = task[key]

    when = {'every': settings['output'][task_settings['when']]}
    if when['every'] < 1:
        logger.info('Skipping output task %s (freq < 1)', task_type)
        return None

    target = task_settings['target']
    if target == 'screen' and progress:
        logger.info(
            'Disabling output to screen %s since progress bar is ON',
            task['name'],
        )
        return None
    formatter = None
    # Initiate the formatter, note that here we can customize the formatter
    # by supplying arguments to it. This was supported in a previous version
    # of PyRETIS, but for now, none of the formatters needs settings to be
    # created.
    if task_settings.get('formatter', None) is not None:
        formatter = initiate_instance(
            task_settings['formatter'],
            task_settings.get('formatter-settings', {}),
        )
    # Create writer:
    writer = None
    if target == 'screen':
        klass = task_settings.get('writer', ScreenOutput)
        writer = klass(formatter)
    if target in ('file', 'file-archive'):
        filename = generate_file_name(task_settings['filename'], directory,
                                      settings)
        file_mode = get_file_mode(settings)
        if target == 'file':
            klass = task_settings.get('writer', FileIO)
            writer = klass(filename, file_mode, formatter, backup=True)
        if target == 'file-archive':
            klass = task_settings.get('writer', PathStorage)
            writer = klass()
    # Finally make the output task:
    if writer is not None:
        return OutputTask(task['name'], task_settings['result'], writer, when)
    logger.warning('Unknown target "%s". Ignoring task: %s',
                   target, task_type)
    return None

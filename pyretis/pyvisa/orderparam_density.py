# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Compiler of PyRETIS simulation data.

This module is part of the PyRETIS library and can be used both for compiling
the simulation data into a compressed file and/or load the data for later
visualization.


Important classes defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Trajectory (:py:class:`.PathBase`)
    A base class to store trajectories composed by only orderp,
    collective variables and energies.

PathDensity (:py:class:`.PathBase`)
    A base class to assemble the data.

PathVisualize (:py:class:`.PathVisualize`)
    A base class to prepare for the visualization.


Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pyvisa_zip (:py:func: `.pyvisa_zip`)
    Compress the PyVisA output file in a .zip format.

pyvisa_unzip (:py:func: `.pyvisa_unzip`)
    Decompress the zip into PyVisA output.

remove_nan (:py:func: `.remove_nan`)
    Checks for the presence of nan values and replace them with a local,
    if available.

pyvisa_compress (:py:func: `.pyvisa_compress`)
    Compress PyRETIS outputs to a .hdf5 file.

"""
import os
import timeit
import warnings
import zipfile
import logging

import numpy as np
import pandas as pd
from scipy import stats  # pylint: disable=import-error

from pyretis.inout.settings import parse_settings_file
from pyretis.pyvisa.common import where_from_to, get_cv_names

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Hard-coded labels for energies and time/cycle steps
ENERGYLABELS = ['time', 'cycle', 'potE', 'kinE', 'totE']


def pyvisa_zip(input_file):
    """Zip compress file of simulation data.

    Parameters
    ----------
    input_file : string
        The file to compress.

    """
    with zipfile.ZipFile(input_file + '.zip', 'w') as zipped_file:
        zipped_file.write(input_file, compress_type=zipfile.ZIP_DEFLATED)

    os.remove(input_file)


def pyvisa_unzip(origin, destination=None):
    """Unzip compressed file before load in visualizer.

    Parameters
    ----------
    origin : string
        Zipped file to unzip.
    destination : string, optional
        Unzipped file name.

    """
    msg = '###################################################\n'
    msg += '# File type recognized as `.zip`, unzipping to tmp file \n'
    msg += f'# {destination} before loading.\n'
    msg += '###################################################\n'
    logger.info(msg)
    with zipfile.ZipFile(origin) as zipped:
        zipped.extractall(path=os.path.dirname(os.path.abspath(origin)))
        # We here assume that only a file is stored in 'zipped'.
        if destination is not None:
            os.rename(zipped.namelist()[0], destination)


def remove_nan(data):
    """Remove nan from data.

    The function shall remove initial nan, assuming that they are originated
    by incomplete initial conditions (e.g. no energy file). In the case that
    nan appears as last cycle, it will not be fixed and an error shall rise
    up later in the code.

    Parameters
    ----------
    data : list, dict, object like pandas.DataFrame or pandas.Series
        Input data structure. If `np.nan` are present, they are
        replaced by the following entry. The method accounts for
        multiple consecutive `np.nan` occurrences.

    """
    nan = True
    inan = -1
    while nan:
        nan = False
        if type(data) in (dict, pd.DataFrame):
            for keys in data:
                if isinstance(data, pd.DataFrame):
                    # Pandas 2.0+ Copy-on-Write: data[keys] returns a copy,
                    # so we must modify a copy and assign it back explicitly.
                    col = data[keys].copy()
                    remove_nan(col)
                    data[keys] = col
                else:
                    remove_nan(data[keys])
            break
        for idx, data_point in reversed(list(enumerate(data))):
            if type(data_point) in (list, np.ndarray):
                remove_nan(data_point)
            else:
                # Check if data_point is a numeric type before calling np.isnan
                if (isinstance(data_point,
                               (int, float, np.floating, np.integer))
                        and np.isnan(data_point)):
                    nan = True
                    inan = idx
                    break
        if nan and inan == len(data) - 1:
            nan = False
        if nan and isinstance(data, pd.Series):
            data.iloc[inan] = data.iloc[inan + 1]
        elif nan:
            data[inan] = data[inan + 1]


def pyvisa_compress(runpath, input_file, pyvisa_dict):
    """Compress simulation data.

    Parameters
    ----------
    runpath : string
        The execution folder where the input files are.
    input_file : string
        The input file for compression.
    pyvisa_dict : dict
        It determines the section of pyvisa to use.

    """
    if not input_file:
        raise ValueError('No input file')

    if input_file.endswith(('hdf5', 'zip')):
        msg = 'Cannot compress an already compressed file.'
        raise ImportError(msg)
    p_data = PathDensity(runpath, input_file)
    p_data.walk_dirs(only_ops=pyvisa_dict['only_order'])
    p_data.hdf5_data()


class Trajectory:  # pylint: disable=too-few-public-methods
    """Class representing a simulation trajectory.

    This class defines the trajectories from the completed
    simulations, with all information available. The labels of
    the order parameter and collective variables will either be
    labeled in the fashion of opx, or from the names in the input file
    if the number of names given and number of descriptors in the system
    are equal.
    """

    def __init__(self, frames, info):
        """Initialize the class.

        Parameters
        ----------
        frames : pandas Dataframe or dict
            Dataframe/dict containing order parameter and energy-data
            for the trajectory. It contains:

            * `OP1`: string
              Order parameter.
            * `OP2`: string
              Collective variables, all other CV's will be named in
              increasing fashion, OP3, OP4 etc.
            * `PotE`: string
              Potential energy.
            * `KinE`: string
              Kinetic energy.
            * `TotE`: string
              Total energy.

        info : dict
            Dictionary containing information about the trajectory.
            It contains:

            * `ensemble_name`: string
              The name of the ensemble.
            * `cycle`: integer
              The cycle number.
            * `status`: string
              Status of accepted/rejected etc.
            * `MC-move`: string
              Generation move of the trajectory.
            * `MC-start`: string
              Generation starting point of the trajectory.
            * `ordermax`: float
              Max value of OP.
            * `ordermin`: float
              Min value of OP.
            * `stored`: boolean
              True if the trajectory has existing trajectory files.

        """
        self.info = info
        self.frames = frames


class PathDensity:
    """Perform the path density analysis.

    This class defines the path density analysis for completed simulations
    with several order parameters.

    """

    def __init__(self, basepath='.', iofile=None):
        """Initialize the class.

        Parameters
        ----------
        basepath : string, optional
            The path to the input file.
        iofile : string, optional
            The input file.

        Attributes
        ----------
        traj_dict : dict
            Values of order params and energy in all ensembles and
            info about trajectories.
            To each key, ensemble_name (e.g. 000, 001, etc.)
            the value is the list of respective Trajectory objects.

        infos : dict
            Information about the simulation,
            it contains:

            * `ensemble_names`: list
              List of ensemble names.
            * `interfaces`: list
              List of interface positions.
            * `num_op`: integer
              Number of order parameters.
            * `op_labels`: list
              List of order parameter names.
            * `energy_labels`: list
              List of energy entry labels.

        """
        self.basepath = os.path.abspath(basepath)

        if isinstance(iofile, PathDensity):
            self.traj_dict = iofile.traj_dict
            self.infos = iofile.infos
            return

        self.iofile = os.path.join(self.basepath, iofile)

        if iofile.endswith('.hdf5'):
            self.initialize_compressed(iofile)
            return

        # Getting interfaces from iofile
        settings = parse_settings_file(self.iofile)

        interfaces = settings.get('simulation', {}).get('interfaces', [])
        if not interfaces:
            raise ValueError('Input file does not contain interface setting')

        intnames = ['0$^{-}$', '0$^{+}$']
        for i in range(1, len(interfaces) - 1):
            intnames.append(str(i) + '$^{+}$')
        path = []
        # Getting ensembles/folders from directory
        for fol in sorted(os.listdir(basepath)):
            if os.path.isdir(os.path.join(basepath, fol)) and fol.isdigit():
                path.append(fol)

        if not path:
            raise ValueError('No files to analyse.')

        # Getting order parameters from order-file of first folder in path.
        # Careful here because the order parameter function and/or the
        # collective variables ones can produce a LIST. In this case the
        # counting will not correspond.
        filename = os.path.join(basepath, path[0], 'order.txt')
        with open(filename, 'r', encoding='utf-8') as figa:
            lines = len(figa.readlines())

        if lines <= 3:
            raise ValueError(f'File order.txt in {path[0]} empty.')

        # Crash-proof, num_op is the mode of length of the last lines.
        with open(os.path.join(basepath, path[0], 'order.txt'),
                  encoding='utf-8') as temp:
            lasts = temp.readlines()[-10:]
            num_op = int(stats.mode([len(i.split()) for i in lasts])[0]) - 1

        op_names = get_cv_names(settings)
        op_labels = op_names if len(op_names) == num_op else \
            [f'op{i}' for i in range(1, num_op + 1)]

        self.traj_dict = {}
        self.order_parameter = op_labels[0]
        self.infos = {'ensemble_names': path,
                      'interfaces': interfaces,
                      'intf_names': intnames,
                      'num_op': num_op,
                      'op_labels': op_labels,
                      'energy_labels': ENERGYLABELS}

    def walk_dirs(self, only_ops=False):
        """Create a lists in acc or rej dictionary for all order parameters.

        First generate list of folders/ensembles to iterate through.
        Then search for number of order parameters(columns) in file in one
        of the folders of path, and create lists in acc/rej dictionaries
        for all order parameters.

        Lastly iterate through all folders and files, filling in correct data
        to the lists and dictionaries.

        Parameters
        ----------
        only_ops : boolean, optional
            If true, PathDensity will not collect data from energy files.

        """
        tic = [timeit.default_timer(), None]
        logger.info('###################################################')
        logger.info('# PathDensity performing "walk" in \n# %s/',
                    os.getcwd())
        logger.info('# Number of subfolders (0**) = %s',
                    len(self.infos['ensemble_names']))
        logger.info('# Found %s order parameters in output',
                    self.infos['num_op'])
        logger.info('###################################################')

        # Looping over folders, reading energy and orderP
        cmax, cmin = 0, 0
        for ensemble_name in self.infos['ensemble_names']:

            tic[1] = timeit.default_timer()
            logger.info('Reading data from %s', ensemble_name)
            self.get_traj_op(ensemble_name)
            if not only_ops:
                self.get_traj_energy(ensemble_name)

            # Sort traj_dict by cycle
            self.traj_dict[ensemble_name] = \
                dict(sorted(self.traj_dict[ensemble_name].items()))
            for cyc in self.traj_dict[ensemble_name].keys():
                traj = self.traj_dict[ensemble_name][cyc]
                if int(cyc) > cmax:
                    cmax = cyc
                if self.order_parameter in traj.frames.columns:
                    traj.info['ordermax'] = np.max(
                        traj.frames[self.order_parameter])
                    traj.info['ordermin'] = np.min(
                        traj.frames[self.order_parameter])
                    traj.info['length'] = len(traj.frames.index)
            line = ('Done with folder, time used: '
                    f'{timeit.default_timer() - tic[1]:4.4f}s, proceeding.')
            logger.info('=' * len(line) + '\n' + line)

        self.infos['cycles'] = [cmin, cmax]

        logger.info('###################################################')
        logger.info('# Data successfully retrieved, in cycles:')
        logger.info('# %s to %s',
                    self.infos['cycles'][0], self.infos['cycles'][1])
        logger.info('# Time spent: %.2fs',
                    timeit.default_timer() - tic[0])
        logger.info('###################################################')

    def hdf5_data(self):
        """Compress the data to a .hdf5 file."""
        logger.info('###################################################')
        logger.info('# Compress dictionaries to file')

        data = pd.Series(self.traj_dict)
        infos = pd.Series(self.infos)

        pfile = 'pyvisa_compressed_data.hdf5'
        warnings.simplefilter(action='ignore',
                              category=pd.errors.PerformanceWarning)
        data.to_hdf(pfile, key='data', mode='w')
        infos.to_hdf(pfile, key='infos')
        logger.info('# %s', pfile)
        pyvisa_zip(pfile)
        logger.info('# %s.zip', pfile)
        logger.info('###################################################')

    def get_traj_op(self, ensemble_name):
        """Read order.txt files and collects data OP data.

        Parameters
        ----------
        ensemble_name : string

        Updates
        -------
        traj_dict : dict
            Dictionary containing Trajectory objects from simulation
            with filled in OP-data.

        """
        bad_characters = ["(", ")", ",", '"', "'"]
        self.traj_dict[ensemble_name] = {}
        weight = []
        statw = []
        frames = {}
        cycle = -1
        info = {}
        with open(os.path.join(self.basepath, ensemble_name, 'order.txt'),
                  'r', encoding='utf-8') as order:
            for idx, oline in enumerate(order):
                if oline.startswith('#') and 'Cycle' not in oline:
                    continue
                if oline.startswith('#') and 'Cycle' in oline:
                    if idx > 1:
                        self.fill_op(frames, info, ensemble_name, cycle)

                    text = oline.split()
                    cycle = int(text[2].rstrip(','))
                    frames = {}
                    info = {'ensemble_name': ensemble_name,
                            'cycle': cycle,
                            'status': text[4].rstrip(','),
                            'stored': False}

                    # Check for stored trajectory-files
                    if info['status'] == 'ACC':
                        status_line = 'traj-acc'

                        # Append statistical weights
                        statw.append(1)
                        weight.append(1)
                    else:
                        status_line = 'traj-rej'
                        # Append statistical weights
                        if weight:
                            weight[-1] += 1
                        statw.append(0)

                    path_to_cycles = (
                        f"{os.getcwd()}/{ensemble_name}/traj/{status_line}"
                    )

                    info['stored'] = os.path.exists(path_to_cycles +
                                                    f'/{cycle}/traj')

                    if len(text) > 6:  # For backward comparability
                        # Add MC move
                        info['MC-move'] = ''.join(filter(
                            lambda i: i not in bad_characters, text[6]))
                else:
                    # Read data line
                    odata = oline.split()
                    if not odata:
                        continue
                    frames[odata[0]] = {self.order_parameter: np.nan}

                    # Fill in OP-data to the frame
                    for label, data in zip(self.infos['op_labels'], odata[1:]):
                        frames[odata[0]][label] = np.array(data)

            self.fill_op(frames, info, ensemble_name, cycle)

        # Give weights to the trajectories
        for cycle in self.traj_dict[ensemble_name].keys():
            traj = self.traj_dict[ensemble_name][cycle]
            index = list(self.traj_dict[ensemble_name].keys()).index(cycle)
            if statw[index] == 1:
                traj.info['weight'] = weight[0]
                weight.pop(0)
            else:
                traj.info['weight'] = 0

    def fill_op(self, frames, info, ensemble_name, cycle):
        """Fill in OP-data to a trajectory.

        Function that fills the dictionary traj_dict containing
        Trajectory objects with data from the order parameter
        and the collective variables from the simulation.

        Parameters
        ----------
        frames: dict
            Dictionary containing energy-data for the trajectory.
        info : dict
            Information about the trajectory.
        ensemble_name : string
            The name of the ensemble.
        cycle : integer
            Cycle number.

        Updates
        -------
        traj_dict : Trajectory object
            Trajectory object with filled in OP-data.

        """
        orderp = self.order_parameter
        frames = pd.DataFrame.from_dict(
            frames, orient='index', dtype=float)
        for var in frames.columns:
            remove_nan(frames[var])
        if orderp in frames.columns and \
                frames[orderp].isnull().sum() != len(frames[orderp].index):
            start_end = where_from_to(np.array(frames[orderp]),
                                      self.infos['interfaces'][0],
                                      self.infos['interfaces'][-1])

            info['reactive'] = start_end[0] + start_end[1] == 'LR'
            self.traj_dict[ensemble_name][cycle] = Trajectory(frames, info)

    def get_traj_energy(self, ensemble_name):
        """Read energy.txt files and collects energy-data.

        Function that fills the dictionary traj_dict containing
        Trajectory objects with energy data from the simulation.

        Parameters
        ----------
        ensemble_name : string
            The name of the ensemble.

        Updates
        -------
        traj_dict : dict
            Dictionary containing Trajectory objects from simulation
            with filled in energy-data.

        """
        en_file = os.path.join(self.basepath, ensemble_name, 'energy.txt')
        if not os.path.exists(en_file):
            return
        cycle = 0
        eframes = {}
        with open(en_file, 'r', encoding='utf-8') as energy:
            for idx, eline in enumerate(energy):

                if eline.startswith('#') and 'Time' in eline:
                    continue

                if eline.startswith('#') and 'Cycle' in eline:
                    if idx > 1:
                        # Fill in energy data to the Trajectory
                        if cycle in self.traj_dict[ensemble_name].keys():
                            self.fill_energy(
                                self.traj_dict[ensemble_name][cycle],
                                eframes, ensemble_name, cycle)

                    cycle = int(eline.split()[2].rstrip(','))
                    eframes = {}

                else:
                    edata = eline.split()
                    if not edata:
                        continue
                    eframes[edata[0]] = {}

                    # in case of missing lines in energy.txt
                    if len(edata) < 3:
                        eframes[edata[0]]['potE'] = np.nan
                        eframes[edata[0]]['kinE'] = np.nan
                        continue
                    eframes[edata[0]]['potE'] = float(edata[1])
                    eframes[edata[0]]['kinE'] = float(edata[2])

            # End of file
            # Fill in energy data to the Trajectory
            if cycle in self.traj_dict[ensemble_name].keys():
                self.fill_energy(self.traj_dict[ensemble_name][cycle],
                                 eframes, ensemble_name, cycle)

    def fill_energy(self, traj, eframes, ensemble_name, cycle):
        """Fill in energy-data to a trajectory.

        Parameters
        ----------
        traj: Trajectory object
            Trajectory object to be filled in energy-data.
        eframes: dict
            Dictionary containing energy-data for the trajectory.
        ensemble_name : string
            The name of the ensemble.
        cycle : integer
            Cycle number.

        Updates
        -------
        traj_dict : Trajectory object
            Trajectory object with filled in energy-data.

        """
        if not eframes:
            return
        energy_frame = pd.DataFrame.from_dict(
            eframes, orient='index', dtype=float)
        for energy_term in energy_frame.columns:
            # If dataframe is only nan
            if energy_frame[energy_term].isnull().sum() == \
                    len(energy_frame[energy_term].index):
                return
            traj.frames[energy_term] = energy_frame[energy_term]
            remove_nan(traj.frames[energy_term])
        traj.frames['totE'] = \
            pd.to_numeric(traj.frames['potE']) + \
            pd.to_numeric(traj.frames['kinE'])
        self.traj_dict[ensemble_name][cycle] = traj

    def initialize_compressed(self, input_file):
        """Load PathDensity from a compressed file.

        Parameters
        ----------
        input_file : string
            The input file.

        """
        in_file = os.path.join(self.basepath, input_file)
        if os.path.isfile(in_file) and input_file.endswith('hdf5'):
            self.traj_dict = pd.read_hdf(in_file, key='data')
            self.infos = pd.read_hdf(in_file, key='infos')
        else:
            raise FileNotFoundError(f'File {in_file} not valid.')


class PathVisualize:
    """Class to define the visualization of data with PathDensity.

    Class definition of the visualization of data gathered from simulation
    directory using the PathDensity class.

    """

    def __init__(self, basepath='.', pfile=None):
        """Initialize the PathVisualize class.

        If a supported compressed input file is present, loads the pre-compiled
        data from it. Else, must use specific functions explicitly.

        Parameters
        ----------
        basepath : string, optional
            The path of the input file.
        pfile : string, optional
            The input file.

        """
        self.infos, self.traj_data = {}, {}
        self.op_labels = []
        self.x_list, self.y_list, self.z_list = [], [], []
        self.data_origin = []
        if pfile is None:
            self.pfile = None
        elif pfile.endswith('rst'):
            self.pfile = pfile
        else:
            self.pfile = os.path.join(basepath, pfile)
            self.load_whatever()

    def load_whatever(self):
        """Load all possible supported files.

        This functions directs traffic towards the real loaders.
        Essentially, it does almost nothing.

        """
        if not os.path.isfile(self.pfile):
            raise FileNotFoundError(f'{self.pfile} does not exist.')
        tmp = None
        if self.pfile.endswith('.zip'):
            origin = self.pfile
            self.pfile = self.pfile.rstrip('.zip')
            # If a zip file has been loaded through the file menu
            tmp = os.path.join(os.path.dirname(os.path.abspath(origin)),
                               os.path.splitext(self.pfile)[0] + '.tmp' +
                               os.path.splitext(self.pfile)[1])
            self.pfile = tmp
            pyvisa_unzip(origin, tmp)
        if self.pfile.endswith('.hdf5'):
            self.load_hdf5()
            if tmp:
                os.remove(tmp)
        else:
            raise ValueError(f'Format of {self.pfile} not recognised')

    def load_hdf5(self):
        """Load precompiled data from a hdf5 file.

        Function that loads precompiled data from a .hdf5 file made
        using pandas.

        """
        self.traj_data = pd.read_hdf(self.pfile, key='data')
        self.infos = pd.read_hdf(self.pfile, key='infos')
        self.op_labels = self.infos['op_labels']

    def load_traj(self, criteria):
        """Load relevant data from Trajectories.

        Parameters
        ----------
        criteria : dict
            Dictionary of the selection criteria for which data to load.
            It contains:

            * `x`: string
              Name of parameter to plot.
            * `y`: string
              Name of parameter to plot.
            * `z`: string
              Name of parameter to plot.
            * `ensemble_name`: string
              Name of ensemble to loop through.
            * `cycles`: tuple
              Cycles to loop over.
            * `status`: string
              Status of the path: accepted/rejected.
            * `MC-move`: string, optional
              Generation move of the trajectory.
            * `stored`: bool, optional
              True if the trajectory has available trajectory files.
            * `weights`: bool, optional
              Option to apply statistical weights to the trajectories,
              used in the weighted density plot.

        Returns
        -------
        x_list : list
            List of data from chosen parameter.
        y_list : list
            List of data from chosen parameter.
        z_list : list
            List of data from chosen parameter, if required.
        data_origin: list
            List of ensemble_name and cycle for each point, if required.

        """
        x_list, y_list, z_list, data_origin = [], [], [], []
        for cyc in range(criteria['cycles'][0], criteria['cycles'][1] + 1):
            if cyc in self.traj_data[criteria['ensemble_name']].keys():
                traj = self.traj_data[criteria['ensemble_name']][cyc]
                # If the trajectory is empty
                if 'length' not in traj.info.keys():
                    continue  # pragma: no cover

                stat_weight = traj.info['weight'] if criteria.get('weight',
                                                                  False) else 1

                for key in criteria.keys():
                    # For different frequencies of outputs
                    if key in ['x', 'y', 'z']:
                        if criteria[key] not in ['time', 'cycle', 'None'] \
                                and criteria[key] not in traj.frames.columns:
                            break  # pragma: no cover
                    if key in traj.info.keys():
                        if key == 'status':
                            if criteria[key] == 'BOTH':
                                continue
                            if criteria[key] == 'REJ' and \
                                    traj.info[key] != 'ACC':
                                continue
                        elif ((key == 'stored' and not criteria[key]) or
                              (key == 'MC-move' and criteria[key] == 'All') or
                              (key in ['x', 'y', 'z', 'cycles', 'weight'])):
                            continue

                        if criteria[key] != traj.info[key]:
                            break
                else:
                    for _ in range(stat_weight):
                        data_origin.extend([[traj.info['ensemble_name'],
                                             traj.info['cycle']]] *
                                           traj.info['length'])
                        for xyz in {'x': x_list,
                                    'y': y_list,
                                    'z': z_list}.items():

                            if criteria[xyz[0]] == 'time':
                                if 'time' in traj.frames.columns:
                                    xyz[1].extend(traj.frames['time'])
                                else:
                                    xyz[1].extend(
                                        range(0, traj.info['length']))
                            elif criteria[xyz[0]] == 'cycle':
                                xyz[1].extend([traj.info['cycle']] *
                                              traj.info['length'])
                            elif criteria[xyz[0]] == 'None':
                                xyz[1].extend([1] * traj.info['length'])
                            else:
                                xyz[1].extend(
                                    traj.frames[criteria[xyz[0]]])
        return x_list, y_list, z_list, data_origin

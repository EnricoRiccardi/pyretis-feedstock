# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""
Common functions for the path density.

Functions used to compare and process data, such as matching similar
lists or attempting periodic shifts of values.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

find_rst_file (:py:func: `.find_rst_file`)
    Search for a rst-file from a chosen subdirectory.

read_traj_txt_file (:py:func: `.read_traj_txt_file`)
    Read the sequence of files in a trajectory from a traj.txt file.

recalculate_all (:py:func:`.recalculate_all`)
    Recalculate order parameter and new collective variables by finding
    all trajectory files from a simulation.

shift_data (:py:func: `.shift_data`)
    Finds the median value of a given list of floats, and shifts the
    lower half of the data by the median.

try_data_shift (:py:func: `.try_data_shift`)
    Takes in two lists of values, ``x`` and ``y``, and calculates a
    linear regression and R**2-correlation of the data set. Attempts a
    shift of each data set by their respective median to increase the
    correlation.

where_from_to (:py:func: `.where_from_to`)
    Check the initial and final steps of a trajectory with respect to
    the provided interfaces.

get_cv_names (:py:func: `.get_cv_names`)
    Outputs a list of the names of the descriptors in the simulation.

recalculate_all (:py:func: `.recalculate_all`)
    Recompute all the order parameters according to the PyRETIS storage
    scheme or for individual files/folders.

find_data  (:py:func: `.find_data`)
    Find suitable frames/trajectories to recompute the order parameter on.

"""
import logging
import os
import timeit

import scipy  # pylint: disable=import-error

from pyretis.initiation.initiate_load import write_order_parameters
from pyretis.inout import settings
from pyretis.inout.common import create_backup, TRJ_FORMATS
from pyretis.inout.formats.path import PathExtFile
from pyretis.orderparameter import expand_order_names
from pyretis.setup.common import create_orderparameter
from pyretis.tools.recalculate_order import recalculate_order

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

__all__ = ['find_rst_file', 'read_traj_txt_file',
           'shift_data', 'try_data_shift', 'where_from_to',
           'get_cv_names', 'recalculate_all', 'find_data']


def try_data_shift(x, y, fixedx):
    """Check if shifting increases correlation.

    Function that checks if correlation of data increases by shifting
    either sets of values, x or y, or both. Correlation is checked by
    doing a simple linear regression on the different sets of data:
    - x and y , x and yshift, xshift and y, xshift and yshift.
    If linear correlation increases (r-squared value), data sets are
    updated.

    As a precaution, no shift
    is performed on x values if they are of the first order parameter 'op1'.

    Parameters
    ----------
    x, y : list
        Floats, data values
    fixedx : bool
        If True, x is main OP and should not be shifted.

    Returns
    -------
    x, y : list
        Floats, updated (or unchanged) data values
        (If changed, returns x_temp or y_temp or both)

    """
    # The unshifted data
    _, _, r_val, _, _ = scipy.stats.linregress(x, y)
    # The Y-shifted data
    y_temp = shift_data(y)
    _, _, r_y, _, _ = scipy.stats.linregress(x, y_temp)
    yshift = r_y**2 > r_val**2
    # The X-shifted data
    x_temp = shift_data(x)
    _, _, r_x, _, _ = scipy.stats.linregress(x_temp, y)
    xshift = r_x**2 > r_val**2 and r_x**2 > r_y**2
    # Comparing effectiveness of both shifts individually, and combined
    _, _, r_xy, _, _ = scipy.stats.linregress(x_temp, y_temp)
    xyshift = r_xy**2 > r_val**2 and r_xy**2 > r_y**2 and r_xy**2 > r_x**2

    # If first op is op1, don't shift data
    if xyshift and not fixedx:
        return x_temp, y_temp
    if xshift and not fixedx:
        return x_temp, y
    if yshift:
        return x, y_temp
    return x, y


def shift_data(x):
    """Shifts the data under the median.

    Function that takes in a list of data, and shifts all values
    below the median value of the data by the max difference,
    effectively shifting parts of the data periodically in order to
    give clusters for visualization.

    Parameters
    ----------
    x : list
        Floats, data values

    Returns
    -------
    xnorm : list
        Floats where some values are shifted values of x,
        and some are left unchanged.

    """
    xmin, xmax = min(x), max(x)
    xnorm = []
    # The max difference in x-data
    diff_x = xmax - xmin
    # The Median of x-data
    medix = xmin + 0.5 * diff_x

    for i in x:
        if i < medix:
            xnorm.append(i + diff_x)
        else:
            xnorm.append(i)
    return xnorm


def read_traj_txt_file(path):
    """Read a traj.txt file.

    Function which reads a traj.txt file and returns a dict containing
    the name of each file in the trajectory and the sign of their velocity.

    Parameters
    ----------
    path : string
        Path to the traj.txt file.

    Returns
    -------
    files : dict
        Dictionary containing each file in the trajectory and
        the sign of their velocity.

    """
    files = {}
    i = 0
    with PathExtFile(path, 'r') as pfile:
        for block in pfile.load():
            for data in block['data']:
                if data[0] == '0':
                    files[i] = [data[1], data[3]]

                if data[1] != files[i][0]:
                    i += 1
                    files[i] = [data[1], data[3]]
    return files


def find_rst_file(search_dir):
    """Search for rst-files.

    Parameters
    ----------
    search_dir : string
        Path to the .rst file.

    Returns
    -------
    out[0] : string
        Path and name of the .rst file.

    """
    for _ in range(len(search_dir.split('/'))):
        for file_name in os.listdir(search_dir):
            if file_name.endswith('rst'):
                os.chdir(search_dir)
                return search_dir + '/' + file_name
        os.chdir('../')
        search_dir = os.getcwd()
    return search_dir


def where_from_to(trj, int_a, int_b=float('-inf')):
    r"""Detect L∕R starts and L / R / \* ends.

    Given a list of order parameters (a trj), the function
    will try to establish where the path started (L or R or \*)
    and where it ended.
    Note: for the 'REJ' paths, this function results might differ
    from PyRETIS.

    Parameters
    ----------
    trj: numpy array
        The order parameters of the trj.
    int_a: float
        The interface that defines state A.
    int_b: float, optional
        The interface that defines state B. If not given, it is assumed
        that the 0^- ensemble is in use without the 0^- L interface.

    Returns
    -------
    start: string\*1
        The initial position of the trajectory in respect to the
        interfaces given (L eft, R ight or \* for nothing).
    end: string\*1
        The final position of the trajectory in respect to the
        interfaces given (L eft, R ight or \* for nothing).

    """
    start, end = '*', '*'
    int_l = min(int_a, int_b)
    int_t = max(int_a, int_b)
    if trj[0] >= int_t:
        start = 'R'
    if trj[0] < int_l:
        start = 'L'
    if trj[-1] >= int_t:
        end = 'R'
    if trj[-1] < int_l:
        end = 'L'

    return start, end


def get_cv_names(input_settings, num_columns=None):
    """Return labels for the order parameter and collective variables.

    The labels follow the same convention for the main order parameter
    and any extra collective variable: each ``[orderparameter]`` /
    ``[collective-variable]`` block may declare ``name`` either as a
    single string or as a list of strings.

    * A list of strings is used as-is for that block (one label per
      value returned by the corresponding ``calculate()``).
    * A single string is expanded with an index suffix
      (``"<name>_1", "<name>_2", ...``) to match the number of values
      the block produces. If the block produces a single value the
      bare string is used.
    * When ``name`` is missing the labels fall back to ``"op_<i>"`` for
      the main order parameter and ``"cv_<i>"`` for each collective
      variable (using the bare prefix if the block produces a single
      value).

    Parameters
    ----------
    input_settings : dict
        Dictionary with the settings from the simulations.
    num_columns : int, optional
        Total number of order-parameter columns observed (e.g. read
        from ``order.txt``). When supplied, this is used to allocate
        unspecified blocks: if there is exactly one section configured
        and ``num_columns`` is greater than one, the single ``name`` is
        expanded with indices to cover all columns. When ``num_columns``
        is given and cannot be reconciled with the configured names,
        a :py:class:`ValueError` is raised.

    Returns
    -------
    names : list of str
        Flat list of column labels, in the same order as the
        concatenated ``calculate()`` outputs.

    Raises
    ------
    ValueError
        If a ``name`` list does not match the number of values
        produced by its block, or if ``num_columns`` cannot be
        reconciled with the configured names.

    """
    blocks = []
    main_block = input_settings.get('orderparameter')
    if main_block is not None:
        blocks.append(('op', main_block, main_block.get('name')))
    cv_blocks = input_settings.get('collective-variable', [])
    for idx, c_v in enumerate(cv_blocks, start=1):
        prefix = 'cv' if len(cv_blocks) == 1 else f'cv{idx}'
        blocks.append((prefix, c_v, c_v.get('name')))

    if not blocks:
        return []

    # When the only information available is the column count from a
    # file (e.g. order.txt) and a single section is configured, that
    # section claims all columns.
    if num_columns is not None and len(blocks) == 1:
        prefix, _, name = blocks[0]
        return expand_order_names(name, num_columns, default_prefix=prefix)

    labels = []
    for prefix, _, name in blocks:
        if isinstance(name, (list, tuple)):
            labels.extend(expand_order_names(name, len(name), prefix))
        else:
            labels.extend(expand_order_names(name, 1, prefix))

    if num_columns is not None and len(labels) != num_columns:
        raise ValueError(
            "Configured order-parameter names produce "
            f"{len(labels)} labels, but the data has {num_columns} "
            "columns. Provide a list ``name`` for any multi-valued "
            "order parameter or collective variable."
        )
    return labels


def _write_cycle(functions, cycles, main_order):
    """Write order parameters for a single trajectory cycle.

    Parameters
    ----------
    functions : list
        Order-parameter functions returned by ``create_orderparameter``.
    cycles : dict
        Cycle entry from the trajectory dict (contains ``'traj'`` key).
    main_order : str or None
        Path to the ensemble-level ``order.txt`` to append to, or
        ``None`` if no aggregation file is needed.

    """
    here = os.path.dirname(os.path.abspath(cycles['traj'][0]))
    local_order = cycles.get('o_txt', os.path.join(here, 'order.txt'))
    results_dict = []
    for trj in cycles['traj']:
        try:
            for order_p in recalculate_order(functions, trj, {}):
                results_dict.append(order_p)  # pragma: no cover
        except (KeyError, AttributeError):
            logger.warning('File %s not valid', trj)
    create_backup(local_order)
    write_order_parameters(local_order, results_dict,
                           cycles.get('header', 'Recalculated data'))
    if main_order is not None:
        with open(local_order, 'rb') as src, open(main_order, 'ab') as dst:
            dst.write(src.read())


def _process_ensemble(functions, ens, runfolder):
    """Process one ensemble: recompute order parameters for every cycle.

    Parameters
    ----------
    functions : list
        Order-parameter functions returned by ``create_orderparameter``.
    ens : dict
        Ensemble entry from the trajectory dict.
    runfolder : str
        Root path of the simulation run.

    """
    main_order = None
    if ens.get('main_o') is not None:
        main_order = os.path.join(runfolder, ens['main_o'])
        create_backup(main_order)
    for _, cycles in sorted(ens['traj'].items()):
        _write_cycle(functions, cycles, main_order)


def recalculate_all(runfolder, iofile, ensemble_names=None, data=None):
    """Recalculate order parameter and collective variables.

    Performs post-processing by analyzing trajectories of old simulations
    to extract data, do new calculations, and write to a new order.txt file.

    Parameters
    ----------
    runfolder: string
        The path of the execution directory.
    iofile: string
        The input file where the settings are collected.
    ensemble_names: list, optional
        List of ensemble names in the simulation to work with.
    data: string, optional
        If given, the function will check only the single file or
        look only in the given directory

    Returns
    -------
    out : boolean
        True if the recomputation was successful, False otherwise.

    """
    if iofile is None:
        raise ValueError('Input file not given')
    input_settings = settings.parse_settings_file(
        os.path.join(runfolder, iofile)
    )
    trj_dict = find_data(runfolder, ensemble_names, data=data)

    logger.progress('Re-computing the collective variables.')
    tic = timeit.default_timer()
    if not trj_dict:
        logger.warning(
            'No data to re-compute from in %s — no trajectory data found. '
            'Ensure the run folder contains numeric ensemble subdirectories '
            '(e.g. 000, 001) or pass the `data` argument pointing to a '
            'file or directory to re-process.', runfolder)
        return False

    try:
        functions = create_orderparameter(input_settings)
    except (ImportError, ValueError):
        logger.warning('Invalid Order Parameter')  # pragma: no cover
        return False  # pragma: no cover

    for _, ens in trj_dict.items():
        _process_ensemble(functions, ens, runfolder)

    logger.progress('# Data successfully recomputed!')
    logger.progress('# Time spent: %.2fs', timeit.default_timer() - tic)
    return True


def find_data(runfolder, ensemble_names=None, data=None):
    """Find the trajectory data used to do post-processing.

    find_traj returns a dict with a structure resembling that
    of the simulation.

    Parameters
    ----------
    runfolder: string, optional
        The path of the execution directory.
    ensemble_names: list, optional
        List of ensemble names in the simulation to work with.
    data: string, optional
        If given, the function will check only the single file or
        look only in the given directory

    Returns
    -------
    trj_dict : dict
        To each key, ensemble_name (e.g. 000, 001, etc)
        the values are: the last accepted trajectories given by the
        `accepted`-key; the generation trajectory or conf files given by
        the `generation`-key, and lastly the dictionary `stored_traj`
        that is given by the `traj`-key. `stored_traj` is split up into
        the dictionaries`traj-acc` and `traj-rej` which have keys for
        all the accepted and rejected cycles respectively,
        where the trajectory files for that cycle is stored.

    """
    trj_dict = {}
    # Specific file/folder
    if data is not None:
        sources = [data] if os.path.isfile(data) else _get_trjs(data)
        if sources:
            trj_dict['000'] = {'traj': {'0': {'traj': sources}}}
        return trj_dict

    # Structured data
    flag_map = {'traj-acc': 'ACC', 'traj-rej': 'REJ'}

    if ensemble_names is None:
        ensemble_names = [i.name for i in os.scandir(runfolder) if
                          (i.is_dir() and i.name.isdigit())]

    for ens_name in sorted(ensemble_names):
        order_txt = os.path.join(os.path.abspath(runfolder), ens_name,
                                 'order.txt')
        trj_dict[ens_name] = {'main_o': order_txt, 'traj': {}}
        for trj_type in ['traj-acc', 'traj-rej']:
            here = os.path.join(os.path.abspath(runfolder),
                                ens_name, 'traj', trj_type)
            if not os.path.isdir(here):
                continue
            for cycle in sorted([i.name for i in os.scandir(here) if
                                i.name.isdigit()]):
                loc = os.path.join(here, cycle, 'traj')
                o_txt = os.path.join(here, cycle, 'order.txt')
                header = f'# Cycle: {cycle},' \
                         f' status: {flag_map[trj_type]}'
                if os.path.exists(o_txt):
                    with open(o_txt, 'r', encoding='utf-8') as file_in:
                        header = file_in.readline().replace('\n', '')

                trj_dict[ens_name]['traj'][cycle] = {'header': header,
                                                     'o_txt': o_txt,
                                                     'traj': _get_trjs(loc)}
    return trj_dict


def _get_trjs(runfolder='.'):
    """Find the trajectory files.

    Parameters
    ----------
    runfolder : string, optional
        the location of the main simulation folder.

    Returns
    -------
    full_name_trj : list
        The trajectory files contained in the folder.

    """
    trj = [i.name for i in os.scandir(runfolder) if i.name[-4:] in TRJ_FORMATS]
    full_name_trj = [os.path.join(runfolder, i) for i in trj]
    return full_name_trj

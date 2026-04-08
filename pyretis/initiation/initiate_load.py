# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Method for initiating paths by loading from previous simulation(s).

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

initiate_load (:py:func:`.initiate_load`)
    A method that will get the initial path from the output of
    a previous simulation.
"""
import collections
import errno
import logging
import os
import shutil
from pyretis.core.common import compute_weight
from pyretis.core.path import Path, paste_paths
from pyretis.core.pathensemble import generate_ensemble_name
from pyretis.initiation.initiate_kick import initiate_path_ensemble_kick
from pyretis.inout import print_to_screen
from pyretis.inout.common import make_dirs, TRJ_FORMATS
from pyretis.inout.formats.order import OrderPathFile, OrderFile
from pyretis.inout.formats.energy import EnergyPathFile
from pyretis.inout.formats.path import PathIntFile, PathExtFile
from pyretis.inout.formats.xyz import read_xyz_file, convert_snapshot
from pyretis.inout.formats.gromacs import read_gromacs_generic
from pyretis.inout.formats.cp2k import read_cp2k_box
from pyretis.tools.recalculate_order import recalculate_order
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


__all__ = ['initiate_load', 'write_order_parameters', 'clean_path',
           'read_path_files', 'read_path_files_ext', 'reorderframes',
           '_check_path', '_do_the_dirty_load_job']


def initiate_load(simulation, settings, cycle, plot_loads=False):
    """Initialise paths by loading already generated ones.

    Parameters
    ----------
    simulation : object like :py:class:`.Simulation`
        The simulation we are setting up.
    cycle : integer
        The simulation cycles we are starting at.
    settings : dictionary
        A dictionary with settings for the initiation.

    """
    folder = settings['initial-path'].get('load_folder', 'load')

    if not os.path.exists(folder):
        logger.critical('Load folder "%s" not found!', folder)
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                folder)

    # Check if there is data in the main folder
    if not os.path.exists(os.path.join(folder, 'traj.txt')):
        listoffiles = [i.name for i in os.scandir(folder) if i.is_file()]
        traj_files = [i for i in listoffiles if i[-4:] in TRJ_FORMATS]

        # If there are files, we shall compute their traj.txt and order.txt
        # to be used, if needed, later.
        if traj_files and 'order.txt' not in listoffiles:
            # This generated the traj.txt in the main load folder
            # We here borrow the 000 ensemble functions. It is a rather
            # safe assumption.
            engine = simulation.ensembles[0]['engine']
            system = simulation.ensembles[0]['system']
            order_function = simulation.ensembles[0]['order_function']
            if 'traj.txt' not in listoffiles:
                _generate_traj_txt_from_ext(folder, system,
                                            accepted=False)

            # This generated the order.txt in the main load folder
            if engine.engine_type == 'internal':
                traj = _load_trajectory(folder)
                _load_order_parameters(traj, folder, system, order_function)
            else:
                traj = _load_external_trajectory(folder, engine, copy=False)
                _load_order_parameters_ext(traj, folder, order_function,
                                           engine)

    for ensemble, set_ens in zip(simulation.ensembles, settings['ensemble']):
        path_ensemble = ensemble['path_ensemble']
        engine = ensemble['engine']
        system = ensemble['system']
        order_function = ensemble['order_function']
        name = path_ensemble.ensemble_name
        logger.info('Loading data for path path ensemble %s:', name)
        print_to_screen(f'Loading data for path path ensemble {name}')
        engine.exe_dir = path_ensemble.directory['generate']
        path = Path(simulation.rgen, maxlen=None)
        path.generated = ('re', 0, 0, 0)  # It remains for formatted loads.
        edir = os.path.join(
            folder,
            generate_ensemble_name(path_ensemble.ensemble_number))
        if not os.path.exists(os.path.join(edir, 'order.txt')):
            logger.info('Input Load missing')
            print_to_screen('No order.txt file found in %, attempt generate!',
                            edir)
            # Make sure the files are where they are supposed to be:
            _do_the_dirty_load_job(folder, edir)
            path.set_move('ld')

        if not os.path.exists(os.path.join(edir, 'traj.txt')):
            # Generate a traj.txt file for the path ensemble:
            _generate_traj_txt_from_ext(edir, system)
            path.set_move('ld')

        if engine.engine_type == 'internal':
            accept, status = read_path_files(
                path,
                edir,
                ensemble
            )
        elif engine.engine_type == 'external':
            accept, status = read_path_files_ext(
                path,
                edir,
                ensemble
            )
        else:
            raise ValueError('Unknown engine type!')

        # If a path is not acceptable, than here we try to rearrange it such
        # to satisfy each ensemble definition
        # Load paths can be accepted while still having many points outside
        # of the defined region. So even if accepted, we need to clean if the
        # path type is 'ld'.
        simtype = set_ens['simulation']['task']
        print("I am this simtype:", simtype)
        if accept and path.get_move() == 'ld':  # Some cleaning if needed
            clean_path(path, path_ensemble, simtype=simtype)
            # The path should still be checked, as something funny might occur
            accept, status = _check_path(path, path_ensemble)
            assert accept, "Path should be accepted after cleaning, but" + \
                           f" it is not! Status: {status}"
        elif set_ens['simulation']['task'] == 'explore':
            set_ens['tis']['shooting_move'] = 'exp'
            clean_path(path, path_ensemble, simtype=simtype)
            status = 'EXP'
        elif not accept:  # Let's try to fix the path
            msg = 'Sorting input files'
            logger.info(msg)
            print_to_screen(msg)
            path, (accept, status) = reorderframes(path, path_ensemble)
            clean_path(path, path_ensemble, simtype=simtype)
            path.set_move('ld')
            if not accept:
                msg = f'Path stored in load folder {edir} does not satisfy the'
                msg += f' relative ensemble ({path_ensemble.ensemble_number}) '
                msg += 'definition. \n The given path has a min of '
                msg += f'{path.ordermin[0]} and a max of {path.ordermax[0]} '
                msg += 'while the relative interface are located at '
                msg += '{}, {}, and {}.\n'.format(*path_ensemble.interfaces)
                if set_ens['initial-path'].get('load_and_kick', False):
                    msg += 'Initialising load_and_kick procedure.'
                    logger.info(msg)
                    print_to_screen(msg)
                    best_frame = path.phasepoints[0]
                    target = ensemble['path_ensemble'].interfaces[1]
                    for phasepoint in path.phasepoints:
                        if abs(best_frame.order[0] - target) > \
                                abs(phasepoint.order[0] - target):
                            best_frame = phasepoint

                    ensemble['system'] = best_frame
                    accept, path, status = initiate_path_ensemble_kick(
                        ensemble,
                        set_ens['tis'],
                        cycle)
                else:
                    msg += 'You might try to add the load_and_kick = True' \
                           ' option in the initial path section to try to fix'\
                           ' this. \n' \
                           'NOTE: points within state B are not used.'
                    raise ValueError(msg)

        path.status = status
        if set_ens['tis'].get('high_accept', False):
            move = set_ens['tis'].get('shooting_move', 'sh')
            if move == 'wf':
                intf = ensemble['path_ensemble'].interfaces
                int_cap = set_ens['tis'].get('interface_cap', intf[2])
                intf = [intf[0], intf[1], int_cap]
            if set_ens['tis'].get('shooting_move', 'sh') in ('ss', 'wf'):
                path.weight = compute_weight(path, intf, move)
        path_ensemble.add_path_data(path, status, cycle)

        # As the loading is sometimes messy, I think it's not a bad idea to
        # just plot the loaded paths, such that the user KNOWS what is loaded.
        # ignore this for coverage as it is solely for debugging purposes.
        if plot_loads:  # pragma: no cover
            left, center, right = path_ensemble.interfaces
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots()
            ax.plot([p.order[0] for p in path.phasepoints], lw=1, marker='o',
                    ms=3)
            ax.axhline(y=left, ls='--', c='k')
            ax.axhline(y=center, ls='--', c='k')
            ax.axhline(y=right, ls='--', c='k')
            ax.set_xlabel('Phasepoint index')
            ax.set_ylabel('Order parameter')
            fig.savefig(f'clean_ensemble{path_ensemble.ensemble_number}.png')

        yield accept, path, status, path_ensemble


def reorderframes(path, path_ensemble):
    """Re-order the phase points in the input path.

    This function assumes that the path needs to be fixed.
    It checks for the validity of the path and keeps only the frames
    that are useful for the region to explore e.g. [0^-].

    Parameters
    ----------
    path : object like :py:class:`.PathBase`
        The path we try to fix.
    path_ensemble : object like :py:class:`.PathEnsembleExt`
        The path ensemble the path could be added to.

    """
    slave_path = path.copy()
    left, center, right = path_ensemble.interfaces
    # We first check that sorting doesn't suffice to fix the path
    path.sorting('order', reverse=True)
    accept, status = _check_path(path, path_ensemble, warning=False)
    if accept:
        return path, (accept, status)
    path.sorting('order', reverse=False)
    accept, status = _check_path(path, path_ensemble, warning=False)
    if accept:
        return path, (accept, status)

    # If the path is from L to R, we order the phase points according
    # to the order parameter.
    # If the path is from R to L (e.g. for the [0^-] ensemble),
    # a reversed path is generated.
    # If the path goes from L and stops, it is doubled to create a
    # L to L path.
    # If the path goes from R and stops, it is doubled to create a
    # R to R path.
    # NOTE: the path is now sorted.
    if 'R' == path_ensemble.start_condition and path.ordermin[0] < center:
        slave_path.sorting('order', reverse=False)
        path = paste_paths(path, slave_path)

    elif "L" == path_ensemble.start_condition and\
            path.get_end_point(left, right) not in {'L', 'R'} and\
            path.ordermax[0] >= center and path.ordermin[0] < left:
        slave_path = path.copy()
        path.sorting('order', reverse=True)
        slave_path.sorting('order', reverse=True)
        path = paste_paths(path, slave_path)
    elif set(path_ensemble.start_condition) == set(['L', 'R']):
        if path.ordermax[0] >= center and path.ordermin[0] < left:
            slave_path = path.copy()
            path.sorting('order', reverse=True)
            slave_path.sorting('order', reverse=True)
            path = paste_paths(path, slave_path)

    path.set_move('ld')

    return path, _check_path(path, path_ensemble)


def clean_path(path, path_ensemble, simtype='retis'):
    """Remove useless frames from a path.

    This function removes phasepoints from a path that are
    not needed (e.g. the one in stable states)

    Parameters
    ----------
    path : object like :py:class:`.PathBase`
        The path we try to fix.
    path_ensemble : object like :py:class:`.PathEnsembleExt`
        The path ensemble the path could be added to.
    simtype : string
        The simulation type. 'retis', 'repptis', 'pptis'

    """
    left, center, right = path_ensemble.interfaces

    # Initialize
    dists = {'left': float('inf'),
             'right': float('inf')}
    ph_min = path.phasepoints[path.ordermin[1]]
    ph_max = path.phasepoints[path.ordermax[1]]
    keep_ph = None

    # First, find the two nearest crossing points to the 0 interface.
    if path_ensemble.start_condition == "L":
        int_opt = (left, left)
        int_a, int_b = center, right
        keep_ph = path.phasepoints[path.ordermax[1]].copy()
    elif path_ensemble.start_condition == "R":
        int_opt = (center, center)
        int_a, int_b = left, center
        keep_ph = path.phasepoints[path.ordermin[1]].copy()
    else:
        int_opt = (left, right)
        int_a, int_b = left, right

    for i, phasepoint in enumerate(path.phasepoints):
        # From left, searching for closest point left of int_opt[0]
        if phasepoint.order[0] < int_opt[0]:
            if int_opt[0] - phasepoint.order[0] < dists['left']:
                dists['left'] = int_opt[0] - phasepoint.order[0]
                ph_min = phasepoint.copy()
        elif simtype == 'retis':  # RETIS behavior
            if phasepoint.order[0] - int_opt[1] < dists['right']:
                dists['right'] = phasepoint.order[0] - int_opt[1]
                ph_max = phasepoint.copy()
        elif simtype in ['pptis', 'repptis']:
            if set(["L", "R"]) == set(path_ensemble.start_condition):
                if phasepoint.order[0] > int_opt[1]:
                    if phasepoint.order[0] - int_opt[1] < dists['right']:
                        dists['right'] = phasepoint.order[0] - int_opt[1]
                        ph_max = phasepoint.copy()
            else:
                if phasepoint.order[0] - int_opt[1] < dists['right']:
                    dists['right'] = phasepoint.order[0] - int_opt[1]
                    ph_max = phasepoint.copy()
    # Prepare for deleting unnecessary frames.
    erase_list = []
    for i, phasepoint in enumerate(path.phasepoints):
        if not int_a <= phasepoint.order[0] <= int_b:
            erase_list.append(i)

    if not erase_list:
        return

    # Erase the largest elements first to avoid counting problems.
    for i in erase_list[::-1]:
        path.delete(i)

    # The path has to cope with the ensemble requirements.
    if path_ensemble.start_condition == "R":
        path.phasepoints.insert(0, ph_min.copy())
        path.phasepoints.insert(0, ph_max.copy())
        # Add the best point there is.
        if path.ordermin[0] >= center:
            path.phasepoints.append(keep_ph)
        path.phasepoints.append(ph_min.copy())
        path.phasepoints.append(ph_max.copy())

    elif path_ensemble.start_condition == 'L':
        path.phasepoints.insert(0, ph_max.copy())
        path.phasepoints.insert(0, ph_min.copy())
        # Add the best point there is.
        if path.ordermax[0] < center:
            path.phasepoints.append(keep_ph)
        path.phasepoints.append(ph_max.copy())
        path.phasepoints.append(ph_min.copy())

    else:
        path.phasepoints.insert(0, ph_min.copy())
        path.phasepoints.append(ph_max.copy())

    if simtype in ['pptis', 'repptis']:
        if ph_max.order[0] > right:
            if set(["L", "R"]) == set(path_ensemble.start_condition):
                path.sorting(key="order")
                logger.info("The load path has been sorted")


def _do_the_dirty_load_job(mainfolder, edir):
    """Copy the files to a place where PyRETIS expects them to be.

    The function checks if in the destination folder some suitable files
    are already present. If not, it tries to copy them from the main load
    folder.

    Parameters
    ----------
    mainfolder : string
        The path of the main load directory with frame/trajectory
        files.
    edir : string
        The path of the path ensemble load directory with
        frame/trajectory files.

    """
    listoffiles = [i.name for i in os.scandir(mainfolder) if i.is_file()]
    traj_files = [i for i in listoffiles if i[-4:] in TRJ_FORMATS]

    # If there is an ensemble folder, let's check what there is inside.
    if os.path.exists(edir):
        # First in the main ensemble folder
        listoffiles = [i.name for i in os.scandir(edir) if i.is_file()]
        traj_files_ens = [i for i in listoffiles if i[-4:] in TRJ_FORMATS]

        # Then in an eventual accepted folder
        traj_ens_acc = []
        if os.path.exists(os.path.join(edir, 'accepted')):
            listoffiles = [i.name for i in os.scandir(os.path.join(edir,
                                                                   'accepted'))
                           if i.is_file()]
            traj_ens_acc = [i for i in listoffiles if i[-4:] in TRJ_FORMATS]

        # If we pass the following, it means no input exist.
        if not traj_files_ens and not traj_ens_acc and not traj_files:
            logger.critical('No files to load in %s', edir)
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                    edir)

        # Here we are sure that there are inputs, we shall act being carefull
        # to not overwrite eventual info.
        # Starting from the outermost directory:
        if traj_ens_acc or traj_files_ens:
            # If there are files in the accepted folder, we just keep them.
            # If they are in the edir folder,
            # they shall be moved to the accepted folder.
            make_dirs(os.path.join(edir, 'accepted'))
            for filei in traj_files_ens:
                if filei not in {'order.txt', 'traj.txt', 'energy.txt',
                                 'pathensemble.txt'}:
                    shutil.move(
                        os.path.join(edir, filei),
                        os.path.join(edir, 'accepted', filei)
                    )
            return

    # If we arrived here, it means that we have to use the files
    # present in the main folder. So let's do the dirty copy.
    if not traj_files:
        logger.critical('No files to load in %s', mainfolder)
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                mainfolder)

    make_dirs(os.path.join(edir, 'accepted'))
    for filei in traj_files:
        if filei in {'order.txt', 'traj.txt', 'energy.txt'}:
            shutil.copy(
                os.path.join(mainfolder, filei),
                os.path.join(edir, filei)
            )
        else:
            shutil.copy(
                os.path.join(mainfolder, filei),
                os.path.join(edir, 'accepted', filei)
            )


def _generate_traj_txt_from_ext(dirname, system, accepted=True):
    """Generate a traj.txt file in the path ensemble folder.

    Parameters
    ----------
    dirname : string
        The path of the path ensemble directory with
        trajectory files.
    system : object like :py:class:`.System`
        A system object we can use when calculating the order
        parameter(s).
    accepted : boolean, optional
        It controls where the traj files shall be located,
        if True in the accepted folder, if False in the main folder.

    """
    load_file = os.path.join(dirname, 'traj.txt')

    if accepted:
        dirname = os.path.join(dirname, 'accepted')

    traj_files = sorted([i for i in os.listdir(dirname)
                         if i[-4:] in TRJ_FORMATS])

    path = Path()
    # Load Sparse works also from multiple initial trajectories
    # both internally and externally
    if system.particles.particle_type == 'internal':
        for filei in traj_files:
            if filei[-4:] == '.xyz':
                for snap in read_xyz_file(os.path.join(dirname, filei)):
                    _, pos, vel, _ = convert_snapshot(snap)
                    snapshot = {'pos': pos, 'vel': vel}
                    phase_point = system_from_snapshot(system, snapshot)
                    path.append(phase_point)
            if filei[-4:] == '.txt':
                traj = _load_trajectory(dirname, filei)
                for snapshot in traj['data']:
                    snapshot = {'pos': snapshot['pos'],
                                'vel': snapshot['vel']}
                    phase_point = system_from_snapshot(system, snapshot)
                    path.append(phase_point)

        with PathIntFile(load_file, 'w', backup=True) as pathfile:
            pathfile.output(0, (path, 'ACC'))
    else:
        for tfile in traj_files:
            if tfile[-4:] == '.txt':
                continue
            all_tfile = os.path.join(dirname, tfile)
            if tfile[-4:] == '.xyz':
                for i, snap in enumerate(read_xyz_file(all_tfile)):
                    snapshot = {'pos': (tfile, i), 'vel': False}
                    phase_point = system_from_snapshot(system, snapshot)
                    path.append(phase_point)
            else:
                for i, snap in enumerate(read_gromacs_generic(all_tfile)):
                    snapshot = {'pos': (tfile, i), 'vel': False}
                    phase_point = system_from_snapshot(system, snapshot)
                    path.append(phase_point)
        with PathExtFile(load_file, 'w', backup=True) as pathfile:
            pathfile.output(0, (path, 'ACC'))


def _load_order_parameters(traj, dirname, system, order_function):
    """Load or recalculate the order parameters.

    Parameters
    ----------
    traj : dictionary
        The trajectory we have loaded. Used here if we are
        re-calculating the order parameter(s).
    dirname : string
        The path to the directory with the input files.
    system : object like :py:class:`.System`
        A system object we can use when calculating the order
        parameter(s).
    order_function : object like :py:class:`.OrderParameter`
        This can be used to re-calculate order parameters in case
        they are not given.

    Returns
    -------
    out : list
        The order parameters, each item in the list corresponds to a
        time frame.

    """
    order_file_name = os.path.join(dirname, 'order.txt')
    try:
        with OrderPathFile(order_file_name, 'r') as orderfile:
            print_to_screen('Loading order parameters.')
            order = next(orderfile.load())
            return order['data'][:, 1:]
    except FileNotFoundError:
        print_to_screen(f'Could not read file: {order_file_name}')
    orderdata = []
    print_to_screen('Recalculating order parameters for input path.')
    logger.info('Recalculating order parameters for input path.')
    for snapshot in traj['data']:
        system.particles.pos = snapshot['pos']
        system.particles.vel = snapshot['vel']
        orderdata.append(order_function.calculate(system))
    return orderdata


def _load_order_parameters_ext(traj, dirname, order_function, engine):
    """Load or re-calculate the order parameters.

    For external trajectories, dumping of specific frames from a
    trajectory might be expensive and we here do slightly more work
    than just dumping the frames.

    Parameters
    ----------
    traj : dictionary
        The trajectory we have loaded. Used here if we are
        re-calculating the order parameter(s).
    dirname : string
        The path to the directory with the input files.
    order_function : object like :py:class:`.OrderParameter`
        This can be used to re-calculate order parameters in case
        they are not given.
    system : object like :py:class:`.System`
        A system object we can use when calculating the order
        parameter(s).

    Returns
    -------
    out : list
        The order parameters, each item in the list corresponds to
        a time frame.

    """
    order_file_name = os.path.join(dirname, 'order.txt')
    if os.path.isfile(order_file_name):
        with OrderPathFile(order_file_name, 'r') as orderfile:
            print_to_screen('Loading order parameters.')
            order = next(orderfile.load())
            return order['data'][:, 1:]
    else:
        print_to_screen('Could not read file: {order_file_name}')

    orderdata = []
    print_to_screen('Recalculating order parameters for input path.')
    logger.info('Recalculating order parameters for input path.')
    # First get unique files and indices for them:
    files = collections.OrderedDict()
    for snapshot in traj['data']:
        filename = snapshot[1]
        idx = int(snapshot[2])
        if filename not in files:
            files[filename] = {'minidx': idx-1, 'maxidx': idx+1,
                               'reverse': snapshot[3]}
        if idx < files[filename].get('minidx', idx+1):
            files[filename]['minidx'] = idx
        if idx > files[filename].get('maxidx', idx-1):
            files[filename]['maxidx'] = idx
    # Ok, now we have the files, calculate the order parameters:
    for filename, info in files.items():
        if filename[-4:] == '.xyz':           # CP2K specific NVE.
            info['box'], _ = read_cp2k_box(engine.input_files['template'])
        for order in recalculate_order(order_function, filename, info):
            orderdata.append(order)

    # Store the re-calculated order parameters so we don't have
    # to re-calculate again later:
    write_order_parameters(order_file_name, orderdata)
    return orderdata


def write_order_parameters(order_file_name, orderdata,
                           header="# Re-calculated order parameters."):
    """Store re-calculated order parameters to a file."""
    with OrderFile(order_file_name, 'w') as orderfile:
        orderfile.write(header)
        for step, data in enumerate(orderdata):
            orderfile.output(step, data)


def _load_energies_for_path(path, dirname):
    """Load energy data for a path.

    Parameters
    ----------
    path : object like :py:class:`.PathBase`
        The path we are to set up/fill.
    dirname : string
        The path to the directory with the input files.

    Returns
    -------
    None, but may add energies to the path.

    """
    # Get energies if any:
    energy_file_name = os.path.join(dirname, 'energy.txt')
    try:
        with EnergyPathFile(energy_file_name, 'r') as energyfile:
            print_to_screen('Loading energy data.')
            energy = next(energyfile.load())
            path.update_energies(energy['data']['ekin'],
                                 energy['data']['vpot'])
    except FileNotFoundError:
        print_to_screen(f'Could not read file: {energy_file_name}')
        print_to_screen('Energies are not set.')


def _check_path(path, path_ensemble, warning=True):
    """Run some checks for the path.

    Parameters
    ----------
    path : object like :py:class:`.PathBase`
        The path we are to set up/fill.
    path_ensemble : object like :py:class:`.PathEnsemble`
        The path ensemble the path could be added to.
    warning : boolean, optional
        If True, it output warnings, else only debug info.

    """
    start, end, _, cross = path.check_interfaces(path_ensemble.interfaces)
    accept = True
    status = 'ACC'
    messages = []

    if start is None or start not in path_ensemble.start_condition:
        messages.append("Initial path for %s starts at the wrong interface!")
        status = 'SWI'
        accept = False
    if end not in ('R', 'L'):
        messages.append("Initial path for %s ends at the wrong interface!")
        status = 'EWI'
        accept = False
    if not cross[1]:
        messages.append(
            "Initial path for %s does not cross the middle interface!")
        status = 'NCR'
        accept = False

    if not accept:
        msg = ' '.join(messages)
        name = path_ensemble.ensemble_name
        if warning:
            logger.critical(msg, *([name] * len(messages)))
        else:
            logger.debug(msg, *([name] * len(messages)))

    path.status = status
    return accept, status


def _load_trajectory(dirname, filename='traj.txt'):
    """Load a trajectory from a file.

    Parameters
    ----------
    dirname : string
        The directory where we can find the trajectory file(s).
    filename : string
        The name of the trajectory file.

    Returns
    -------
    traj : dict
        A dictionary containing the trajectory information.
        That is, position and velocity for each frame.

    """
    with PathIntFile(os.path.join(dirname, filename), 'r') as trajfile:
        # Just get the first trajectory:
        traj = next(trajfile.load())
    return traj


def read_path_files(path, dirname, ensemble):
    """Read data needed for a path from a directory.

    Parameters
    ----------
    path : object like :py:class:`.PathBase`
        The path we are to set up/fill.
    dirname : string
        The path to the directory with the input files.
    ensemble : dictionary of objects
        It contains:

        * `path_ensemble`: object like :py:class:`.PathEnsemble`
          This is the path ensemble to perform the initialization.
        * `order_function`: object like :py:class:`.OrderParameter`
          The class used for obtaining the order parameter(s).
        * `engine`: object like :py:class:`.EngineBase`
          The engine to use for propagating a path.
        * `system`: object like :py:class:`.System`
          A system object we can use when calculating the order parameter(s).

    """
    path_ensemble = ensemble['path_ensemble']
    engine = ensemble['engine']
    system = ensemble['system']
    order_function = ensemble['order_function']
    left, _, right = path_ensemble.interfaces
    traj = _load_trajectory(dirname)
    orderdata = _load_order_parameters(traj, dirname, system, order_function)

    # Add to path:
    print_to_screen('Creating path from files.')
    logger.debug('Creating path from files.')
    for snapshot, orderi in zip(traj['data'], orderdata):
        snapshot = {'order': orderi, 'pos': snapshot['pos'],
                    'vel': snapshot['vel']}
        phase_point = system_from_snapshot(system, snapshot)
        engine.add_to_path(
            path,
            phase_point,
            left,
            right
        )
    _load_energies_for_path(path, dirname)
    return _check_path(path, path_ensemble)


def _load_external_trajectory(dirname, engine, copy=True):
    """Load an external trajectory.

    Here, we also do some moving of files to set up for a path
    simulation (if copy=True).

    Parameters
    ----------
    dirname : string
        The directory where we can find the trajectory file(s).
    engine : object like :py:class:`.ExternalMDEngine`
        The engine we use, here it is used to access the directories
        for the new simulation.
    copy : boolean, opt
        If true, it copies the files to the main ensemble folders.

    Returns
    -------
    traj : dict
        A dictionary containing the trajectory information. Here,
        the trajectory information is a list of files with indices and
        information about velocity direction.

    """
    traj_file_name = os.path.join(dirname, 'traj.txt')
    with PathExtFile(traj_file_name, 'r') as trajfile:
        # Just get the first trajectory:
        traj = next(trajfile.load())
        # Hard-copy the files. Here we assume that they are stored
        # in an folder called "accepted" and that we can get the
        # **correct** names from the traj.txt file!
        files = set([])
        for snapshot in traj['data']:
            filename = os.path.join(os.path.abspath(dirname),
                                    'accepted', snapshot[1])
            files.add(filename)
        print_to_screen('Copying trajectory files.')

        # If the files are not copied, the location has to be consistent (loc)
        loc = dirname
        if copy:
            loc = engine.exe_dir
            for filename in files:
                logger.debug('Copying %s -> %s', filename, engine.exe_dir)
                shutil.copy(filename, engine.exe_dir)
        # Update trajectory to use full path names:
        for i, snapshot in enumerate(traj['data']):
            config = os.path.join(loc, snapshot[1])
            traj['data'][i][1] = config
            reverse = int(snapshot[3]) == -1
            idx = int(snapshot[2])
            traj['data'][i][2] = idx
            traj['data'][i][3] = reverse
        return traj


def read_path_files_ext(path, dirname, ensemble):
    """Read data needed for a path from a directory.

    Parameters
    ----------
    path : object like :py:class:`.Path`
        The path we are to set up/fill.
    path_ensemble : object like :py:class:`.PathEnsembleExt`
        The path ensemble the path could be added to.
    dirname : string
        The path to the directory with the input files.
    ensemble : dictionary of objects
        It contains:

        * `path_ensemble`: object like :py:class:`.PathEnsemble`
          This is the path ensemble to perform the initialization.
        * `order_function`: object like :py:class:`.OrderParameter`
          The class used for obtaining the order parameter(s).
        * `engine`: object like :py:class:`.EngineBase`
          The engine to use for propagating a path.
        * `system`: object like :py:class:`.System`
          A system object we can use when calculating the order parameter(s).

    """
    path_ensemble = ensemble['path_ensemble']
    engine = ensemble['engine']
    system = ensemble['system']
    order_function = ensemble['order_function']
    left, _, right = path_ensemble.interfaces
    traj = _load_external_trajectory(dirname, engine)
    orderdata = _load_order_parameters_ext(traj, dirname,
                                           order_function, engine)
    # Add to path:
    print_to_screen('Creating path from files.')
    logger.debug('Creating path from files.')
    for snapshot, order in zip(traj['data'], orderdata):
        snapshot = {'order': order,
                    'pos': (snapshot[1], snapshot[2]),
                    'vel': snapshot[3],
                    'vpot': None,
                    'ekin': None}
        phase_point = system_from_snapshot(system, snapshot)
        engine.add_to_path(
            path,
            phase_point,
            left,
            right
        )
    _load_energies_for_path(path, dirname)
    return _check_path(path, path_ensemble)


def system_from_snapshot(system, snapshot):
    """Create a system from a given snapshot."""
    system_copy = system.copy()
    system_copy.particles.ekin = snapshot.get('ekin', None)
    system_copy.particles.vpot = snapshot.get('vpot', None)
    system_copy.order = snapshot.get('order', None)
    system_copy.particles.set_pos(snapshot.get('pos', None))
    system_copy.particles.set_vel(snapshot.get('vel', None))
    return system_copy

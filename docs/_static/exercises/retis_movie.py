# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This is a simple RETIS example animating the algorithm.

You can play with the interfaces, the potential parameters,
the temperature, the ratio of the different RETIS moves etc.

Have fun!
"""
import sys
import colorama
import numpy as np
import matplotlib as mpl
from matplotlib import pylab as plt
from matplotlib import animation
from matplotlib import gridspec
from pyretis.initiation import initiate_path_simulation
from pyretis.core.properties import Property
from pyretis.inout.settings import (fill_up_tis_and_retis_settings,
                                    _add_default_settings,
                                    _add_specific_default_settings)
from pyretis.setup import create_simulation
from pyretis.analysis.path_analysis import _pcross_lambda_cumulative
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


INTERFACES = [-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, 1.0]
PCROSS_LOG = False
# or: PCROSS_LOG = True
# Let us define the simulation:
SETTINGS = {}
# Basic settings for the simulation:
SETTINGS['simulation'] = {
    'task': 'retis',
    'steps': 150,
    'interfaces': INTERFACES
}
# Basic settings for the simulation box:
SETTINGS['box'] = {'cell': [3]}

# Basic settings for the system:
SETTINGS['system'] = {'units': 'reduced',
                      'temperature': 0.09,
                      'dimensions': 1}
# Basic settings for the Langevin integrator:
SETTINGS['engine'] = {
    'class': 'Langevin',
    'gamma': 0.3,
    'high_friction': False,
    'seed': 0,
    'timestep': 0.002,
}
# Particles:
SETTINGS['particles'] = {'mass': {'Ar': 1.0},
                         'name': ['Ar'],
                         'ptype': [0],
                         'position': {'generate': '1d'}}

# Potential parameters:
# The potential is: `V_\text{pot} = a x^4 - b (x - c)^2`
SETTINGS['potential'] = [
    {'a': 1.0, 'b': 2.0, 'c': 0.0, 'class': 'DoubleWell'}
]
# Settings for the order parameter:
SETTINGS['orderparameter'] = {
    'class': 'PositionVelocity',
    'dim': 'x',
    'index': 0,
    'periodic': False,
}
# TIS specific settings:
SETTINGS['tis'] = {
    'freq': 0.5,
    'maxlength': 20000,
    'aimless': True,
    'allowmaxlength': False,
    'sigma_v': -1,
    'seed': 0,
    'zero_momentum': False,
    'rescale_energy': False,
}
SETTINGS['initial-path'] = {'method': 'kick'}
# RETIS specific settings:
SETTINGS['retis'] = {
    'swapfreq': 0.5,
    'relative_shoots': None,
    'nullmoves': True,
    'swapsimul': True,
}

# For convenience:
TIMESTEP = SETTINGS['engine']['timestep']
ANALYSIS = {'ngrid': 100, 'nblock': 5}

_add_default_settings(SETTINGS)
_add_specific_default_settings(SETTINGS)
fill_up_tis_and_retis_settings(SETTINGS)

# Set up for plotting:
mpl.rc('font', size=14, family='serif')
mpl.rc('lines', color='#262626')
mpl.rc('savefig', directory=None)
NINT = len(INTERFACES)
CMAP = plt.get_cmap('Set1')
COLORS = [CMAP(float(i)/float(NINT)) for i in range(NINT)]
TXTCOLOR = {
    'SW': '#006BA4',
    'NU': '#FF800E',
    'TR': '#ABABAB',
    'SH': '#595959',
    'IN': '#808080',
}
FTOT = 0


def get_path(path):
    """Get points on a trajectory, useful for plotting.

    Parameters
    ----------
    path : object like :py:class:`.PathBase`
        The path/trajectory we are collecting points from.

    Returns
    -------
    out[0] : numpy.array
        The order parameter(s) as a function of time.
    out[1] : numpy.array
        The positions as a function of time.
    out[2] : numpy.array
        The velocities as a function of time.

    """
    order = []
    pos = []
    vel = []
    leng = path.length
    if leng <= 50:
        freq = 1
    elif 50 < leng <= 600:
        freq = 6
    elif 600 < leng < 2000:
        freq = 10
    else:
        freq = 20
    for i, point in enumerate(path.phasepoints):
        if i == 0 or i == leng-1 or i % freq == 0:
            order.append(point.order)
            pos.append(point.particles.pos[0])
            vel.append(point.particles.vel[0])
    return np.array(order), np.array(pos), np.array(vel)


def step_txt(ensembles, retis_result, prun):
    """Return some text information about the RETIS step.

    Parameters
    ----------
    ensembles : list of dicts, it contains:
       *  path_ensemble: objects like :py:class:`.PathEnsemble`
               The different path ensembles we are simulating.
    retis_result : list of lists
        The results from a RETIS simulation step.
    prun : list of floats
        The running average for the ensemble crossing probabilities.

    Returns
    -------
    out : list of strings
        Each string contains information about the move performed
        in an ensemble.

    """
    txt = []
    force = 0  # counter for force evaluations
    # The force counter will just count the number of
    # force evaluations for generating a path, it will not include
    # the force evaluation that might be performed prior to propagation,
    # i.e. the initial force evaluation since this can be stored in memory,
    # along the path -> i.e. we can avoid it.
    for ensemble in ensembles:
        path_ensemble = ensemble['path_ensemble']
        name = path_ensemble.ensemble_name
        idx = path_ensemble.ensemble_number
        name_of_move = retis_result['move-{}'.format(idx)]
        accepted = retis_result['accept-{}'.format(idx)]
        line = []
        if name_of_move == 'swap':
            idx2 = retis_result['all-{}'.format(idx)]['swap-with']
            name2 = ensembles[idx2]['path_ensemble'].ensemble_name
            move = '{} {},'.format(name_of_move, name2)
            if idx == 0 or (idx == 1 and idx2 == 0):
                # Evaluate forces when swapping [0^-] <-> [0^+]:
                force += path_ensemble.paths[-1]['length'] - 2
        elif name_of_move == 'tis':
            trial_path = retis_result['path-{}'.format(idx)]
            tis_move = trial_path.generated[0]
            move = '{} ({}),'.format(name_of_move, tis_move)
            if tis_move == 'sh':
                force += path_ensemble.paths[-1]['length'] - 1
        else:
            move = '{},'.format(name_of_move)
        line.append('{}: {:11s}'.format(name, move))
        line.append('{},'.format(accepted))
        if idx > 0:
            line.append('p = {:<8.6g}'.format(prun[idx]))
        txt.append(' '.join(line))
    return txt, force


def probability_path_ensemble(path_ensemble, step, prun, orderp):
    """Update running estimates of probabilities for an ensemble.

    Parameters
    ----------
    path_ensemble : objects like :py:class:`.PathEnsemble`
        The ensemble to analyse.
    step : int
        The current simulation step.
    prun : float
        Current running estimate for the crossing probability
        for this ensemble.
    orderp : numpy.array
        The maximum order parameters for all accepted paths.

    Returns
    -------
    prun : float
        Updated running average for the crossing probability.
    ordernew : numpy.array
        Updated list of all order parameters seen.
    lamb : numpy.array
        Coordinates used for obtaining the crossing probability
        as a function of the order parameter.
    pcross : numpy.array
        The crossing probability as a function of the order
        parameter.

    """
    ordermax = path_ensemble.last_path.ordermax[0]
    success = 1 if ordermax > path_ensemble.detect else 0
    if prun is None:
        prun = success
    else:
        npath = step + 1
        prun = float(success + prun * (npath - 1)) / float(npath)
    if orderp is None:
        ordernew = np.array([ordermax])
    else:
        ordernew = np.append(orderp, ordermax)
    ordermax = min(max(ordernew), max(path_ensemble.interfaces))
    ordermin = path_ensemble.interfaces[1]
    pcross, lamb = _pcross_lambda_cumulative(ordernew, ordermin, ordermax,
                                             ANALYSIS['ngrid'])
    return prun, ordernew, lamb, pcross


def simple_block(data, nblock):
    """Block error analysis for running average."""
    length = len(data)
    if length < nblock:
        return float('inf')
    skip, _ = divmod(len(data), nblock)
    run_avg = []
    block_avg = []
    for i in range(nblock):
        idx = (i + 1) * skip - 1
        run_avg.append(data[idx])
        if i == 0:
            block_avg.append(run_avg[i])
        else:
            block_avg.append((i + 1)*run_avg[i] - i*run_avg[i-1])
    var = sum([(x - data[-1])**2 for x in block_avg]) / (nblock - 1)
    err = np.sqrt(var / nblock)
    return err


def analyse_path_ensembles(ensembles, step, variables):
    """Calculate some properties for the path ensembles.

    Here we obtain the crossing probabilities and the initial flux.

    Parameters
    ----------
    ensembles : a list of dicts. They contain:
      *  objects like :py:class:`.PathEnsemble`
             These objects contain information about accepted paths for the
             different path ensembles.
    step : integer
        The current step number.
    variables : dict of objects
        This dict contains variables we can use for updating derived data
        for the different path ensembles.

    Returns
    -------
    out : dict of floats
        A dictionary with computed initial flux, crossing probability and
        errors in these.

    """
    calc = {'flux': 0, 'fluxe': 0,
            'pcross': 0, 'pcrosse': float('inf'),
            'kab': 0, 'kabe': float('inf')}
    length0 = variables['length0']
    length1 = variables['length1']
    length0.add(ensembles[0]['path_ensemble'].last_path.length)
    length1.add(ensembles[1]['path_ensemble'].last_path.length)
    flux = 1.0 / ((length0.mean + length1.mean - 4.0) * TIMESTEP)
    calc['flux'] = flux
    variables['flux'].append(flux)
    fluxe = simple_block(variables['flux'], ANALYSIS['nblock'])
    calc['fluxe'] = fluxe
    accprob = 1.0
    matched_prob = []
    matched_lamb = []
    for i, ensemble in enumerate(ensembles):
        if i == 0:
            continue
        prun = variables['prun'][i]
        orderp = variables['orderp'][i]
        prun, ordernew, lamb, pcross = probability_path_ensemble(
            ensemble['path_ensemble'], step, prun, orderp)
        variables['orderp'][i] = ordernew
        variables['prun'][i] = prun
        variables['pcross'][i] = pcross
        variables['lamb'][i] = lamb
        idx = np.where(lamb <= ensemble['path_ensemble'].detect)[0]
        matched_lamb.extend(lamb[idx])
        matched_prob.extend(pcross[idx] * accprob)
        variables['pcross2'][i] = pcross*accprob
        accprob *= prun
    variables['mlamb'] = matched_lamb
    variables['mpcross'] = matched_prob
    if matched_lamb[-1] < INTERFACES[-1]:
        pcross = 0
        pcrosse = float('inf')
    else:
        pcross = matched_prob[-1]
        variables['match'].append(pcross)
        pcrosse = simple_block(variables['match'], ANALYSIS['nblock'])
    calc['pcross'] = pcross
    calc['pcrosse'] = pcrosse
    calc['kab'] = flux * pcross
    if pcross == 0 or flux == 0:
        calc['kabe'] = float('inf')
    else:
        calc['kabe'] = calc['kab'] * np.sqrt(pcrosse**2 / pcross**2 +
                                             fluxe**2 / flux**2)
    return calc


def matplotlib_setup():
    """Set up matplotlib & prepare animation."""
    new_fig = plt.figure(figsize=(12, 6))
    grid = gridspec.GridSpec(3, 4)
    ax1 = new_fig.add_subplot(grid[:, :2])
    ax2 = new_fig.add_subplot(grid[0, 2:])
    ax3 = new_fig.add_subplot(grid[1:, 2:])
    if not PCROSS_LOG:
        axp = ax3.twinx()
    axes = (ax1, ax2, ax3)
    plot_patches = {'paths': [],
                    'prob': [],
                    'prob2': [],
                    'matched': None,
                    'fluxline': None,
                    'txtmove': [],
                    'txtcycle': None,
                    'start': [],
                    'end': []}
    k = 0
    for i, pos in enumerate(INTERFACES):
        ax1.axvline(x=pos, lw=2, ls=':', color='#262626')
        ax3.axvline(x=pos, lw=2, ls=':', color='#262626')
        newline, = ax1.plot([], [], lw=3, ls='-', color=COLORS[i])
        plot_patches['paths'].append(newline)
        if PCROSS_LOG:
            newlinep2, = ax3.plot([], [], lw=3, ls='-',
                                  color=newline.get_color())
            plot_patches['prob2'].append(newlinep2)
        else:
            newlinep, = ax3.plot([], [], lw=3, ls='-',
                                 color=newline.get_color())
            plot_patches['prob'].append(newlinep)

        newscat = ax1.scatter(None, None, s=75, marker='o',
                              color=newline.get_color())
        plot_patches['start'].append(newscat)
        newscat = ax1.scatter(None, None, s=75, marker='s',
                              color=newline.get_color())
        plot_patches['end'].append(newscat)
        # Note, we also add a line for [0^-], this is just to get
        # the colors right, plines[0] is not used for anything else...
        if i == 0:
            txt_x = pos-0.4
            txt_y = 1.8
        else:
            txt_x = INTERFACES[i-1]
            txt_y = 1.8 - (k-1) * 0.2
        txtbox = ax1.text(txt_x, txt_y, '', transform=ax1.transData,
                          backgroundcolor='w', fontsize=12)
        plot_patches['txtmove'].append(txtbox)
        if k % 3 == 0 and k > 0:
            k = 0
        k += 1
    plot_patches['txtcycle'] = ax1.text(-1.4, -1.8, '',
                                        transform=ax1.transData,
                                        backgroundcolor='w', fontsize=14)
    ax1.set_xlabel(r'Order parameter ($\lambda$)')
    ax1.set_ylabel(r'Velocity ($\dot{\lambda}$)')
    if PCROSS_LOG:
        ax3.set_yscale('log')
        plot_patches['matched'] = ax3.plot([], [], lw=6, ls='-',
                                           color='#262626',
                                           zorder=0, alpha=0.7)[0]
        ax3.set_xlim(-1, 1.05)
        ax3.set_ylim(1e-7, 1)
    else:
        axp.set_yscale('log')
        plot_patches['matched'] = axp.plot([], [], lw=6, ls='-',
                                           color='#262626',
                                           zorder=0, alpha=0.7)[0]
        ax3.set_xlim(-1, 1.05)
        axp.set_ylim(1e-7, 1)
        ax3.set_ylim(0, 1)
    ax3.set_xlabel(r'$\lambda$')
    ax3.set_ylabel(r'Probability')

    ax1.set_ylim(-2, 2)
    ax1.set_xlim(-1.5, 1.5)

    plot_patches['fluxline'] = ax2.plot([0], [0], lw=3, ls='-',
                                        color='#4C72B0')[0]
    ax2.set_ylim(0, 1)
    ax2.set_xlim(0, 1)
    ax2.set_xticklabels([])
    ax2.set_ylabel('Flux')
    ax2.set_xlabel('Cycles completed')
    ax2.locator_params(axis='y', nbins=4)
    if PCROSS_LOG:
        new_fig.subplots_adjust(left=0.1, bottom=0.15, right=0.95, top=0.95,
                                wspace=0.5, hspace=0.5)
    else:
        new_fig.set_tight_layout(True)
    return new_fig, plot_patches, axes


def analyse_prob(ensemble, props, idx, step):
    """Update running estimates for probability.

    Parameters
    ----------
    ensemble : object like :py:class:`.PathEnsemble`
        The ensemble to analyse.
    props : dict
        A dictionary for storing properties we calculate.
    idx : int
        An index for the path ensemble.
    step : int
        The current simulation step.

    """
    orderp = ensemble.last_path.ordermax[0]
    success = 1 if orderp > ensemble.detect else 0
    prun = props['prun'][idx]
    if not prun:
        prun.append(success)
    else:
        npath = step + 1
        prun.append(float(success + prun[-1] * (npath-1)) / float(npath))
    props['orderp'][idx].append(orderp)
    orderparam = np.array(props['orderp'][idx])
    ordermax = min(orderparam.max(), max(ensemble.interfaces))
    ordermin = ensemble.interfaces[1]
    pcross, lamb = _pcross_lambda_cumulative(orderparam, ordermin, ordermax,
                                             ANALYSIS['ngrid'])
    props['pcross'][idx] = (lamb, pcross)


def update(frame, simulation, plot_patches, variables, axes):
    """Run the simulation and update plots.

    Parameters
    ----------
    frame : int
        The current frame number, supplied by animation.FuncAnimation.
    simulation : object like :py:class:`.Simulation`
        The simulation we are running.
    plot_patches : dict of objects
        This dict contains the lines, text boxes etc. from matplotlib
        which we will use to display our data.
    variables : dict of objects
        A dict with results from the simulation.
    axes : tuple
        This tuple contains the axes used for plotting.

    Returns
    -------
    out : list
        list of the patches to be drawn.

    """
    patches = []
    if not simulation.is_finished():
        result = simulation.step()
        step = result['cycle']['step']
        ensembles = simulation.ensembles
        logger.info('# Current cycle: {}'.format(step))
        anr = analyse_path_ensembles(ensembles, step, variables)
        retis_txt, force = step_txt(ensembles, result, variables['prun'])
        global FTOT
        FTOT += force
        for line in retis_txt:
            logger.info('# {}'.format(line))
        logger.info('# Flux: {flux:<8.6g} +- {fluxe:<8.6g}'.format(**anr))
        logger.info(('# Crossing probability: {pcross:<8.6g} +-'
                     '{pcrosse:<8.6g}').format(**anr))
        logger.info('# K_AB: {kab:<8.6g} +- {kabe:<8.6g}'.format(**anr))
        logger.info('# No. of force evaluations: {:g}'.format(force))

        for i, ensemble in enumerate(ensembles):
            _, pos, vel = get_path(ensemble['path_ensemble'].last_path)
            plot_patches['paths'][i].set_data(pos, vel)
            patches.append(plot_patches['paths'][i])
            plot_patches['start'][i].set_offsets((pos[0][0], vel[0][0]))
            plot_patches['start'][i].set_visible(True)
            patches.append(plot_patches['start'][i])
            plot_patches['end'][i].set_offsets((pos[-1][0], vel[-1][0]))
            plot_patches['end'][i].set_visible(True)
            patches.append(plot_patches['end'][i])
            if i > 0:
                if not PCROSS_LOG:
                    plot_patches['prob'][i].set_data(variables['lamb'][i],
                                                     variables['pcross'][i])
                    patches.append(plot_patches['prob'][i])

        fluxx, fluxy = plot_patches['fluxline'].get_data()
        fluxx = np.append(fluxx, fluxx[-1] + 1)
        fluxy = np.append(fluxy, anr['flux'])
        plot_patches['fluxline'].set_data(fluxx, fluxy)
        axes[1].set_xlim(0, step+1)
        axes[1].set_ylim(0, 1.1*fluxy.max())
        patches.append(plot_patches['fluxline'])
        for key in result:
            if key.startswith('move'):
                i = int(key.split('-')[1])
                if result[key] == 'tis':
                    trial = result['path-{}'.format(i)]
                    txtmove = trial.generated[0].upper()
                else:
                    txtmove = result[key][:2].upper()
                plot_patches['txtmove'][i].set_text(txtmove)
                plot_patches['txtmove'][i].set_color(TXTCOLOR[txtmove])
                patches.append(plot_patches['txtmove'][i])
        plot_patches['txtcycle'].set_text('Cycle: {}'.format(step))
        patches.append(plot_patches['txtcycle'])

        # match probabilities:
        for i, ensemble in enumerate(simulation.ensembles):
            if i == 0:
                continue
            if PCROSS_LOG:
                prob2 = variables['pcross2'][i]
                plot_patches['prob2'][i].set_data(variables['lamb'][i], prob2)
                patches.append(plot_patches['prob2'][i])
        plot_patches['matched'].set_data(variables['mlamb'],
                                         variables['mpcross'])
        patches.append(plot_patches['matched'])
        return patches
    # Just return without updating:
    logger.info('# Simulation is done (frame = {})'.format(frame))
    logger.info('# Total number of force evaluations: {}'.format(FTOT))
    patches = []
    for key in ('paths', 'prob', 'prob2', 'txtmove', 'start', 'end'):
        patches.extend(plot_patches[key])
    for key in ('matched', 'fluxline', 'txtcycle'):
        patches.append(plot_patches[key])
    return patches


def main():
    """Run the simulation."""
    colorama.init(autoreset=True)
    logger.info('# CREATING SYSTEM')
    simulation = create_simulation(SETTINGS)
    logger.info(simulation)
    logger.info('# GENERATING INITIAL PATHS')

    for i, _ in enumerate(initiate_path_simulation(simulation, SETTINGS)):
        path_ensemble = simulation.ensembles[i]['path_ensemble']
        name = path_ensemble.ensemble_name
        logger.info('Info about ensemble {}:'.format(name))
        logger.info(path_ensemble)
        logger.info('Info about the initial path:')
        logger.info(path_ensemble.last_path)
    # We make a dictionary of these variable for easier access:
    # Set up some variables for storing results:
    fig, plot_patches, axes = matplotlib_setup()
    variables = {'length0': Property('Path length in [0^-]'),
                 'length1': Property('Path length in [0^+]'),
                 'flux': [],
                 'match': [],
                 'orderp': [None for _ in INTERFACES],
                 'prun': [None for _ in INTERFACES],
                 'lamb': [None for _ in INTERFACES],
                 'mlamb': [None for _ in INTERFACES],
                 'pcross': [None for _ in INTERFACES],
                 'mpcross': [None for _ in INTERFACES],
                 'pcross2': [None for _ in INTERFACES]}

    _ = animation.FuncAnimation(fig, update,
                                frames=SETTINGS['simulation']['steps']+1,
                                fargs=[simulation, plot_patches, variables,
                                       axes],
                                repeat=False,
                                interval=10,
                                blit=use_blit())
    plt.show()


def use_blit():
    """Check if we should use blitting."""
    return not sys.platform.startswith('darwin')


if __name__ == '__main__':
    main()

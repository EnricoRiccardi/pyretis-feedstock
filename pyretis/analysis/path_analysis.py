# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Methods for analysis of path ensembles.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

analyse_path_ensemble (:py:func:`.analyse_path_ensemble`)
    Method to analyse a path ensemble, it will calculate crossing
    probabilities and information about moves etc. This method
    can be applied to files as well as path ensemble objects.

analyse_path_ensemble_object (:py:func:`.analyse_path_ensemble_object`)
    Method to analyse a path ensemble, it will calculate crossing
    probabilities and information about moves etc. This method is
    intended to work directly on path ensemble objects.

match_probabilities (:py:func:`.match_probabilities`)
    Match probabilities from several path ensembles and calculate
    efficiencies and the error for the matched probability.

retis_flux (:py:func:`.retis_flux`)
    Calculate the initial flux with errors for a RETIS simulation.

retis_rate (:py:func:`.retis_rate`)
    Calculate the rate constant with errors for a RETIS simulation.
"""
import os
import logging
import numpy as np
from pyretis.analysis.analysis import running_average, block_error_corr
from pyretis.analysis.histogram import histogram, histogram_and_avg
from pyretis.inout.formats import OrderPathFile

logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())

np.set_printoptions(legacy='1.21')

__all__ = ['analyse_path_ensemble', 'analyse_repptis_ensemble',
           'match_probabilities', 'retis_flux', 'retis_rate',
           'perm_calculations', 'skip_paths']


SHOOTING_MOVES = {'sh', 'ss', 'wt', 'wf', 'ts'}
INITIAL_MOVES = {'ki', 'ld', 'is', 're'}


def _pcross_lambda_cumulative(orderparam, ordermin, ordermax, ngrid,
                              weights=None, ha_weights=None):
    """Obtain crossing probability as a function of the order parameter.

    It will do the actual calculation of the crossing probability as
    a function of order parameter. It is split off from `pcross_lambda`
    since the analysis is intended to be backward compatible with the
    output/results from the old ``TISMOL FORTRAN`` program.

    Parameters
    ----------
    orderparam : numpy.array
        Array containing the order parameters.
    ordermin : float
        Minimum allowed order parameter.
    ordermax : float
        Maximum allowed order parameter.
    ngrid : int
        This is the number of grid points.
    weights : numpy.array, optional
        The weight of each order parameter. This is used in order to
        count a specific order parameter more than once. If not given,
        the values in `orderparam` will be weighted equally.
    ha_weights : numpy.array, optional
        The weights of the high acceptance scheme (e.g. Stone Skipping)

    Returns
    -------
    lamb: numpy.array
        The grid on which pcross is based.
    pcross : numpy.array
        Estimated crossing probability distribution.

    """
    lamb = np.linspace(ordermin, ordermax, ngrid)
    pcross = np.zeros(ngrid)
    delta_l = lamb[1] - lamb[0]
    sumw = 0.0
    for i, orderp in enumerate(orderparam):
        idx = np.floor((orderp - ordermin) / delta_l)
        idx = max(0, int(idx) + 1)
        # +1: idx is here defined so that lamb[idx-1] <= orderp < lamb[idx]
        # further this lambda will contribute up to and including lamb[idx]
        # this is accomplished by the idx+1 when summing weights below
        weight = 1. if weights is None else weights[i]
        ha_weight = 1. if ha_weights is None else ha_weights[i]

        sumw += weight*ha_weight
        if idx >= ngrid:
            pcross += weight*ha_weight
        else:
            pcross[:idx + 1] += weight*ha_weight  # +1 to include up to idx
    pcross /= sumw  # normalisation
    return lamb, pcross


def _get_path_length(path, ensemble_number):
    """Return the path length for different moves.

    Different moves may have a different way of obtaining the path
    length. (Example: time-reversal vs. shooting move).

    Parameters
    ----------
    path : dict
        This is the dict containing the information about the path.
        It can typically be obtained by iterating over the path
        ensemble object, e.g. with a
        `for path in path_ensemble.get_paths():`.
    ensemble_number : int
        This integer identifies the ensemble. This is used for
        the swapping moves in [0^-] and [0^+].

    Returns
    -------
    out : int
        The path length

    """
    move = path['generated'][0]
    return_table = {'tr': 0, 's+': 0, 's-': 0, '00': 0, 'mr': 0}
    if move in return_table:
        if move == 's+' and ensemble_number == 0:
            return path['length'] - 2
        if move == 's-' and ensemble_number == 1:
            return path['length'] - 2
        return return_table[move]
    if move in SHOOTING_MOVES:
        return path['length'] - 1
    if move in INITIAL_MOVES:
        logger.info('Skipped initial path (move "%s")', move)
        return None
    logger.warning('Skipped path with unknown mc move: %s', move)
    return None


def _update_shoot_stats(shoot_stats, path):
    """Update the shooting statistics with the status of the given path.

    Parameters
    ----------
    shoot_stats : dict
        This dict contains the results from the analysis of shooting
        moves. E.g. `shoot_stats[key]` contain the order parameters
        for the status `key` which can be the different statuses
        defined in `pyretis.core.path._STATUS` or 'REJ' (for rejected).
    path : dict
        This is the path information, represented as a dictionary.

    Returns
    -------
    out : None
        Returns `None` but will update `shoot_stats` for shooting moves.

    """
    move = path['generated'][0]
    if move in SHOOTING_MOVES:
        orderp = path['generated'][1]
        status = path['status']
        if status not in shoot_stats:
            shoot_stats[status] = []
        shoot_stats[status].append(orderp)
        if status != 'ACC':
            shoot_stats['REJ'].append(orderp)
        shoot_stats['ALL'].append(orderp)


def _create_shoot_histograms(shoot_stats, bins):
    """Create histograms and scale for the shoot analysis.

    Parameters
    ----------
    shoot_stats : dict
        This dict contains the results from the analysis of shooting
        moves. E.g. `shoot_stats[key]` contain the order parameters
        for the status `key` which can be the different statuses
        defined in `pyretis.core.path._STATUS` or 'REJ' (for rejected).
    bins : int
        The number of bins to use for the histograms.

    Returns
    -------
    out[0] : dict
        For each possible status ('ACC, 'BWI', etc.) this dict will
        contain a histogram as returned by the histogram function.
        It will also contain a 'REJ' key which is the concatenation of
        all rejections and a 'ALL' key which is simply all the values.
    out[1] : dict
        For each possible status ('ACC, 'BWI', etc.) this dict will
        contain the scale factors for the histograms. The scale factors
        are obtained by dividing with the 'ALL' value.

    See Also
    --------
    :py:func:`.histogram` in :py:mod:`.histogram`.

    """
    histograms = {}
    scale = {}
    for key in shoot_stats:
        if not shoot_stats[key]:
            logger.warning('No shoots data found for %s (empty histogram)',
                           key)
            shoot_stats[key] = np.array([])
            mind = 0.0
            maxd = 0.1
        else:
            shoot_stats[key] = np.array(shoot_stats[key])
            mind = shoot_stats[key].min()
            maxd = shoot_stats[key].max()
        histograms[key] = histogram(shoot_stats[key], bins=bins,
                                    limits=(mind, maxd), density=True)
        # Guard against division by zero when there are no 'ALL' samples
        all_len = float(len(shoot_stats.get('ALL', [])))
        if all_len == 0.0:
            scale[key] = 0.0
        else:
            scale[key] = float(len(shoot_stats[key])) / all_len
    return histograms, scale


def analyse_path_ensemble(path_ensemble, settings):
    """Analyse a path ensemble.

    This function will make use of the different analysis functions and
    analyse a path ensemble. This function is more general than the
    `analyse_path_ensemble_object` function in that it should work on
    both `PathEnsemble` and `PathEnsembleFile` objects. The running
    average is updated on-the-fly, see Wikipedia for
    details [wikimov]_.

    Parameters
    ----------
    path_ensemble : object like :py:class:`.PathEnsemble`
        This is the path ensemble to analyse.
    settings : dict
        This dictionary contains settings for the analysis.
        We make use of the following keys:

        * `ngrid`: The number of grid points for calculating the
          crossing probability as a function of the order parameter.
        * `maxblock`: The max length of the blocks for the block error
          analysis. Note that this will maximum be equal the half of the
          length of the data, see `block_error` in `.analysis`.
        * `blockskip`: Can be used to skip certain block lengths.
          A `blockskip` equal to `n` will consider every n'th block up
          to `maxblock`, i.e. it will use block lengths equal to `1`,
          `1+n`, `1+2n`, etc.
        * `bins`: The number of bins to use for creating histograms.

    Returns
    -------
    out : dict
        This dictionary contains the main results for the analysis which
        can be used for plotting or other kinds of output.

    See Also
    --------
    :py:func:`._update_shoot_stats`, :py:func:`.pcross_lambda_cumulative`
    and :py:func:`._create_shoot_histograms`.

    References
    ----------
    .. [wikimov] Wikipedia, "Moving Average",
       https://en.wikipedia.org/wiki/Moving_average

    """
    detect = settings['tis']['detect']
    ens_number = path_ensemble.ensemble_number
    if ens_number == 0:
        return analyse_path_ensemble0(path_ensemble, settings)
    result = {'cycle': [],
              'detect': detect,
              'ensemble': path_ensemble.ensemble_name,
              'ensembleid': ens_number,
              'interfaces': list(path_ensemble.interfaces)}
    orderparam = []  # list of all accepted order parameters
    weights = []
    pdata = []
    length_acc = []
    length_all = []
    shoot_stats = {'REJ': [], 'ALL': []}
    ha_w = 0
    nacc = 0
    nacc_swap = 0
    npath = 0
    npath_swap = 0
    skip_lines = settings.get("analysis", {}).get("skip", 0)
    skip_weight = 0
    production_run = False  # True from the first ACC path not in I.M.
    for i, path in enumerate(path_ensemble.get_paths()):  # loop over all paths
        if not production_run and path['status'] == 'ACC' and \
                path['generated'][0] in SHOOTING_MOVES:
            production_run = True
        if not production_run:
            continue

        npath += 1
        if path['generated'][0] in ('s-', 's+'):
            npath_swap += 1
            if path['status'] == 'ACC':
                nacc_swap += 1
        if path['status'] == 'ACC':
            nacc += 1
            weights.append(1)
            ha_w = 1./path.get('weight', 1.)
            orderparam.append(path['ordermax'][0])
            length_acc.append([path['length'], ha_w])
            if ens_number != 0:
                success = 1 if path['ordermax'][0] > detect else 0.
                pdata.append([success, ha_w])  # Store data for block analysis
        elif nacc != 0:  # just increase the weights
            weights[-1] += 1
        if i < skip_lines:
            skip_weight += 1

        # we also update the running average of the probability here:
        result['cycle'].append(path['cycle'])
        # get the length - note that this length depends on the type of move
        # see the `_get_path_length` function.
        length = _get_path_length(path, ens_number)
        if length is not None:
            length_all.append([length, ha_w])
        # update the shoot stats, this will only be done for shooting moves
        _update_shoot_stats(shoot_stats, path)

    # Update the weights for skipped paths
    weights, nacc = skip_paths(weights, nacc, skip_weight)
    # Check if we have enough data to do analysis:
    if nacc == 0:
        logger.warning('No accepted paths to analyse in ensemble %s',
                       ens_number)
        return result
    # When restarting a simulations, or by stacking together different
    # simulations, the cycles might not be in order. We thus reset the counter.
    result['cycle'] = np.arange(len(result['cycle']))
    # 1) result['prun'] is already calculated.
    pdata = np.array(pdata)
    result['pdata'] = np.array(pdata)
    # 2) lambda pcross:
    analysis = settings['analysis']
    if ens_number != 0:
        # 1) result['prun'] is calculated here with the correct path weights.
        result['prun'] = running_average(np.repeat(pdata, weights, axis=0))
        result['pdata'] = np.array(pdata)
        # 2) lambda pcross:
        orderparam = np.array(orderparam)
        ordermax = min(orderparam.max(), max(path_ensemble.interfaces))
        result['pcross'] = _pcross_lambda_cumulative(
            orderparam,
            path_ensemble.interfaces[1],
            ordermax,
            analysis['ngrid'],
            weights=weights,
            ha_weights=result['pdata'][:, 1])

        # 3) block error analysis:
        result['blockerror'] = block_error_corr(
                data=np.repeat(pdata, weights, axis=0),
                maxblock=analysis['maxblock'],
                blockskip=analysis['blockskip'])
    # 4) length analysis:
    length_acc = np.array(length_acc)
    length_all = np.array(length_all)
    hist1 = histogram_and_avg(np.repeat(length_acc[:, 0], weights),
                              analysis['bins'], density=True,
                              weights=np.repeat(length_acc[:, 1], weights))
    hist2 = histogram_and_avg(np.array(length_all[:, 0]),
                              analysis['bins'], density=True,
                              weights=length_all[:, 1])
    result['pathlength'] = (hist1, hist2)
    # 5) shoots analysis:
    result['shoots'] = _create_shoot_histograms(shoot_stats,
                                                analysis['bins'])
    # 6) Add some simple efficiency metrics:
    result['tis-cycles'] = npath
    result['efficiency'] = [float(nacc - nacc_swap) /
                            float(npath - npath_swap)]
    if npath_swap == 0:
        result['efficiency'].append(np.nan)
    else:
        result['efficiency'].append(float(nacc_swap) / float(npath_swap))

    # extra analysis for the [0^- and 0^+] ensembles in case we will determine
    # the initial flux:
    if ens_number in [0, 1]:
        result['lengtherror'] = block_error_corr(
                data=np.repeat(length_acc, weights, axis=0),
                maxblock=analysis['maxblock'],
                blockskip=analysis['blockskip'])
        lenge2 = result['lengtherror'][4] * hist1[2][0] / (hist1[2][0]-2.)
        result['fluxlength'] = [hist1[2][0]-2.0, lenge2,
                                lenge2 * (hist1[2][0]-2.)]
        result['fluxlength'].append(result['efficiency'][1] * lenge2**2)
    return result


def analyse_path_ensemble0(path_ensemble, settings):
    """Analyse the [0^-] ensemble.

    Parameters
    ----------
    path_ensemble : object like :py:class:`.PathEnsemble`
        This is the path ensemble to analyse.
    settings : dict
        This dictionary contains settings for the analysis.
        We make use of the following keys:

        * `ngrid`: The number of grid points for calculating the
          crossing probability as a function of the order parameter.
        * `maxblock`: The max length of the blocks for the block error
          analysis. Note that this will maximum be equal the half of the
          length of the data, see `block_error` in `.analysis`.
        * `blockskip`: Can be used to skip certain block lengths.
          A `blockskip` equal to `n` will consider every n'th block up
          to `maxblock`, i.e. it will use block lengths equal to `1`,
          `1+n`, `1+2n`, etc.
        * `bins`: The number of bins to use for creating histograms.

    Returns
    -------
    result : dict
        The results from the analysis for the ensemble.

    """
    simulation_settings = settings.get('simulation', {})
    permeability = simulation_settings.get('permeability', False)
    detect = settings['tis']['detect']
    ensemble_number = path_ensemble.ensemble_number
    sub_cycles = settings.get('engine', {}).get('subcycles', 1)
    timestep = settings.get('engine', {}).get('timestep', 1)*sub_cycles
    if permeability:
        path_ensemble.interfaces = (settings['simulation'].get('zero_left',
                                                               float('-inf')
                                                               ),
                                    path_ensemble.interfaces[1],
                                    path_ensemble.interfaces[2])

    result = {'cycle': [],
              'detect': detect,
              'ensemble': path_ensemble.ensemble_name,
              'ensembleid': ensemble_number,
              'interfaces': list(path_ensemble.interfaces),
              'permeability': permeability}
    length_acc, length_all, weights = [], [], []
    # Regular setting[analysis] is overwritten by the code that generates
    # out.rst, but this one is complete
    analysis = settings['analysis']

    # Permeability analysis
    bin_edges = analysis.get('tau_ref_bin', None)
    end_acc, tau_acc, tau_ref_acc = [], [], []
    tau_ref_bins = make_tau_ref_bins(list(path_ensemble.interfaces), nbins=10,
                                     tau_region=bin_edges)
    shoot_stats = {'REJ': [], 'ALL': []}
    nacc, npath = 0, 0
    nacc_swap, npath_swap = 0, 0
    production_run = True
    fn_op = path_ensemble.filename.rsplit("/", 1)[0]+"/order.txt"
    o_ps = None
    op_data = None
    skip_lines = settings.get("analysis", {}).get("skip", 0)
    skip_weight = 0

    if os.path.isfile(fn_op):
        op_file = OrderPathFile(fn_op, 'r')
        o_ps = op_file.load()
        op_file.close()
    for num, path in enumerate(path_ensemble.get_paths()):  # loop over paths
        if o_ps is not None:
            err = False
            try:
                o_p = next(o_ps)
            except StopIteration:
                o_ps = None
                err = True
            if not err:
                cycle_comment = [i for i in o_p['comment'] if "Cycle" in i][-1]
                op_cycle = int(
                    cycle_comment.split(',')[0].split(':')[1].strip())
                err = op_cycle != path['cycle']
                if not err:
                    op_data = o_p['data'].T
            if err:
                op_data = None
                msg = ("Order-file does not align with path-ensemble setting "
                       "op data to None")
                logger.warning(msg)

        if path['generated'][0] in INITIAL_MOVES:
            production_run = False
        if not production_run and path['status'] == 'ACC' and \
                path['generated'][0] in SHOOTING_MOVES:
            production_run = True
        if not production_run:
            continue
        npath += 1
        if path['generated'][0] in ('s-', 's+'):
            npath_swap += 1
            if path['status'] == 'ACC':
                nacc_swap += 1
        if path['status'] == 'ACC':
            nacc += 1
            ha_w = 1./path.get('weight', 1.)
            weights.append(1)

            length_acc.append([path['length'], ha_w])
            # Get the end point
            end_acc.append(path['interface'][-1])
            # Calculate tau
            tau_acc.append(_calc_tau(op_data, bin_edges, timestep=timestep))
            tau_ref_acc.append(
                [_calc_tau(op_data, edges, timestep=timestep)
                 for edges in tau_ref_bins]
            )
        elif nacc != 0:  # just increase the weights
            weights[-1] += 1
        if num < skip_lines:
            skip_weight += 1
        result['cycle'].append(path['cycle'])
        length = _get_path_length(path, ensemble_number)
        if length is not None:
            ha_w = 1./path.get('weight', 1.)
            length_all.append([length, ha_w])
        # update the shoot stats, this will only be done for shooting moves
        _update_shoot_stats(shoot_stats, path)
    # Update the weights to skip properly
    weights, nacc = skip_paths(weights, nacc, skip_weight)
    length_acc = np.array(length_acc)
    length_all = np.array(length_all)

    # Check if we have enough data to do analysis:
    ens_number = 0
    if nacc == 0:
        logger.warning('No accepted paths to analyse in ensemble %s',
                       ens_number)
        return result

    # Perform the different analysis tasks:

    result['cycle'] = np.array(result['cycle'])
    # 1) length analysis:
    hist1 = histogram_and_avg(np.repeat(length_acc[:, 0], weights),
                              analysis['bins'], density=True)
    hist2 = histogram_and_avg(np.array(length_all[:, 0]),
                              analysis['bins'], density=True)
    result['pathlength'] = (hist1, hist2)
    # 2) block error of lengths:
    result['lengtherror'] = block_error_corr(data=np.repeat(length_acc,
                                                            weights, axis=0),
                                             maxblock=analysis['maxblock'],
                                             blockskip=analysis['blockskip'])
    # 3) shoots analysis:
    result['shoots'] = _create_shoot_histograms(shoot_stats,
                                                analysis['bins'])
    # 4) Add some simple efficiency metrics:
    result['efficiency'] = [float(nacc - nacc_swap) /
                            float(npath - npath_swap)]
    if npath_swap == 0:
        result['efficiency'].append(np.nan)
    else:
        result['efficiency'].append(float(nacc_swap) / float(npath_swap))
    lenge2 = result['lengtherror'][4] * hist1[2][0] / (hist1[2][0]-2.)
    result['fluxlength'] = [hist1[2][0]-2.0, lenge2, lenge2 * (hist1[2][0]-2.)]
    result['fluxlength'].append(result['efficiency'][1] * lenge2**2)
    result['tis-cycles'] = npath
    # 5) optional, do the permeability analysis
    if permeability:
        end_acc = [1 if i == 'R' else 0 for i in end_acc]
        # Xi is 1/running average
        result['xi'] = running_average(np.repeat(end_acc, weights))
        result['xierror'] = block_error_corr(result['xi'],
                                             maxblock=analysis['maxblock'],
                                             blockskip=analysis['blockskip'])
        result['tau'] = running_average(np.repeat(tau_acc, weights))
        result['tauerror'] = block_error_corr(result['tau'],
                                              maxblock=analysis['maxblock'],
                                              blockskip=analysis['blockskip'])
        result['tau_bin'] = bin_edges
        result['tau_ref'] = running_average(np.repeat(tau_ref_acc,
                                                      repeats=weights,
                                                      axis=0))[-1]
        result['tau_ref_bins'] = tau_ref_bins

    return result


def analyse_repptis_ensemble(path_ensemble, settings):
    """Analyse a repptis path  ensemble.

    This function is a modification of analyse_path_ensemble(), which
    make use of the different analysis functions and analyse a repptis
    path ensemble.

    Parameters
    ----------
    path_ensemble : object like :py:class:`.PathEnsemble`
        This is the path ensemble to analyse.
    settings : dict
        This dictionary contains settings for the analysis.
        We make use of the following keys:

        * `ngrid`: The number of grid points for calculating the
          crossing probability as a function of the order parameter.
        * `maxblock`: The max length of the blocks for the block error
          analysis. Note that this will maximum be equal the half of the
          length of the data, see `block_error` in `.analysis`.
        * `blockskip`: Can be used to skip certain block lengths.
          A `blockskip` equal to `n` will consider every n'th block up
          to `maxblock`, i.e. it will use block lengths equal to `1`,
          `1+n`, `1+2n`, etc.
        * `bins`: The number of bins to use for creating histograms.

    Returns
    -------
    out : dict
        This dictionary contains the main results for the analysis which
        can be used for plotting or other kinds of output.

    """
    detect = settings['tis']['detect']
    ens_number = path_ensemble.ensemble_number
    if ens_number == 0:
        return analyse_path_ensemble0(path_ensemble, settings)
    result = {'cycle': [],
              'detect': detect,
              'ensemble': path_ensemble.ensemble_name,
              'ensembleid': ens_number,
              'interfaces': list(path_ensemble.interfaces)}
    orderparam = []  # list of all accepted order parameters
    weights = []
    nacc_sl, nacc_sr = 0, 0
    weights_sl, weights_sr = [], []
    pdata = []
    pdata_sl = []  # start left
    pdata_sr = []  # start right
    length_acc = []
    length_all = []
    shoot_stats = {'REJ': [], 'ALL': []}
    ha_w = 0
    nacc = 0
    nacc_swap = 0
    npath = 0
    npath_swap = 0
    skip_lines = settings.get("analysis", {}).get("skip", 0)
    skip_weight = 0
    production_run = False  # True from the first ACC path not in I.M.
    lmrs = []  # keep track of a path's left, middle and right
    ordermin = []
    intfs000 = path_ensemble.interfaces
    ptypes = [[0, 0, 0, 0]]  # LML, LMR, RMR, RML amounts
    ptype = [0, 0, 0, 0]  # initialise before loop; updated on each ACC path
    for i, path in enumerate(path_ensemble.get_paths()):  # loop over all paths
        # we skip the first number of lines that should not be part of
        # the analysis.
        if i < skip_lines:
            continue
        if not production_run and path['status'] == 'ACC' and \
                path['generated'][0] in SHOOTING_MOVES:
            production_run = True
        if not production_run:
            continue

        npath += 1
        if path['generated'][0] in ('s-', 's+'):
            npath_swap += 1
            if path['status'] == 'ACC':
                nacc_swap += 1
        if path['status'] == 'ACC':
            nacc += 1
            weights.append(1)
            ha_w = 1./path.get('weight', 1.)
            orderparam.append(path['ordermax'][0])
            ordermin.append(path['ordermin'][0])
            length_acc.append([path['length'], ha_w])
            lmrs.append(''.join(path['interface']))
            if ens_number != 0:
                success = 1 if path['ordermax'][0] > detect else 0.
                pdata.append([success, ha_w])  # Store data for block analysis
                if lmrs[-1][0] == 'L':
                    nacc_sl += 1
                    weights_sl.append(1)
                    success_sl = 1 if path['ordermax'][0] > intfs000[2] else 0.
                    pdata_sl.append([success_sl, ha_w])
                    ptype = [0, 1, 0, 0] if success_sl == 1 else [1, 0, 0, 0]
                elif lmrs[-1][0] == 'R':
                    weights_sr.append(1)
                    nacc_sr += 1
                    success_sr = 1 if path['ordermin'][0] < intfs000[0] else 0.
                    pdata_sr.append([success_sr, ha_w])
                    ptype = [0, 0, 0, 1] if success_sr == 1 else [0, 0, 1, 0]

        elif nacc != 0:  # just increase the weights
            weights[-1] += 1
            if ens_number != 0:
                if lmrs[-1][0] == 'L':
                    weights_sl[-1] += 1
                elif lmrs[-1][0] == 'R':
                    weights_sr[-1] += 1
        if i < skip_lines:
            skip_weight += 1

        # Update ptypes with the last ptype
        ptypes.append(ptype)

        # we also update the running average of the probability here:
        result['cycle'].append(path['cycle'])
        # get the length - note that this length depends on the type of move
        # see the `_get_path_length` function.
        length = _get_path_length(path, ens_number)
        if length is not None:
            length_all.append([length, ha_w])
        # update the shoot stats, this will only be done for shooting moves
        _update_shoot_stats(shoot_stats, path)

    # save the ptypes, which will be used for recursive block analysis
    result['ptypes'] = np.cumsum(np.array(ptypes), axis=0)

    # Check if we have enough data to do analysis:
    if nacc == 0:
        logger.warning('No accepted paths to analyse in ensemble %s',
                       ens_number)
        return result
    # When restarting a simulations, or by stacking together different
    # simulations, the cycles might not be in order. We thus reset the counter.
    result['cycle'] = np.arange(len(result['cycle']))
    # 1) result['prun'] is already calculated.
    pdata = np.array(pdata)
    result['pdata'] = np.array(pdata)
    # 2) lambda pcross:
    analysis = settings['analysis']
    if ens_number != 0:
        # 1) result['prun'] is calculated here with the correct path weights.
        result['prun_sl'] = running_average(np.repeat(pdata_sl,
                                                      weights_sl,
                                                      axis=0))
        result['prun_sr'] = running_average(np.repeat(pdata_sr,
                                                      weights_sr,
                                                      axis=0))
        result['pdata'] = np.array(pdata)
        # 2) lambda pcross:
        orderparam = np.array(orderparam)
        ordermax = min(orderparam.max(), max(path_ensemble.interfaces))
        result['pcross'] = _pcross_lambda_cumulative(
            orderparam,
            path_ensemble.interfaces[1],
            ordermax,
            analysis['ngrid'],
            weights=weights,
            ha_weights=result['pdata'][:, 1])
        # 3) block error analysis:
        # if we get here, it means that either pdata_sl and/or pdata_sr is not
        # empty. So we set the non-empty one first, and then equate the empty
        # one to zeros and copy the blockerror of the non-empty one.
        if len(pdata_sl) > 0:
            result['blockerror_sl'] = block_error_corr(
                    data=np.repeat(pdata_sl, weights_sl, axis=0),
                    maxblock=analysis['maxblock'],
                    blockskip=analysis['blockskip'])
        if len(pdata_sr) > 0:
            result['blockerror_sr'] = block_error_corr(
                    data=np.repeat(pdata_sr, weights_sr, axis=0),
                    maxblock=analysis['maxblock'],
                    blockskip=analysis['blockskip'])
        if len(pdata_sl) == 0:
            result['blockerror_sl'] = result['blockerror_sr']
            result['prun_sl'] = np.zeros_like(result['prun_sr'])
        if len(pdata_sr) == 0:
            result['blockerror_sr'] = result['blockerror_sl']
            result['prun_sr'] = np.zeros_like(result['prun_sl'])

        # PPRETIS analysis
        pp_dic = {'lmrs': np.array(lmrs),
                  'interfaces': path_ensemble.interfaces,
                  'weights': np.array(weights),
                  'ordermin': np.array(ordermin),
                  'ordermax': np.array(orderparam)}
        pptup = cross_dist_distr(pp_dic)
        result['repptis'] = pptup

    # 4) length analysis:
    length_acc = np.array(length_acc)
    length_all = np.array(length_all)
    hist1 = histogram_and_avg(np.repeat(length_acc[:, 0], weights),
                              analysis['bins'], density=True,
                              weights=np.repeat(length_acc[:, 1], weights))
    hist2 = histogram_and_avg(np.array(length_all[:, 0]),
                              analysis['bins'], density=True,
                              weights=length_all[:, 1])
    result['pathlength'] = (hist1, hist2)
    # 5) shoots analysis:
    result['shoots'] = _create_shoot_histograms(shoot_stats,
                                                analysis['bins'])
    # 6) Add some simple efficiency metrics:
    result['tis-cycles'] = npath
    result['efficiency'] = [float(nacc - nacc_swap) /
                            float(npath - npath_swap)]
    swap_acc_ratio = float(nacc_swap) / float(npath_swap) if npath_swap > 0 \
        else np.nan
    result['efficiency'].append(swap_acc_ratio)
    # extra analysis for the [0^- and 0^+] ensembles in case we will determine
    # the initial flux:
    if ens_number in [0, 1]:
        result['lengtherror'] = block_error_corr(
                data=np.repeat(length_acc, weights, axis=0),
                maxblock=analysis['maxblock'],
                blockskip=analysis['blockskip'])
        lenge2 = result['lengtherror'][4] * hist1[2][0] / (hist1[2][0]-2.)
        result['fluxlength'] = [hist1[2][0]-2.0, lenge2,
                                lenge2 * (hist1[2][0]-2.)]
        result['fluxlength'].append(result['efficiency'][1] * lenge2**2)
    return result


def match_probabilities(path_results, detect, settings=None):
    """Match probabilities from several path ensembles.

    It calculates the efficiencies and errors for the matched
    probability.

    Parameters
    ----------
    path_results : list
        These are the results from the path analysis. `path_results[i]`
        contains the output from `analyse_path_ensemble` applied to
        ensemble no. `i`. Here we make use of the following keys from
        `path_results[i]`:
        * `pcross`: The crossing probability.
        * `prun`: The running average of the crossing probability.
        * `blockerror`: The output from the block error analysis.
        * `efficiency`: The output from the efficiency analysis.
    detect : list of floats
        These are the detect interfaces used in the analysis.
    settings : dict, optional
        This dictionary contains settings for the analysis.
        Here we make use of the following keys from the
        analysis section:

        * `ngrid`: The number of grid points for calculating the
          crossing probability as a function of the order parameter.
        * `maxblock`: The max length of the blocks for the block error
          analysis. Note that this will maximum be equal the half of the
          length of the data, see `block_error` in `.analysis`.
        * `blockskip`: Can be used to skip certain block lengths.
          A `blockskip` equal to `n` will consider every n'th block up
          to `maxblock`, i.e. it will use block lengths equal to `1`,
          `1+n`, `1+2n`, etc.
        * `bins`: The number of bins to use for creating histograms.

    Returns
    -------
    results : dict
        These are results for the over-all probability and error
        and also some over-all TIS efficiencies.

    """
    results = {'matched-prob': [],
               'overall-prun': [],
               'overall-err': [],
               'overall-prob': [[], []]}
    accprob = 1.0
    accprob_err = 0.0
    minlen_pdata, minlen_prun = float('inf'), float('inf')
    active_results = []
    for idet, result in zip(detect, path_results):
        if 'pcross' not in result:
            logger.warning('Skipping ensemble %s in probability matching: '
                           'no accepted paths', result.get('ensemble', '?'))
            continue
        active_results.append(result)
        # do match only in part left of idetect:
        idx = np.where(result['pcross'][0] <= idet)[0]
        results['overall-prob'][0].extend(result['pcross'][0][idx])
        results['overall-prob'][1].extend(result['pcross'][1][idx] * accprob)
        # update probabilities, error and efficiency:
        mat = np.column_stack((result['pcross'][0], result['pcross'][1]))
        mat[:, 1] *= accprob
        results['matched-prob'].append(mat)
        accprob *= result['prun'][-1]
        accprob_err += result['blockerror'][4]**2

        # Find the maximum number of cycles for ensemble
        minlen_pdata = min(minlen_pdata, len(result['pdata']))
        minlen_prun = min(minlen_prun, len(result['prun']))

    if not active_results:
        results['overall-cycle'] = np.array([])
        results['overall-prun'] = np.array([])
        results['overall-pdata'] = np.array([])
        results['overall-prob'] = np.array([[], []])
        results['prob'] = np.nan
        results['relerror'] = np.nan
        return results

    # Finally Construct the cumulative output now
    results['overall-cycle'] = np.arange(minlen_prun)
    results['overall-prun'] = [1]*minlen_prun
    results['overall-pdata'] = [[1, 1]]*minlen_pdata
    for result in active_results:
        results['overall-prun'] = np.multiply(result['prun'][-minlen_prun:],
                                              results['overall-prun'])
        results['overall-pdata'] = np.multiply(result['pdata'][-minlen_pdata:],
                                               results['overall-pdata'])

    if settings is not None:
        analysis = settings['analysis']
        results['overall-error'] = block_error_corr(results['overall-pdata'],
                                                    analysis['maxblock'],
                                                    analysis['blockskip'])

    results['overall-prob'] = np.transpose(results['overall-prob'])
    results['prob'] = accprob
    results['relerror'] = np.sqrt(accprob_err)
    return results


def retis_flux(results0, results1, timestep):
    """Calculate the initial flux for RETIS.

    Parameters
    ----------
    results0 : dict
        Results from the analysis of ensemble [0^-]
    results1 : dict
        Results from the analysis of ensemble [0^+]
    timestep : float
        The simulation timestep.

    Returns
    -------
    flux : float
        The initial flux.
    flux_error : float
        The relative error in the initial flux.

    """
    flux0 = results0['fluxlength']
    flux1 = results1['fluxlength']
    tsum = flux0[0] + flux1[0]
    flux = 1.0 / (tsum * timestep)
    flux_error = (np.sqrt((flux0[1]*flux0[0])**2 + (flux1[1]*flux1[0])**2) /
                  tsum)
    return flux, flux_error


def retis_rate(pcross, pcross_relerror, flux, flux_relerror):
    """Calculate the rate constant for RETIS.

    Parameters
    ----------
    pcross : float
        Estimated crossing probability
    pcross_relerror : float
        Relative error in crossing probability.
    flux : float
        The initial flux.
    flux_relerror : float
        Relative error in the initial flux.

    Returns
    -------
    rate : float
        The rate constant
    rate_error : float
        The relative error in the rate constant.

    """
    rate = pcross * flux
    rate_error = np.sqrt(pcross_relerror**2 + flux_relerror**2)
    return rate, rate_error


def perm_calculations(results0, pcross, pcross_err):
    """Calculate the permeability.

    This function collects the relevant permeability properties from the output
    data. It then calculates the permeability and its error based on the values
    of xi, tau, and the crossing probability pcross.

    Parameters
    ----------
    results0 : dict
        Results from the analysis of ensemble [0^-]
    pcross : float
        Estimated crossing probability
    pcross_err : float
        Relative error in crossing probability.

    Returns
    -------
    xxi : float
        The value of xi
    tau : float
        The value of tau
    perm : float
        The value of the permeability
    xi_err : float
        The error in xi
    tau_err : float
        The error in tau
    perm_err : float
        The error in permeability

    """
    out = results0['out']
    permeability = out.get("permeability", False)
    xxi, tau, perm, xi_err, tau_err, perm_err = (
        None, None, None, float('inf'), float('inf'), float('inf'))
    if permeability:
        xxi = out['xi'][-1]
        xi_err = out['xierror'][4]
        tau = out['tau'][-1]
        tau_err = out['tauerror'][4]
        perm = xxi * (1 / tau) * pcross
        perm_err = np.sqrt(xi_err**2+tau_err**2+pcross_err**2)
    return xxi, tau, perm, xi_err, tau_err, perm_err


def _calc_tau(op_data, bin_edges=None, timestep=1):
    """Calculate the value of tau.

    This function calculates tau, a measure of time spent within a certain
    interval (bin) in the [0^-'] ensemble. It does so by counting the number
    of times the order parameter falls within this interval defined by the
    bin_edges. It then multiplies this count by the timestep and divides by
    the width of the interval. If bin_edges is None or not a 2-element list,
    or if the order parameter data op_data is None, it just returns 0.

    Parameters
    ----------
    op_data : list
        The order parameter data from the [0^-'] ensemble.
    bin_edges : list, optional
        The edges of the bin to calculate tau for. If None, it returns 0.
    timestep : float, optional
        The simulation timestep. Default is 1.

    Returns
    -------
    tau : float
        The value of tau.

    """
    if bin_edges is None or len(bin_edges) != 2 or op_data is None:
        return 0
    order = op_data[1]
    bin_l, bin_r = tuple(bin_edges)
    delta = bin_r-bin_l
    count = sum(1 for i in order if bin_l < i <= bin_r)
    return count*timestep/delta


def make_tau_ref_bins(interfaces, nbins=10, tau_region=None):
    """Make the reference bins for calculating tau.

    This function makes reference bins for calculating tau based on the [0^-']
    interfaces and a number of bins. If the left interface is -inf and no
    tau region is provided, it returns a list with None. If a tua region is
    provided, it adjusts the left interface accordingly. It then divides the
    range from left to right interface into nbins equal parts and returns
    these as a list of tuples.

    Parameters
    ----------
    interfaces : list
        The interfaces from the [0^-'] ensemble.
    nbins : int, optional
        The number of bins to divide the range into. Default is 10.
    tau_region : list, optional
        The region to calculate tau for. If None, it uses the whole range.

    Returns
    -------
    ref_bins : list of tuples
        The reference bins for calculating tau.

    """
    left, _, right = interfaces
    if left == float('-inf') and (tau_region is None or tau_region == []):
        return [None]
    if left == float('-inf'):
        # Take the left of the tau_region to be the middle
        left = right-2*(right-tau_region[0])

    diff = (right-left)/nbins
    return [(left+diff*i, left+diff*(i+1)) for i in range(nbins)]


def skip_paths(weights, nacc, nskip):
    """Skips not counted paths.

    Parameters
    ----------
    weights : array of integers
        The weights of the paths.
    nacc : integer
        Number of accepted paths.
    nskip : integer
        Number of paths to skip.

    Returns
    -------
    new_weights : array of integers
        The new weights of the paths once skipped the nskip ones.
    new_acc : integer
        Recomputed number of accepted paths.

    """
    new_weights = []
    new_acc = nacc
    i = 0
    while nskip > 0:
        try:
            weight = weights[i]
        except IndexError as err:
            raise RuntimeError("Skipping more trajs than available") from err
        new_weight = max(weight-nskip, 0)
        diff = weight-new_weight
        new_acc -= 1
        nskip -= diff
        new_weights.append(new_weight)
        i += 1
    # Slicing on an index that does not exist gives an empty list
    # and always works
    new_weights += weights[i:]
    return new_weights, new_acc


def cross_dist_distr(pp_dic):
    """Return the distribution of lambda values for the LMR and RML paths.

    It calculates the distribution of lambda_max values for LM* paths, and
    the distribution of lambda_min values for RM* paths.

    Parameters
    ----------
    pp_dic : dict
        The partial-path dictionary. It contains the following keys:
        'lmrs' : array of strings (type of paths: "LML", "LMR", ...)
        'weights' : array of integers (weights of the paths)
        'interfaces' : array of floats (interfaces of the pathensemble)
        'ordermin' : array of floats (lambda_mins of the paths)
        'ordermax' : array of floats (lambda_maxs of the paths)

    Returns
    -------
    left : float
        The left interface of the pathensemble.
    middle : float
        The middle interface of the pathensemble.
    right : float
        The right interface of the pathensemble.
    percents : array of floats
        The distribution of lambda_max values for LM* paths, given at lambda
        values 'lambs' between middle and right interfaces.
    lambs : array of floats
        The lambda values at which the lambda_max distribution is given.
    percents2 : array of floats
        The distribution of lambda_min values for RM* paths, given at lambda
        values 'lambs2' between left and middle interfaces.
    lambs2 : array of floats
        The lambda values at which the lambda_min distribution is given.

    """
    lmrs = pp_dic['lmrs']
    weights = pp_dic['weights']
    interfaces = pp_dic['interfaces']
    ordermin = pp_dic['ordermin']
    ordermax = pp_dic['ordermax']

    # LM*:
    paths = select_with_or_masks(ordermax, [lmrs == "LML", lmrs == "LMR"])
    repeats = select_with_or_masks(weights, [lmrs == "LML", lmrs == "LMR"])
    # Ensure repeats are integer counts and guard empty cases
    if repeats.size != 0:
        repeats = np.array(repeats, dtype=int)
    else:
        repeats = np.array([], dtype=int)
    if paths.size == 0 or repeats.size == 0:
        repeat = np.array([])
    else:
        repeat = np.repeat(paths, repeats)
    left, middle, right = interfaces[0], interfaces[1], interfaces[2]
    percents = []

    if paths.size != 0:
        lambs = np.linspace(middle, np.max(paths), 100)
    else:
        lambs = np.linspace(middle, right, 100)
    for i in lambs:
        percents.append(np.sum(repeat >= i))

    # Convert to numpy array and guard against division by zero
    percents = np.array(percents, dtype=float)
    if percents.size == 0:
        percents = np.array([])
    elif np.allclose(percents[0], 0):
        percents = np.zeros_like(percents)
    else:
        percents = percents / percents[0]

    # RM*:
    paths2 = select_with_or_masks(
        ordermin,
        [lmrs == "RMR", lmrs == "RML"],
    )
    repeats2 = select_with_or_masks(
        weights,
        [lmrs == "RMR", lmrs == "RML"],
    )
    if repeats2.size != 0:
        repeats2 = np.array(repeats2, dtype=int)
    else:
        repeats2 = np.array([], dtype=int)
    if paths2.size == 0 or repeats2.size == 0:
        repeat2 = np.array([])
    else:
        repeat2 = np.repeat(paths2, repeats2)
    percents2 = []
    if paths2.size != 0:
        lambs2 = np.linspace(np.min(paths2), middle, 100)
    else:
        lambs2 = np.linspace(left, middle, 100)
    for i in lambs2:
        percents2.append(np.sum(repeat2 <= i))
    percents2 = np.array(percents2, dtype=float)
    if percents2.size == 0:
        percents2 = np.array([])
    elif np.allclose(percents2[-1], 0):
        percents2 = np.zeros_like(percents2)
    else:
        percents2 = percents2 / percents2[-1]
    return left, middle, right, percents, lambs, percents2, lambs2


def select_with_or_masks(array, masks):
    """Select elements of array that are True in at least one of the masks.

    Parameters
    ----------
    array : np.array
        Array to select elements from
    masks : list of masks
        where each mask is a boolean array with the same shape as array.

    Returns
    -------
    np.array
        A new array with the elements of A that are True in at least one of the
        masks.

    """
    # first check whether masks have the same shape as A
    for mask in masks:
        if mask.shape != array.shape:
            raise ValueError('Mask shape must match the input array shape')
    # now we can use the masks to select the elements of A
    union_mask = np.any(masks, axis=0).astype(bool)
    return array[union_mask]

#! /usr/bin/env python3
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Methods that will output results from the analysis functions.

The Methods defined here will also run the analysis and output
according to given settings.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

analyse_file (:py:func:`.analyse_file`)
    Method to analyse a file. For example, it can be used as

    >>> from pyretis.inout.analysisio import analyse_file
    >>> result, raw_data = analyse_file('cross', 'cross.txt', settings)

    To output the results, for instance by plotting, one can do:

    >>> figures = output_results(file_type, plotter, result, raw_data)
    >>> txt_file = output_results(file_type, txt_plotter, result, raw_data)


run_analysis_files (:py:func:`.run_analysis_files`)
    Methods to the analysis on a set of files. It will create output
    that can be used for reporting.
"""
import logging
import os
import numpy as np
from pyretis.core.units import create_conversion_factors
from pyretis.core.pathensemble import generate_ensemble_name
from pyretis.analysis import (analyse_flux, analyse_energies, analyse_orderp,
                              analyse_path_ensemble, analyse_repptis_ensemble,
                              match_probabilities, retis_flux, retis_rate,
                              perm_calculations)
from pyretis.inout.formats.formatter import format_number
from pyretis.inout.plotting import create_plotter
from pyretis.inout.plotting import TxtPlotter
from pyretis.inout.report import generate_report
from pyretis.inout.settings import SECTIONS, copy_settings
from pyretis.inout.formats.pathensemble import PathEnsembleFile
from pyretis.inout.formats.energy import EnergyFile
from pyretis.inout.formats.order import OrderFile
from pyretis.inout.formats.cross import CrossFile
from pyretis.inout.simulationio import OUTPUT_TASKS
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())

np.set_printoptions(legacy='1.21')

__all__ = ['analyse_file', 'run_analysis_files',
           'run_analysis', 'read_first_block',
           'repptis_running_pcross_analysis']


_ANALYSIS_FUNCTIONS = {'cross': analyse_flux,
                       'order': analyse_orderp,
                       'energy': analyse_energies,
                       'pathensemble': analyse_path_ensemble,
                       'pathensemble_repptis': analyse_repptis_ensemble}


# Input files for analysis
_FILES = {'md-flux': {},
          'md-nve': {},
          'md': {},
          'tis': {},
          'make-tis-files': {},
          'retis': {},
          'repptis': {}}
# Add files
for key in ('cross', 'energy', 'order'):
    _FILES['md-flux'][key] = OUTPUT_TASKS[key]['filename']
for key in ('energy',):
    _FILES['md-nve'][key] = OUTPUT_TASKS[key]['filename']
    _FILES['md'][key] = OUTPUT_TASKS[key]['filename']
for key in ('pathensemble',):
    for key2 in ('make-tis-files', 'tis', 'retis', 'repptis'):
        _FILES[key2][key] = OUTPUT_TASKS[key]['filename']


def run_analysis(settings):
    """Run a predefined analysis task.

    Parameters
    ----------
    settings : dict
        Simulation settings and settings for the analysis.

    Returns
    -------
    out : dict
        A dictionary with the results from the analysis. This dict
        can be used to generate a report.

    """
    runners = {'retis': run_retis_analysis,
               'tis': run_tis_analysis,
               'make-tis-files': run_tis_analysis,
               'md-flux': run_mdflux_analysis,
               'repptis': run_repptis_analysis}
    sim_task = settings['simulation']['task'].lower()
    report_dir = settings['analysis'].get('report-dir', None)
    plotter = create_plotter(settings['analysis']['plot'],
                             out_dir=report_dir)
    backup = settings['output'].get('backup', SECTIONS['output']['backup'])
    txt_plotter = TxtPlotter(settings['analysis']['txt-output'],
                             backup=backup,
                             out_dir=report_dir)
    if sim_task in runners:
        runner = runners[sim_task]
        return runner(settings, plotter, txt_plotter)
    msgtxt = f'Unknown analysis task "{sim_task}" requested!'
    logger.error(msgtxt)
    raise ValueError(msgtxt)


def get_path_ensemble_files(ensemble_number, settings, detect,
                            interfaces):
    """Return files for a single path ensemble.

    Here, we will return the files needed to analyse a single path
    ensemble and we will also return settings which can be used for
    the analysis.

    Parameters
    ----------
    ensemble_number : int
        This is the integer representing the ensemble.
    settings : dict
        The settings to use for an analysis/simulation
    detect : float or None
        The interface use for detecting if a path is successful for not.
    interfaces : list of floats
        The interfaces used for this particular path simulation.

    Returns
    -------
    local_settings : dict
        This dict contains settings which can be used for a initial
        flux analysis.
    files : list of tuples
        The tuples in this list are the files which can be analysed
        further, using the settings in `out[0]`.

    """
    sim_task = settings['simulation']['task'].lower()
    lsetting = copy_settings(settings)
    lsetting['simulation']['interfaces'] = interfaces
    lsetting['tis']['detect'] = detect
    lsetting['tis']['ensemble_number'] = ensemble_number
    run_path = settings['simulation']['exe-path']
    files = []
    for file_type in _FILES[sim_task]:
        filename = os.path.join(run_path,
                                generate_ensemble_name(ensemble_number),
                                _FILES[sim_task][file_type])
        if os.path.isfile(filename):
            files.append((file_type, filename))
            logger.debug('Adding file "%s" for analysis', filename)
        else:
            logger.debug('Could not find file %s', filename)
    return lsetting, files


def get_path_simulation_files(sim_settings):
    """Set up for analysis of TIS and RETIS simulations.

    For these kinds of simulations, we expect to analyse several
    directories with raw-data. For TIS we expect to find a directory
    with raw-data for the initial flux (named ``flux``) and the
    directories for the path simulations (named ``001``, ``002`` etc.
    for ensembles [0^+], [1^+] and so on). For RETIS, we expect to
    find the same directories but with a ``000`` (for the ``[0^-]``
    ensemble) rather than the ``flux`` directory.

    Parameters
    ----------
    sim_settings : dict
        The settings used for the (RE)TIS simulation(s). These settings
        are used to get the interfaces used and the path ensembles
        defined by these interfaces.

    Returns
    -------
    all_settings : list of dict
        This dict in this list contains settings which can be used for
        the analysis.
    all_files : list of lists of tuples
        `all_files[i]` is a list of files from path ensemble
        simulation `i`. For TIS, `all_files[0]` should be the files
        obtained in the initial flux simulation. These files can be
        analysed using the settings in `all_settings[i]`.

    """
    # Check if we can do flux analysis:
    all_files, all_settings = [], []
    interfaces = sim_settings['simulation']['interfaces']
    reactant = interfaces[0]
    product = interfaces[-1]
    task = sim_settings['simulation']['task']
    # for PPTIS or REPPTIS: memory
    mem = sim_settings.get(task, {}).get('memory', 1)

    if task == 'make-tis-files':
        all_files.append(None)
        all_settings.append(None)
    else:  # just add the 0 ensemble
        detect = None
        interface = [-float('inf'), reactant, reactant]
        setts, files = get_path_ensemble_files(0, sim_settings, detect,
                                               interface)
        all_files.append(files)
        all_settings.append(setts)
    for i, middle in enumerate(interfaces[:-1]):
        detect = interfaces[i + 1] if i < len(interfaces) else product
        if task in ['pptis', 'repptis']:
            reactant = interfaces[max(0, i - mem)]
            product = interfaces[min(i + mem, len(interfaces) - 1)]
            interface = [reactant, middle, product]
        else:
            interface = [reactant, middle, product]
        setts, files = get_path_ensemble_files(i + 1, sim_settings, detect,
                                               interface)
        all_files.append(files)
        all_settings.append(setts)
    return all_settings, all_files


def print_value_error(heading, value, rel_error, level=None):
    """Print out the matched probabilities."""
    log_fn = logger.progress if level == 'success' else logger.info
    if value is None or not np.isfinite(value):
        log_fn(f'{heading}: insufficient data to compute')
        if rel_error is None or not np.isfinite(rel_error):
            log_fn('(Relative error: nan %)')
        else:
            fmt_scale = format_number(rel_error * 100, 0.1, 100)
            msgtxt = f'(Relative error: {fmt_scale.rstrip()} %)'
            log_fn(msgtxt)
        return
    val = format_number(value, 0.1, 100)
    msgtxt = f'{heading}: {val}'
    log_fn(msgtxt.strip())
    fmt_scale = format_number(rel_error * 100, 0.1, 100)
    msgtxt = f'(Relative error: {fmt_scale.rstrip()} %)'
    log_fn(msgtxt)


def run_single_tis_analysis(settings, plotter, txt_plotter):
    """Run the analysis for a single TIS ensemble.

    Parameters
    ----------
    settings : dict
        The settings to use for an analysis/simulation.
    plotter : object like :py:class:`.Plotter`
        This is the object that handles the plotting.
    txt_plotter : object like :py:class:`.Plotter`
        This is the object that handles the text output.

    Returns
    -------
    out : list or dict or similar
        The output from the analysis.

    """
    sim = settings['simulation']
    tis = settings['ensemble'][0]['tis']
    sett, files = get_path_ensemble_files(
        tis['ensemble_number'], settings,
        tis.get('detect', sim['interfaces'][-1]),
        sim['interfaces'])
    msgtxt = f"Analysing ensemble {tis['ensemble_number']}"
    logger.info(msgtxt)
    result = run_analysis_files(sett, files, plotter, txt_plotter)
    report_txt = generate_report('tis', result, output='txt')[0]
    logger.progress(''.join(report_txt))
    return result


def run_tis_analysis(settings, plotter, txt_plotter):
    """Run the analysis for TIS.

    Parameters
    ----------
    settings : dict
        The settings to use for an analysis/simulation.
    plotter : object like :py:class:`.Plotter`
        This is the object that handles the plotting.
    txt_plotter : object like :py:class:`.Plotter`
        This is the object that handles the text output.

    Returns
    -------
    out : list or dict or similar
        The output from the analysis.

    """
    list_ens = settings['ensemble']
    if len(list_ens) == 1 and settings['simulation']['task'] == 'tis':
        return run_single_tis_analysis(settings, plotter, txt_plotter)
    results = {'cross': None, 'pathensemble': [], 'matched': None}
    all_settings, all_files = get_path_simulation_files(settings)
    nens = len(all_settings) - 1
    for i, (sett, files) in enumerate(zip(all_settings, all_files)):
        if i == 0:
            msgtxt = ('Initial flux is not calculated here.\n'
                      'Remember to calculate this separately!')
            logger.info(msgtxt)
        else:
            msgtxt = f'Analysing ensemble {i} of {nens}'
            logger.info(msgtxt)
            result = run_analysis_files(sett, files, plotter, txt_plotter)
            results['pathensemble'].append(result['pathensemble'])
            report_txt = generate_report('tis', result,
                                         output='txt')[0]
            logger.progress(''.join(report_txt))
    # match probabilities:
    out, fig, txt = analyse_and_output_matched(results['pathensemble'],
                                               plotter, txt_plotter,
                                               settings)
    results['matched'] = {'out': out, 'figures': fig, 'txtfile': txt}
    logger.info('Overall results')
    logger.info('===============')
    print_value_error('TIS Crossing probability',
                      out['prob'], out['relerror'], level='success')
    return results


def run_retis_analysis(settings, plotter, txt_plotter):
    """Run the analysis for RETIS.

    Parameters
    ----------
    settings : dict
        The settings to use for an analysis/simulation.
    plotter : object like :py:class:`.Plotter`
        This is the object that handles the plotting.
    txt_plotter : object like :py:class:`.Plotter`
        This is the object that handles the text output.

    """
    units = settings['system']['units']
    permeability = settings['simulation'].get('permeability', False)
    create_conversion_factors(units, **settings.get('unit-system', {}))
    all_settings, all_files = get_path_simulation_files(settings)
    results = {'cross': None,
               'pathensemble': [],
               'matched': None,
               'permeability': permeability}
    nens = len(all_settings) - 1
    for i, (sett, files) in enumerate(zip(all_settings, all_files)):
        msgtxt = f'Analysing ensemble {i} of {nens}'
        logger.info(msgtxt)
        if i == 0:
            result = run_analysis_files(sett, files, plotter, txt_plotter)
            results['pathensemble0'] = result['pathensemble']
            report_txt = generate_report('retis0', result, output='txt')[0]
            logger.progress(''.join(report_txt))
        else:
            result = run_analysis_files(sett, files, plotter, txt_plotter)
            results['pathensemble'].append(result['pathensemble'])
            report_txt = generate_report('tis', result,
                                         output='txt')[0]
            logger.progress(''.join(report_txt))
    # flux first:
    time_subcycles = settings['engine'].get('subcycles', 1)
    timestep = settings['engine']['timestep'] * time_subcycles
    flux, flux_error = retis_flux(results['pathensemble0']['out'],
                                  results['pathensemble'][0]['out'],
                                  timestep)
    results['flux'] = {'value': flux, 'error': flux_error,
                       'unit': units}
    # match probabilities:
    out, fig, txt = analyse_and_output_matched(results['pathensemble'],
                                               plotter, txt_plotter,
                                               settings,
                                               flux=flux)
    results['matched'] = {'out': out, 'figures': fig, 'txtfile': txt}
    rate, rate_error = retis_rate(out['prob'], out['relerror'],
                                  flux, flux_error)
    results['rate'] = {'value': rate, 'error': rate_error,
                       'unit': units}

    xi, tau, perm, xi_err, tau_err, perm_err = perm_calculations(
        results['pathensemble0'], out['prob'], out['relerror'])
    xi_figs, tau_figs = perm_figures(results['pathensemble0'], plotter)
    results['xi'] = {'value': xi,
                     'error': xi_err,
                     'fig': xi_figs}
    results['tau'] = {'value': tau,
                      'error': tau_err,
                      'fig': tau_figs}
    results['perm'] = {'value': perm,
                       'error': perm_err}

    logger.info('Overall results')
    logger.info('===============')
    print_value_error('RETIS Crossing probability',
                      out['prob'], out['relerror'], level='success')
    print_value_error(f'Initial flux (units 1/{units})', flux,
                      flux_error, level='success')
    print_value_error(f'Rate constant (units 1/{units})', rate,
                      rate_error, level='success')
    if permeability:
        print_value_error('Xi', xi, xi_err, level='success')
        print_value_error(f'Tau/dz (unit: {units}/OP-unit)',
                          tau, tau_err, level='success')
        print_value_error(f'Permeability (unit: OP-unit/{units})',
                          perm, perm_err, level='success')
    return results


def run_repptis_analysis(settings, plotter, txt_plotter):
    """Run the analysis for PPRETIS.

    Parameters
    ----------
    settings : dict
        The settings to use for an analysis/simulation.
    plotter : object like :py:class:`.Plotter`
        This is the object that handles the plotting.
    txt_plotter : object like :py:class:`.Plotter`
        This is the object that handles the text output.

    """
    units = settings['system']['units']
    permeability = settings['simulation'].get('permeability', False)
    create_conversion_factors(units, **settings.get('unit-system', {}))
    all_settings, all_files = get_path_simulation_files(settings)
    results = {'cross': None,
               'pathensemble': [],
               'matched': None,
               'permeability': permeability,
               'reptis': {}}
    nens = len(all_settings) - 1
    probs = []
    ptypes = []
    for i, (sett, files) in enumerate(zip(all_settings, all_files)):
        msgtxt = f'Analysing ensemble {i} of {nens}'
        # We change filetype to repptis instead of retis for anaylsis
        repptis_str = 'pathensemble_repptis'
        files0 = [(repptis_str, files[0][1])]
        logger.info(msgtxt)
        if i == 0:
            result = run_analysis_files(sett, files0, plotter, txt_plotter)
            results['pathensemble0'] = result['pathensemble_repptis']
            report_txt = generate_report('retis0', result, output='txt')[0]
            logger.progress(''.join(report_txt))
        else:
            result = run_analysis_files(sett, files0, plotter, txt_plotter)
            results['pathensemble'].append(result['pathensemble_repptis'])
            pp_lmr = result[repptis_str]['out']['prun_sl'][-1]
            pp_rml = result[repptis_str]['out']['prun_sr'][-1]
            ptype = result[repptis_str]['out']['ptypes']
            ptypes.append(ptype)
            pp_lml = 1 - pp_lmr
            pp_rmr = 1 - pp_rml
            pp_lm_err = result[repptis_str]['out']['blockerror_sl'][2]
            pp_rm_err = result[repptis_str]['out']['blockerror_sr'][2]
            probs.append({"LMR": {"mean": pp_lmr,
                                  "std": pp_lm_err},
                          "RML": {"mean": pp_rml,
                                  "std": pp_rm_err},
                          "RMR": {"mean": pp_rmr,
                                  "std": pp_rm_err},
                          "LML": {"mean": pp_lml,
                                  "std": pp_lm_err}, })

            report_txt = generate_report('pptis', result, output='txt')[0]
            logger.progress(''.join(report_txt))

    # flux first:
    time_subcycles = settings['engine'].get('subcycles', 1)
    timestep = settings['engine']['timestep'] * time_subcycles
    flux, flux_error = retis_flux(results['pathensemble0']['out'],
                                  results['pathensemble'][0]['out'],
                                  timestep)
    results['flux'] = {'value': flux, 'error': flux_error,
                       'unit': units}

    results.update(repptis_running_pcross_analysis(ptypes))
    figures = plotter.output_matched_probability(
            [],
            settings['simulation']['interfaces'],
            results, reptis=True
    )
    txt_plotter.output_matched_probability(
            [],
            settings['simulation']['interfaces'],
            results, reptis=True
    )
    results['matched_fig'] = figures
    results['cross'] = {'value': results['pcrossrun'][-1],
                        'error': results['pcrossrun_hav_rele']}

    rate, rate_error = retis_rate(results['pcrossrun'][-1],
                                  results['pcrossrun_hav_rele'],
                                  flux, flux_error)
    results['rate'] = {'value': rate, 'error': rate_error,
                       'unit': units}

    xi, tau, perm, xi_err, tau_err, perm_err = perm_calculations(
        results['pathensemble0'], results['pcrossrun'][-1],
        results['pcrossrun_hav_rele'])
    xi_figs, tau_figs = perm_figures(results['pathensemble0'], plotter)
    results['xi'] = {'value': xi,
                     'error': xi_err,
                     'fig': xi_figs}
    results['tau'] = {'value': tau,
                      'error': tau_err,
                      'fig': tau_figs}
    results['perm'] = {'value': perm,
                       'error': perm_err}

    logger.info('Overall results')
    logger.info('===============')
    print_value_error('REPPTIS Crossing probability',
                      results['pcrossrun'][-1],
                      results['pcrossrun_hav_rele'], level='success')
    print_value_error(f'Initial flux (units 1/{units})', flux,
                      flux_error, level='success')
    print_value_error(f'Rate constant (units 1/{units})', rate,
                      rate_error, level='success')
    if permeability:
        print_value_error('Xi', xi, xi_err, level='success')
        print_value_error(f'Tau/dz (unit: {units}/OP-unit)',
                          tau, tau_err, level='success')
        print_value_error(f'Permeability (unit: OP-unit/{units})',
                          perm, perm_err, level='success')
    return results


def run_mdflux_analysis(settings, plotter, txt_plotter):
    """Run the analysis for the md flux simulation.

    Parameters
    ----------
    settings : dict
        The settings to use for an analysis/simulation.
    plotter : object like :py:class:`.Plotter`
        This is the object that handles the plotting.
    txt_plotter : object like :py:class:`.Plotter`
        This is the object that handles the text output.

    Returns
    -------
    out : list or dict or similar
        The output from the analysis.

    """
    sim = settings['simulation']
    sim_task = sim['task']
    exe_path = settings['simulation'].get('exe-path', os.getcwd())
    files = []
    for file_type in _FILES[sim_task]:
        filename = _FILES[sim_task][file_type]
        if os.path.isfile(os.path.join(exe_path, filename)):
            files.append((file_type, os.path.join(exe_path, filename)))
    msgtxt = 'Running analysis of a MD flux simulation...'
    logger.info(msgtxt)
    result = run_analysis_files(settings, files, plotter, txt_plotter)
    report_txt = generate_report('md-flux', result, output='txt')[0]
    logger.progress(''.join(report_txt))
    return result


def run_analysis_files(settings, files, plotter, txt_plotter):
    """Run the analysis on a collection of files.

    Parameters
    ----------
    settings : dict
        This dict contains settings which dictate how the
        analysis should be performed and it should also contain
        information on how the simulation was performed.
    files : list of tuples
        This list contains the raw files to be analysed. The
        tuples are on format ('filetype', 'filename').
    plotter : object like :py:class:`.Plotter`
        This is the object that handles the plotting.
    txt_plotter : object like :py:class:`.Plotter`
        This is the object that handles the text output.

    Returns
    -------
    out: dict
        The results from the analysis.

    """
    results = {}
    for (file_type, file_name) in files:
        result, raw_data = analyse_file(file_type, file_name, settings)
        figures = output_results(file_type, plotter, result, raw_data)
        txt_file = output_results(file_type, txt_plotter, result, raw_data)
        results[file_type] = {'out': result,
                              'figures': figures,
                              'txtfile': txt_file}
    return results


def read_first_block(file_type, file_name):
    """Read the first block of data from a file.

    This method will read the first block of data from a file and
    immediately return.

    Parameters
    ----------
    file_type : string
        A string used for selecting the file object to use here.
    file_name : string
        The file to open.

    Returns
    -------
    out : numpy.array
        The raw data read from the file.

    """
    first_block = None
    class_map = {'energy': EnergyFile, 'order': OrderFile, 'cross': CrossFile}
    if file_type not in class_map:
        raise ValueError('Unknown file type requested in analysis.')

    klass = class_map[file_type]

    with klass(file_name, 'r') as fileobj:
        for block in fileobj.load():
            if first_block is None:
                first_block = block
            else:
                logger.warning(
                    (
                        'Noticed a second block in the input file "%s".\n'
                        'This will be ignored by the flux analysis.\n'
                        'Are you are running the analysis with '
                        'the correct input?'
                    ),
                    file_name
                )
                break
    if first_block is None:
        return None
    return first_block['data']


def output_results(file_type, plotter, result, rawdata):
    """Plot the results.

    Plotting might here refer to both figures and text.

    Parameters
    ----------
    file_type : string
        This determines what we are going to plot.
    plotter : object like :py:class:`.Plotter`
        This is the object that handles the plotting.
    result : list of lists or dicts
        This contains the results from a specific analysis.
    rawdata : list of floats, lists or objects
        The raw data with analysis results.

    Returns
    -------
    out : list of strings
        These are the files that the plotter created.

    """
    if file_type == 'cross':
        return plotter.output_flux(result)
    if file_type == 'order':
        return plotter.output_orderp(result, rawdata)
    if file_type == 'energy':
        return plotter.output_energy(result, rawdata)
    if file_type == 'pathensemble':
        if 'pathlength' not in result:
            return None
        return plotter.output_path(result, rawdata)
    if file_type == 'pathensemble_repptis':
        if 'pathlength' not in result:
            return None
        return plotter.output_pppath(result, rawdata)
    return None


def analyse_file(file_type, file_name, settings):
    """Run analysis on the given file.

    This function is included for convenience so that we can call an
    analysis like ``analyse_file('cross', 'cross.txt')`` i.e. it should
    automatically open the file and apply the correct analysis
    according to a given file type. Here we return a function to do the
    analysis, so we are basically wrapping one of the analysis
    functions. This is done in case we wish to rerun the analysis but
    with different settings for instance.

    Parameters
    ----------
    file_type : string
        This is the type of file we are to analyse.
    file_name : string
        The file name to open.
    settings : dict
        This dict contains settings which dictate how the
        analysis should be performed and information on how the
        simulation was performed.

    Returns
    -------
    results : list or dict
        The output from the analysis
    raw_data : list, numpy.array or another type of object
        The raw data used in the analysis.

    """
    function = _ANALYSIS_FUNCTIONS.get(file_type, None)
    if function is None:
        msgtxt = f'Unknown file type "{file_type}" requested!'
        logger.error(msgtxt)
        raise ValueError(msgtxt)
    results, raw_data = None, None
    if file_type in ('energy', 'order', 'cross'):
        # Read the first block of raw data and analyse it:
        raw_data = read_first_block(file_type, file_name)
        results = function(raw_data, settings)
    elif file_type in ('pathensemble', 'pathensemble_repptis'):
        # Create a PathEnsembleFile object to use for the analysis.
        interfaces = settings['simulation']['interfaces']
        ensemble_settings = {
            'ensemble_number': settings['tis']['ensemble_number'],
            'interfaces': interfaces,
            'detect': settings['tis'].get('detect', interfaces[1])
        }
        raw_data = PathEnsembleFile(file_name, 'r', ensemble_settings)
        results = function(raw_data, settings)
    return results, raw_data


def analyse_and_output_matched(raw_data, plotter, txt_plotter,
                               settings, flux=None):
    """Analyse and output matched probability.

    This will calculate the over-all crossing probability by combining
    results from many path simulations.

    Parameters
    ----------
    raw_data : list, numpy.array or another type of object
        The raw data used in the analysis.
    plotter : object like :py:class:`.Plotter`
        This is the object that handles the plotting.
    txt_plotter : object like :py:class:`.Plotter`
        This is the object that handles the text output.
    settings : dict
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

    flux : float, optional
        The computed flux to be used only in RETIS if the running average
        of the computed rate is desired.

    Returns
    -------
    out[0] : dict
        This dict contains the results from the analysis
    out[1] : list of dicts
        A dictionary with the figure files created (if any).
    out[2] : list of strings
        A list with the text files created (if any).

    """
    path_results, ensemble_names, detect = [], [], []
    interface_left = None
    for ensemble in raw_data:
        path_results.append(ensemble['out'])
        ensemble_names.append(ensemble['out']['ensemble'])
        detect.append(ensemble['out']['detect'])
        if interface_left is None:
            interface_left = ensemble['out']['interfaces'][0]
    result = match_probabilities(path_results, detect, settings)
    if flux is not None:
        result['overall-rrun'] = flux*result['overall-prun']
    # for the figure, we add the A interface:
    detect_plot = [interface_left] + detect
    figures = plotter.output_matched_probability(ensemble_names,
                                                 detect_plot,
                                                 result)
    outtxt = txt_plotter.output_matched_probability(ensemble_names,
                                                    detect,
                                                    result)
    return result, figures, outtxt


def perm_figures(results0, plotter):
    """Generate permanence (xi/tau) figures for repptis analysis.

    Parameters
    ----------
    results0 : dict
        Analysis results for the [0^-] ensemble, expected to contain an
        'out' key with 'xi', 'tau', and related arrays.
    plotter : object like :py:class:`.Plotter`
        The plotter used to create output figures.

    Returns
    -------
    xi_figures : list of str or None
        Paths to the xi figures, or None if xi data is absent.
    tau_figures : list of str or None
        Paths to the tau figures, or None if tau data is absent.

    """
    out = results0['out']
    # Quick check to see if we need to do something
    if 'xi' not in out:
        return None, None
    xi = out['xi']
    xierror = out['xierror']

    bin_edges = out['tau_bin']
    tau = out['tau']
    tauerror = out['tauerror']
    tau_ref = out['tau_ref']
    tau_ref_bins = out['tau_ref_bins']
    xi_figures = plotter.output_xi(xi, xierror)
    tau_figures = plotter.output_tau(tau, tauerror, bin_edges,
                                     tau_ref, tau_ref_bins)
    return xi_figures, tau_figures


def repptis_running_pcross_analysis(l_ptypes):
    """Recursive block analysis for repptis.

    Parameters
    ----------
    l_ptypes : list
        List of cumsum ptypes for the ensembles [0+-'], [1+-], ... [N-1+-].

    Returns
    -------
    results : dict
        The results from the analysis.

    """
    l_locps = []
    for ptypes in l_ptypes:  # loop over ensembles
        # if a division by zero occurs, we set the probability to NaN
        # casting to floats for true_divide
        ptypes = ptypes.astype(float)
        pmps = np.divide(ptypes[:, 0], ptypes[:, 0]+ptypes[:, 1],
                         out=np.full_like(ptypes[:, 0], np.nan),
                         where=(ptypes[:, 0]+ptypes[:, 1]) != 0)
        pmms = np.divide(ptypes[:, 1], ptypes[:, 0]+ptypes[:, 1],
                         out=np.full_like(ptypes[:, 0], np.nan),
                         where=(ptypes[:, 0]+ptypes[:, 1]) != 0)
        ppps = np.divide(ptypes[:, 2], ptypes[:, 2]+ptypes[:, 3],
                         out=np.full_like(ptypes[:, 2], np.nan),
                         where=(ptypes[:, 2]+ptypes[:, 3]) != 0)
        ppms = np.divide(ptypes[:, 3], ptypes[:, 2]+ptypes[:, 3],
                         out=np.full_like(ptypes[:, 2], np.nan),
                         where=(ptypes[:, 2]+ptypes[:, 3]) != 0)
        # we have to_list them, as we will append them until all have N_max cyc
        l_locps.append([pmms.tolist(), pmps.tolist(),
                        ppps.tolist(), ppms.tolist()])
    # So, l_locps: list of lists of lists: n_ensembles x 4 x n_blocks
    # extend the l_locps to the same length. We copy the last element of the
    # list until the length is the same as the longest list
    max_len = max(len(ll[0]) for ll in l_locps)
    for ll in l_locps:
        while len(ll[0]) < max_len:
            for i in range(4):
                ll[i].append(ll[i][-1])
    # We can now calculate the global crossing probabilities, for each cycle
    pcross_run = []
    for icyc in range(max_len):
        # We have to get the local crossing probabilities for each ensemble
        l_pcross = [[ll[i][icyc] for i in range(4)] for ll in l_locps]
        # so l_pcross is shape n_ensembles x 4, but we need 4 lists of
        # length n_ensembles
        l_pcross = list(map(list, zip(*l_pcross)))
        # Now we calculate the global crossing probabilities
        _, _, pcross = get_global_probz(*l_pcross)
        # We append the last element of pcross to the running list
        pcross_run.append(pcross[-1] if pcross is not np.nan else np.nan)
    # Now do the recursive block analysis
    rel_errors, nstatineff, half_av_err = recursive_block_analysis(pcross_run)

    results = {'pcross': pcross,
               'pcrossrun': pcross_run,
               'pcrossrun_rele': rel_errors,
               'nstatineff': nstatineff,
               'pcrossrun_hav_rele': half_av_err}

    # make into retis format
    results['overall-cycle'] = list(range(len(pcross_run)))
    results['overall-prun'] = pcross_run
    results['overall-error'] = [list(range(len(rel_errors))), 0, 0,
                                rel_errors, rel_errors[-1]]
    return results


def get_global_probz(pmps, pmms, ppps, ppms):
    """Return the global crossing probabilities for the PPTIS simulation.

    This follows the recursive relation
    P_plus[j] = (pLMR[j-1]*P_plus[j-1]) / (pLMR[j-1]+pLML[j-1]*P_min[j-1])
    P_min[j] = (pRML[j-1]P_min[j-1])/(pLMR[j-1]+pLML[j-1]*P_min[j-1])
    where P_plus[0] = 1, P_min[0] = 1, and j is the ensemble number.

    Parameters
    ----------
    pmps : list of floats
        The local probability (P) of having type LMR (MinusPlus) for each
        ensemble (S).
    pmms : list of floats
        Local LML
    ppps : list of floats
        Local RMR
    ppms : list of floats
        Local RML

    Returns
    -------
    pmin : list of floats
        Float per ensemble i. Represents the probability of crossing A earlier
        than i+1, given that the path crossed i.
    pplus : list of floats
        Float per ensemble i. Represents the probability of crossing i+1
        earlier than A, given that the path crossed i.
    pcross : list of floats
        Float per ensemble i. Represents the TIS probability of crossing i+1.
    """
    # if there is any NaN in pmps, pmms, ppps, ppms, return NaN
    nan_arr = [np.nan] * (len(pmps) + 1)
    if np.isnan(pmps).any() or np.isnan(pmms).any() or \
            np.isnan(ppps).any() or np.isnan(ppms).any():
        return nan_arr, nan_arr, nan_arr
    pplus, pmin, pcross = [1.], [1.], [1., pmps[0]]
    for i, pmp, pmm, _, ppm in zip(range(len(pmps)), pmps, pmms, ppps, ppms):
        if i == 0:  # This is [0^{\pm}'], so skip
            continue
        # Calculate the global crossing probabilities
        # if divide by zero, or divide by NaN, return NaN
        if pmp + pmm * pmin[-1] == 0 or\
                (pmp is np.nan) or\
                (pmm * pmin[-1] is np.nan):
            return nan_arr, nan_arr, nan_arr
        pplus.append((pmp*pplus[-1])/(pmp+pmm*pmin[-1]))
        pmin.append((ppm*pmin[-1])/(pmp+pmm*pmin[-1]))
        # Calculate the TIS probabilities
        pcross.append(pmps[0]*pplus[-1])
    return pmin, pplus, pcross


def recursive_block_analysis(flist, minblocks=5):
    """Run a recursive block analysis on (time-series) data f.

    Parameters
    ----------
    flist : list of floats
        The data to analyse.
    minblocks : int, optional
        The minimum number of blocks to use for the analysis.

    Returns
    -------
    rel_errors : list of floats
        The relative errors for each block length.
    nstatineff : float
        The statistical inefficiency.
    half_av_err : float
        The average error for the second half of the relative errors.

    """
    # flist is a list of floats, but we need to remove the np.nan's and 0's
    flist = np.array(flist)
    flist = flist[~np.isnan(flist)]
    flist = flist[flist != 0]
    # We now run the recursive blocking analysis
    # If we don't have enough data, we don't want the entire analysis to
    # fail. So just return NaNs, such that we still get those beautiful local
    # plots in the report PDF.
    if len(flist) < minblocks:
        return [np.nan], np.nan, np.nan
    bestavg = flist[-1]
    maxblocklen = int(len(flist)/minblocks)
    rel_errors = []
    for blocklen in range(1, maxblocklen+1):
        runav_red = flist[blocklen-1::blocklen]
        blocks = recursive_blocks(runav_red)
        sum_quad_diff = np.sum(np.fromiter(((b-bestavg)**2 for b in blocks),
                               dtype=float))
        nb = len(blocks)
        aerr2 = sum_quad_diff/(nb*(nb-1))
        aerr = np.sqrt(aerr2)
        rerr = aerr/bestavg
        rel_errors.append(rerr)

    # We will now take the mean of the second half of the relative errors.
    # The early return above guarantees len(rel_errors) >= 1 so the slice
    # below is non-empty.
    half_av_err = float(np.mean(rel_errors[int(len(rel_errors) / 2):]))

    # Guard division by zero / non-finite first sample when computing the
    # statistical inefficiency (a zero or NaN first relative error would
    # otherwise trigger a numpy RuntimeWarning).
    if (rel_errors[0] == 0
            or not np.isfinite(rel_errors[0])
            or np.isnan(half_av_err)):
        nstatineff = np.nan
    else:
        nstatineff = (half_av_err / rel_errors[0]) ** 2

    # Now we have to add the NaN's back in, as 0's
    return rel_errors, nstatineff, half_av_err


def recursive_blocks(recu):
    """Convert cumulative block averages to per-block values.

    Parameters
    ----------
    recu : list of floats
        Cumulative block averages where index ``i`` represents the average
        over ``i+1`` blocks.

    Returns
    -------
    blocks : list of floats
        Per-block values reconstructed from the cumulative averages.

    """
    blocks = []
    for i, val in enumerate(recu):
        if i == 0:
            blocks.append(val)
        else:
            blocks.append((i + 1) * val - i * recu[i - 1])
    return blocks

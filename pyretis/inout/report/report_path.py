# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Method for generating reports from path sampling simulations.

The reports are useful for displaying results from the analysis.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

generate_report_tis (:py:func:`.generate_report_tis`)
    Generate a report for the overall results from a TIS simulation.

generate_report_tis_path (:py:func:`.generate_report_tis_path`)
    Generate a report for a single TIS simulation.

generate_report_retis (:py:func:`.generate_report_retis`)
    Generate a report for the overall results from a RETIS simulation.
"""
import logging
from pyretis.inout.report.markup import (generate_rst_table,
                                         generate_latex_table)
from pyretis.inout.formats.formatter import apply_format, format_number
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


__all__ = ['generate_report_tis', 'generate_report_tis_path',
           'generate_report_retis']


def generate_report_tis_path(analysis, output='rst'):
    """Generate a report for a single TIS simulation.

    Parameters
    ----------
    analysis : dict
        This is the output from the analysis.
    output : string, optional
        This is the desired output format. It must match one of the
        formats defined in
        :py:const:`pyretis.inout.report.report._TEMPLATES`. Default
        is 'rst' (reStructuredText).

    Returns
    -------
    out[0] : string
        The generated report in the desired format.
    out[1] : string
        The file extension (i.e. file type) for the generated report.

    """
    result = analysis['pathensemble']
    if 'prun' not in result['out']:
        return {'ensemble': result['out'].get('ensemble', ''),
                'figures': {},
                'tables': {'interfaces': None,
                           'probability': None,
                           'path': None,
                           'efficiency': None}}
    report = {'ensemble': result['out']['ensemble'],
              'figures': _get_path_figures(result['figures']),
              'tables': {'interfaces': None,
                         'probability': None,
                         'path': None,
                         'efficiency': None}}
    # Create tables
    results = [result]
    tables = report['tables']
    tables['interfaces'] = _table_interface(results, fmt=output)[1]
    tables['probability'] = _table_probability(results, fmt=output)[1]
    tables['efficiency'] = _table_efficiencies(results, fmt=output)[1]
    return report


def generate_report_pptis_path(analysis, output='rst'):
    """Generate a report for a single TIS simulation.

    Parameters
    ----------
    analysis : dict
        This is the output from the analysis.
    output : string, optional
        This is the desired output format. It must match one of the
        formats defined in
        :py:const:`pyretis.inout.report.report._TEMPLATES`. Default
        is 'rst' (reStructuredText).

    Returns
    -------
    out[0] : string
        The generated report in the desired format.
    out[1] : string
        The file extension (i.e. file type) for the generated report.

    """
    result = analysis['pathensemble_repptis']
    report = {'ensemble': result['out']['ensemble'],
              'figures': _get_path_figures(result['figures']),
              'tables': {'interfaces': None,
                         'probability': None,
                         'path': None,
                         'efficiency': None}}
    # Create tables
    results = [result]
    tables = report['tables']
    tables['interfaces'] = _table_interface(results, fmt=output)[1]
    tables['probability'] = _table_probability_repptis(results, fmt=output)[1]
    tables['efficiency'] = _table_efficiencies(results, fmt=output)[1]
    return report


def generate_report_tis(analysis, output='rst'):
    """Generate a report for the over-all results from a TIS simulation.

    Parameters
    ----------
    analysis : dict
        This is the output from the analysis.
    output : string, optional
        This is the desired output format. It must match one of the
        formats defined in
        :py:const:`pyretis.inout.report.report._TEMPLATES`. Default
        is 'rst' (reStructuredText).

    Returns
    -------
    out[0] : string
        The generated report in the desired format.
    out[1] : string
        The file extension (i.e. file type) for the generated report.

    """
    report = {'figures': {'tis': [],
                          'matched': None},
              'tables': {'interfaces': None,
                         'probability': None,
                         'path': None,
                         'efficiency': None},
              'numbers': {'pcross': None, 'perr': None}}
    results = analysis['pathensemble']
    figures = report['figures']
    # Add figures:
    for result in results:
        figures['tis'].append(_get_path_figures(result['figures']))
    # Get matched result:
    matched_fig = analysis['matched']['figures']
    matched_out = analysis['matched']['out']
    figures['matched'] = matched_fig.get('matched-probability', None)
    figures['total'] = matched_fig.get('total-probability', None)
    figures['progress'] = matched_fig.get('overall-prun', None)
    figures['error'] = matched_fig.get('overall-err', None)
    # Get numbers:
    report['numbers']['pcross'] = format_number(matched_out['prob'], 0.1, 1)
    report['numbers']['perr'] = format_number(matched_out['relerror'] * 100,
                                              0.1, 100)
    # Get tables:
    tables = report['tables']
    tables['interfaces'] = _table_interface(results, fmt=output)[1]
    tables['probability'] = _table_probability(results, fmt=output)[1]
    tables['efficiency'] = _table_efficiencies(results, fmt=output)[1]
    return report


def generate_report_retis(analysis, output='rst'):
    """Generate a report for the over-all results from a RETIS simulation.

    Parameters
    ----------
    analysis : dict
        This is the output from the analysis.
    output : string, optional
        This is the desired output format. It must match one of the
        formats defined in
        :py:const:`pyretis.inout.report.report._TEMPLATES`. Default
        is 'rst' (reStructuredText).

    Returns
    -------
    out[0] : string
        The generated report in the desired format.
    out[1] : string
        The file extension (i.e. file type) for the generated report.

    """
    report = {'figures': {'tis': [],
                          'matched': None},
              'tables': {'interfaces': None,
                         'interfaces0': None,
                         'path0': None,
                         'probability': None,
                         'path': None,
                         'efficiency': None,
                         'summary': None},
              'numbers': {'pcross': None, 'perr': None, 'flux': None,
                          'flux-err': None},
              'text': {'flux-unit': None}}
    result0 = analysis['pathensemble0']
    results = [result0] + analysis['pathensemble']
    figures = report['figures']
    permeability = result0['out'].get('permeability', False)
    # Add figures:
    for result in results[1:]:
        figures['tis'].append(_get_path_figures(result['figures']))
    figures['path0'] = _get_path_figures(result0['figures'])
    # Get matched result:
    matched_fig = analysis['matched']['figures']
    matched_out = analysis['matched']['out']
    figures['matched'] = matched_fig.get('matched-probability', None)
    figures['total'] = matched_fig.get('total-probability', None)
    figures['progress'] = matched_fig.get('overall-rrun', None)
    figures['error'] = matched_fig.get('overall-err', None)
    if permeability:
        figures['xi'] = analysis['xi']['fig']['xirun']
        figures['xierror'] = analysis['xi']['fig']['xierror']
        figures['tau'] = analysis['tau']['fig']['taurun']
        figures['tauerror'] = analysis['tau']['fig']['tauerror']
        figures['tauref'] = analysis['tau']['fig']['tauref']

    # Get numbers:
    numbers = report['numbers']
    numbers['pcross'] = format_number(matched_out['prob'], 0.1, 1)
    numbers['perr'] = format_number(matched_out['relerror'] * 100, 0.1, 100)
    numbers['flux'] = format_number(analysis['flux']['value'], 0.1, 100)
    numbers['flux-err'] = format_number(analysis['flux']['error'] * 100,
                                        0.1, 100)
    numbers['rate'] = format_number(analysis['rate']['value'], 0.1, 100)
    numbers['rate-err'] = format_number(analysis['rate']['error'] * 100,
                                        0.1, 100)
    numbers['permeability'] = permeability
    if permeability:
        numbers['xi'] = format_number(analysis['xi']['value'], 0.1, 100)
        numbers['xi-err'] = format_number(analysis['xi']['error'] * 100,
                                          0.1, 100)
        numbers['tau'] = format_number(analysis['tau']['value'], 0.1, 100)
        numbers['tau-err'] = format_number(analysis['tau']['error'] * 100,
                                           0.1, 100)
        numbers['perm'] = format_number(analysis['perm']['value'], 0.1, 100)
        numbers['perm-err'] = format_number(analysis['perm']['error'] * 100,
                                            0.1, 100)

    report['text']['flux-unit'] = analysis['flux']['unit']
    # Get tables:
    tables = report['tables']
    tables['interfaces'] = _table_interface(results, fmt=output)[1]
    tables['probability'] = _table_probability(results[1:], fmt=output)[1]
    tables['efficiency'] = _table_efficiencies(results, fmt=output)[1]
    tables['summary'] = _table_summary(report, fmt=output)[1]
    return report


def generate_report_retis0(analysis, output='txt'):
    """Generate a report for [0^-] in a RETIS/REPTIS simulation.

    This method is primarily intended for displaying results on the
    screen during the analysis.

    Parameters
    ----------
    analysis : dict
        This is the output of the analysis.
    output : string, optional
        This is the desired output format. It must match one of the
        formats defined in
        :py:const:`pyretis.inout.report.report._TEMPLATES`. Default
        is 'rst' (reStructuredText).

    Returns
    -------
    out[0] : string
        The generated report in the desired format.
    out[1] : string
        The file extension (i.e. file type) for the generated report.

    """
    if output != 'txt':
        output = 'txt'
        logger.warning('Report for [0^-] is only indended as "txt"'
                       '\nOutput format "%s" ignored', output)
    if 'pathensemble' in analysis:
        result = analysis['pathensemble']
    else:
        result = analysis['pathensemble_repptis']
    flux = result['out']['fluxlength']
    report = {
        'ensemble': result['out']['ensemble'],
        'numbers': {'fluxlength': format_number(flux[0], 0, 10000),
                    'fluxlengtherror': format_number(flux[1] * 100, 0, 100)},
        'tables': {'interfaces0': _table_interface0([result], fmt=output)[1]},
    }
    return report


def generate_report_repptis(analysis, output='rst'):
    """Generate a report for the over-all results from a PPRETIS simulation.

    Parameters
    ----------
    analysis : dict
        This is the output from the analysis.
    output : string, optional
        This is the desired output format. It must match one of the
        formats defined in
        :py:const:`pyretis.inout.report.report._TEMPLATES`. Default
        is 'rst' (reStructuredText).

    Returns
    -------
    out[0] : string
        The generated report in the desired format.
    out[1] : string
        The file extension (i.e. file type) for the generated report.

    """
    report = {'figures': {'tis': [],
                          'matched': None},
              'tables': {'interfaces': None,
                         'probability': None,
                         'path': None,
                         'efficiency': None,
                         'summary': None},
              'numbers': {'pcross': None, 'perr': None, 'flux': None,
                          'flux-err': None},
              'text': {'flux-unit': None}}

    result0 = analysis['pathensemble0']
    results = [result0] + analysis['pathensemble']
    figures = report['figures']
    permeability = result0['out'].get('permeability', False)
    # Add figures:
    for result in results[1:]:
        figures['tis'].append(_get_path_figures(result['figures']))
    figures['path0'] = _get_path_figures(result0['figures'])
    # Get matched result:
    matched_fig = analysis['matched_fig']
    figures['matched'] = matched_fig.get('matched-probability', None)
    figures['progress'] = matched_fig.get('overall-prun', None)
    figures['error'] = matched_fig.get('overall-err', None)
    if permeability:
        figures['xi'] = analysis['xi']['fig']['xirun']
        figures['xierror'] = analysis['xi']['fig']['xierror']
        figures['tau'] = analysis['tau']['fig']['taurun']
        figures['tauerror'] = analysis['tau']['fig']['tauerror']
        figures['tauref'] = analysis['tau']['fig']['tauref']
    # Get numbers:
    numbers = report['numbers']
    numbers['pcross'] = format_number(analysis['cross']['value'], 0.1, 1)
    numbers['perr'] = format_number(analysis['cross']['error'] * 100,
                                    0.1, 100)
    numbers['flux'] = format_number(analysis['flux']['value'], 0.1, 100)
    numbers['flux-err'] = format_number(analysis['flux']['error'] * 100,
                                        0.1, 100)
    numbers['rate'] = format_number(analysis['rate']['value'], 0.1, 100)
    numbers['rate-err'] = format_number(analysis['rate']['error'] * 100,
                                        0.1, 100)
    numbers['permeability'] = permeability
    if permeability:
        numbers['xi'] = format_number(analysis['xi']['value'], 0.1, 100)
        numbers['xi-err'] = format_number(analysis['xi']['error'] * 100,
                                          0.1, 100)
        numbers['tau'] = format_number(analysis['tau']['value'], 0.1, 100)
        numbers['tau-err'] = format_number(analysis['tau']['error'] * 100,
                                           0.1, 100)
        numbers['perm'] = format_number(analysis['perm']['value'], 0.1, 100)
        numbers['perm-err'] = format_number(analysis['perm']['error'] * 100,
                                            0.1, 100)

    report['text']['flux-unit'] = 'time-unit'  # matched['flux']['unit']
    # Get tables:
    tables = report['tables']
    tables['interfaces'] = _table_interface_repptis(results, fmt=output)[1]
    tables['probability'] = _table_probability_repptis(results[1:],
                                                       fmt=output)[1]
    tables['efficiency'] = _table_efficiencies(results, fmt=output)[1]
    tables['summary'] = _table_summary(report, fmt=output)[1]
    return report


def _table_interface_repptis(results, fmt='rst'):
    """Generate the table for the interfaces.

    This table will display the location of the different interfaces.

    Parameters
    ----------
    results : list of dicts
        These are the results of the analysis.
    fmt : string, optional
        Determines if we create reStructuredText ('rst') or
        latex ('tex').

    Returns
    -------
    out[0] : list of strings
        These are the rows of the table.
    out[1] : string
        This is a string in the desired format which represents
        the table.

    """
    table = []
    for idx, result in enumerate(results):
        interf = result['out']['interfaces']
        # To set the L interface correctly for PPRETIS
        if idx > 1:
            interf0 = results[idx-1]['out']['interfaces'][1]
        else:
            interf0 = interf[0]
        detect = result['out']['detect']
        row = [
            f"{result['out']['ensemble']:^8s}",
            apply_format(interf0, '{:^8.4f}'),
            apply_format(interf[1], '{:^8.4f}'),
        ]
        intf2 = interf[2] if detect is None else detect
        row.append(apply_format(intf2, '{:^8.4f}'))
        table.append(row)
    if fmt in ['tex', 'latex']:
        table_str = generate_latex_table(table, 'Interfaces',
                                         ['Ensemble', 'Left', 'Middle',
                                          'Detect'],
                                         fixnum={0, 1, 2, 3})
    else:
        table_str = generate_rst_table(table, 'Interfaces',
                                       ['Ensemble', 'Left', 'Middle',
                                        'Detect'])
    table_txt = '\n'.join(table_str)
    return table, table_txt


def _table_interface(results, fmt='rst'):
    """Generate the table for the interfaces.

    This table will display the location of the different interfaces.

    Parameters
    ----------
    results : list of dicts
        These are the results of the analysis.
    fmt : string, optional
        Determines if we create reStructuredText ('rst') or
        latex ('tex').

    Returns
    -------
    out[0] : list of strings
        These are the rows of the table.
    out[1] : string
        This is a string in the desired format which represents
        the table.

    """
    table = []
    for result in results:
        interf = result['out']['interfaces']
        row = [
            f"{result['out']['ensemble']:^8s}",
            apply_format(interf[0], '{:^8.4f}'),
            apply_format(interf[1], '{:^8.4f}'),
        ]
        detect = result['out']['detect']
        if detect is None:
            row.append('')
        else:
            row.append(apply_format(detect, '{:^8.4f}'))
        table.append(row)
    if fmt in ['tex', 'latex']:
        table_str = generate_latex_table(table, 'Interfaces',
                                         ['Ensemble', 'Left', 'Middle',
                                          'Detect'],
                                         fixnum={0, 1, 2, 3})
    else:
        table_str = generate_rst_table(table, 'Interfaces',
                                       ['Ensemble', 'Left', 'Middle',
                                        'Detect'])
    table_txt = '\n'.join(table_str)
    return table, table_txt


def _table_interface0(results, fmt='rst'):
    """Generate the table for the [0^-] ensemble.

    This table will display the location of the different interfaces.

    Parameters
    ----------
    results : list of dicts
        These are the results of the analysis.
    fmt : string, optional
        Determines if we create reStructuredText ('rst') or
        latex ('tex').

    Returns
    -------
    out[0] : list of strings
        These are the rows of the table.
    out[1] : string
        This is a string in the desired format which represents
        the table.

    """
    table = []
    for result in results:
        row = [f"{result['out']['ensemble']:^8s}"]
        interf = result['out']['interfaces']
        row.append(apply_format(interf[0], '{:^8.4f}'))
        row.append(apply_format(interf[1], '{:^8.4f}'))
        row.append(apply_format(interf[2], '{:^8.4f}'))
        table.append(row)
    if fmt in ['tex', 'latex']:
        table_str = generate_latex_table(table, 'Interfaces',
                                         ['Ensemble', 'Left', 'Middle',
                                          'Right'],
                                         fixnum={0, 1, 2, 3})
    else:
        table_str = generate_rst_table(table, 'Interfaces',
                                       ['Ensemble', 'Left', 'Middle',
                                        'Right'])
    table_txt = '\n'.join(table_str)
    return table, table_txt


def _table_probability(results, fmt='rst'):
    """Generate the table for the probabilities.

    This table will display the crossing probabilities with
    uncertainties.

    Parameters
    ----------
    results : list of dicts
        The dictionaries are the results obtained from the analysis.
    fmt : string, optional
        Determines if we create reStructuredText ('rst') or
        latex ('tex').

    Returns
    -------
    out[0] : list of strings
        These are the rows of the table.
    out[1] : string
        This is a string in reStructuredText format which represents
        the table.

    """
    table = []
    for result in results:
        if 'prun' not in result['out']:
            continue
        row = [
            f"{result['out']['ensemble']:^8s}",
            apply_format(result['out']['prun'][-1], '{:^10.6f}'),
            apply_format(result['out']['blockerror'][2], '{:^10.6f}'),
            apply_format(result['out']['blockerror'][4] * 100, '{:^10.6f}'),
        ]
        table.append(row)
    if fmt in ['tex', 'latex']:
        table_str = generate_latex_table(table, r'Crossing probabilities',
                                         [r'Ensemble', r'$P_\text{cross}$',
                                          r'Error', r'Rel. error (\%)'],
                                         fixnum={0, 1, 2, 3})
    else:
        table_str = generate_rst_table(table, r'Crossing probabilities',
                                       [r'Ensemble', r'Pcross', r'Error',
                                        r'Rel. error (%)'])
    table_txt = '\n'.join(table_str)
    return table, table_txt


def _table_probability_repptis(results, fmt='rst'):
    """Generate the table for the probabilities.

    This table will display the crossing probabilities with
    uncertainties.

    Parameters
    ----------
    results : list of dicts
        The dictionaries are the results obtained from the analysis.
    fmt : string, optional
        Determines if we create reStructuredText ('rst') or
        latex ('tex').

    Returns
    -------
    out[0] : list of strings
        These are the rows of the table.
    out[1] : string
        This is a string in reStructuredText format which represents
        the table.

    """
    table = []
    for result in results:
        row = [
            f"{result['out']['ensemble']:^8s}",
            apply_format(result['out']['prun_sl'][-1], '{:^10.6f}'),
            apply_format(result['out']['prun_sr'][-1], '{:^10.6f}'),
            apply_format(result['out']['blockerror_sl'][2], '{:^10.6f}'),
            apply_format(result['out']['blockerror_sr'][2], '{:^10.6f}'),
        ]
        table.append(row)
    if fmt in ['tex', 'latex']:
        str1 = r'$\left\langle p^{ \pm}\right\rangle_{\text {full }}$'
        str2 = r'$\left\langle p^{ \mp}\right\rangle_{\text {full }}$'
        table_str = generate_latex_table(table, r'Probabilities',
                                         [r'Ensemble', str1, str2,
                                          r'Error$^{ \pm}$',
                                          r'Error$^{ \mp}$'],
                                         fixnum={0, 1, 2, 3, 4})
    else:
        table_str = generate_rst_table(table, r'Probabilities',
                                       [r'Ensemble', r'Pcross pm',
                                        r'Pcross mp', r'Error',
                                        r'Rel. error (%)'])
    table_txt = '\n'.join(table_str)
    return table, table_txt


def _table_efficiencies(results, fmt='rst'):
    """Generate table for efficiencies.

    This table will display results for the efficiencies, acceptance
    ratios and correlation.

    Parameters
    ----------
    results : list of dicts
        The dictionaries are the results obtained from the analysis.
    fmt : string, optional
        Determines if we create reStructuredText ('rst') or
        latex ('tex').

    Returns
    -------
    out[0] : list of strings
        These are the rows of the table.
    out[1] : string
        This is a string in reStructuredText format which represents
        the table.

    """
    table = []
    for result in results:
        if 'pathlength' not in result['out']:
            continue
        hist1 = result['out']['pathlength'][0]
        row = [
            f"{result['out']['ensemble']:^8s}",
            apply_format(result['out']['tis-cycles'], '{:^10.0f}'),
            apply_format(result['out']['efficiency'][0], '{:^10.6f}'),
            apply_format(result['out']['efficiency'][1], '{:^10.6f}'),
        ]
        row.append(apply_format(hist1[2][0], '{:^10.6f}'))
        table.append(row)
    if fmt in ['tex', 'latex']:
        table_str = generate_latex_table(table, 'Pathensemble data',
                                         ['Ensemble', 'TIS cycles',
                                          'Shoot acc. ratio',
                                          'Swap acc. ratio',
                                          'Avg. path length'],
                                         fixnum={0, 1, 2, 3, 4, 5})
    else:
        table_str = generate_rst_table(table, 'Pathensemble data',
                                       ['Ensemble', 'TIS cycles',
                                        'Shoot acc. ratio',
                                        'Swap acc. ratio',
                                        'Avg. path length'])
    table_txt = '\n'.join(table_str)
    return table, table_txt


def _table_summary(report, fmt='rst'):
    """Generate table with summary of main results.

    This table will display results for the crossing probability,
    the initial flux and rate constant.

    Parameters
    ----------
    report : dict
        Dict with current report. We use the information in this
        dictionary to generate the table.
    fmt : string, optional
        Determines if we create reStructuredText ('rst') or
        latex ('tex').

    Returns
    -------
    out[0] : list of strings
        These are the rows of the table.
    out[1] : string
        This is a string in reStructuredText format which represents
        the table.

    """
    numbers = report['numbers']
    text = report['text']
    table = [['Crossing probability', numbers['pcross'], numbers['perr']],
             [f"Flux (1/{text['flux-unit']})", numbers['flux'],
              numbers['flux-err']],
             [f"Rate constant (1/{text['flux-unit']})",
              numbers['rate'], numbers['rate-err']]]
    if fmt in ['tex', 'latex']:
        table_str = generate_latex_table(table, 'Summary of main results',
                                         ['Property', 'Value',
                                          r'Relative error ($\%$)'],
                                         fixnum={1, 2})
    else:
        table_str = generate_rst_table(table, 'Summary of main results',
                                       ['Property', 'Value',
                                        'Relative error (%)'])
    table_txt = '\n'.join(table_str)
    return table, table_txt


def _get_path_figures(figures):
    """Return path figures from a dict of figures.

    This method extracts figures from results and makes them available
    to the report.

    Parameters
    ----------
    figures : dict
        The figures generated by the analysis.

    Returns
    -------
    path_figures : dict
        A dict which can be used in the report.

    """
    path_figures = {}
    if not figures:
        return path_figures
    for fig in ('pcross', 'prun', 'perror', 'perror_sl',
                'perror_sr', 'lpath', 'prun_sl',
                'prun_sr', 'shoots'):
        for key in figures:
            if key.endswith(fig):
                path_figures[fig] = figures[key]
    return path_figures

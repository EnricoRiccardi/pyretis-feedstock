# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Functions for generating reports.

The reports are useful for displaying results from the analysis.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

generate_report_mdflux (:py:func:`.generate_report_mdflux`)
    Generate a report for a MD flux simulation.
"""
from pyretis.inout.report.markup import (generate_rst_table,
                                         generate_latex_table)
from pyretis.inout.formats.formatter import apply_format


__all__ = ['generate_report_mdflux']


def generate_report_mdflux(analysis, output='rst'):
    """Generate a report for a MD flux simulation.

    Parameters
    ----------
    analysis : dict
        This is the output from the analysis.
    output : string, optional
        This is the desired output format. It must match one of the
        formats defined in
        :py:const:`pyretis.inout.report.report.__TEMPLATES`. Default
        is 'rst' (reStructuredText).

    Returns
    -------
    out : dictionary
        The generated report as a string.

    """
    report = {'figures': {}, 'tables': {}}
    # generate some tables:
    cross_out = analysis['cross']['out']
    report['figures']['flux'] = [
        figure for figure in analysis['cross']['figures']
        if figure.get('runflux') and figure.get('errflux')
    ]
    report['figures']['energy'] = analysis['energy']['figures']
    report['figures']['order'] = analysis['order']['figures']
    report['tables']['md-flux'] = _table_md_flux(cross_out, fmt=output)[1]
    report['tables']['md-cycles'] = _table_md_flux_cycles(cross_out,
                                                          fmt=output)[1]
    report['tables']['md-efficiency'] = _table_md_efficiency(cross_out,
                                                             fmt=output)[1]
    return report


def _table_md_efficiency(results, fmt='rst'):
    """Generate table with MD-flux results for efficiencies and correlations.

    Parameters
    ----------
    results : dict
        These are the results obtained by :py:func:`.analyse_flux`.
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
    # table for interfaces:
    table = []
    for i, _ in enumerate(results['interfaces']):
        pmd = results['pMD'][i]
        teff = results['teffMD'][i]
        corr = results['corrMD'][i]
        prel = results['1-p'][i]
        row = [
            f'{i + 1:^10d}',
            apply_format(pmd, '{:^10.6f}'),
            apply_format(prel, '{:^10.6f}'),
            apply_format(teff, '{:^10.6f}'),
            apply_format(corr, '{:^10.6f}'),
        ]
        table.append(row)

    if fmt in ['tex', 'latex']:
        table_str = generate_latex_table(table, 'Efficiency',
                                         ['Interface',
                                          r'$p_\text{MD}$',
                                          r'$\frac{1-p}{p}$',
                                          'Efficiency time', 'Correlation'],
                                         fixnum={1, 2, 3, 4})
    elif fmt in ['txt']:
        table_str = generate_rst_table(table, 'Efficiency',
                                       ['Interface',
                                        'p_MD',
                                        '(1-p_MD)/p_MD',
                                        'Efficiency time', 'Correlation'])
    else:
        table_str = generate_rst_table(table, 'Efficiency',
                                       ['Interface',
                                        r':math:`p_\text{MD}`',
                                        r':math:`\frac{1-p}{p}`',
                                        'Efficiency time', 'Correlation'])
    table_txt = '\n'.join(table_str)
    return table, table_txt


def _table_md_flux_cycles(results, fmt='rst'):
    """Generate a table for the MD-flux results for cycle numbers.

    The table will display the number of steps in state A, state B,
    overall state A and B and the total number of MD cycles.

    Parameters
    ----------
    results : dict
        These are the results obtained by :py:func:`.analyse_flux`.
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
    # table for interfaces:
    table = [
        ['A', f'{results["times"]["A"]:8d}'],
        ['B', f'{results["times"]["B"]:8d}'],
        ['overall A', f'{results["times"]["OA"]:8d}'],
        ['overall B', f'{results["times"]["OB"]:8d}'],
        ['Total cycles', f'{results["totalcycle"]:8d}'],
    ]
    if fmt in ['tex', 'latex']:
        table_str = generate_latex_table(table, 'Cycles spent in state',
                                         ['State', 'Cycles'],
                                         fixnum={2})
    else:
        table_str = generate_rst_table(table, 'Cycles spent in state',
                                       ['State', 'Cycles'])
    table_txt = '\n'.join(table_str)
    return table, table_txt


def _table_md_flux(results, fmt='rst'):
    """Generate the table for the MD-flux results.

    Parameters
    ----------
    results : dict
        These are the results obtained by :py:func:`.analyse_flux`.
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
    # table for interfaces:
    table = []
    for i, idet in enumerate(results['interfaces']):
        flux = results['runflux'][i][-1]
        error = results['errflux'][i][2]
        relerror = results['errflux'][i][4]
        row = [
            f'{i + 1:^4d}',
            apply_format(idet, '{:^8.4f}'),
            apply_format(flux, '{:^10.6f}'),
            apply_format(error, '{:^10.6f}'),
            apply_format(relerror, '{:^10.6f}'),
            f'{results["ncross"][i]:^8d}',
            f'{results["neffcross"][i]:^8d}',
            apply_format(results['neffc/nc'][i], '{:^10.6f}'),
            apply_format(results['cross_time'][i], '{:^10.6f}'),
        ]
        table.append(row)
    if fmt in ['tex', 'latex']:
        table_str = generate_latex_table(table, 'Flux for interfaces',
                                         ['Int.', 'Position', 'Flux / units',
                                          'Error', 'Rel. error', 'Ncross',
                                          'Neffcross', 'Neffcross/Ncross',
                                          'Steps per cross'],
                                         fixnum={2, 3, 4, 5, 6, 7, 8})
    else:
        table_str = generate_rst_table(table, 'Flux for interfaces',
                                       ['Int.', 'Position', 'Flux / units',
                                        'Error', 'Rel. error', 'Ncross',
                                        'Neffcross', 'Neffcross/Ncross',
                                        'Steps per cross'])
    table_txt = '\n'.join(table_str)
    return table, table_txt

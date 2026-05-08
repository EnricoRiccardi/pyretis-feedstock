# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Some common functions for generating simple tables and numbers.

This module contains some common functions for the generation of
reports. The functions defined here are typically used to format
numbers and generate tables for the reports.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

generate_rst_table (:py:func:`.generate_rst_table`)
    Generate reStructuredText for a table.

generate_latex_table (:py:func:`.generate_latex_table`)
    Generate latex code for a table.

latexify_number (:py:func:`.latexify_number`)
    Change exponential notation into something nicer for latex.

"""


__all__ = ['generate_rst_table', 'generate_latex_table', 'latexify_number']


def generate_rst_table(table, title, headings):
    """Generate reStructuredText for a table.

    This is a general function to generate a table in reStructuredText.
    The table is specified with a title, headings for the columns and
    the contents of the columns and rows.

    Parameters
    ----------
    table : list of lists
        `table[i][j]` is the contents of column `j` of row `i` of the
        table.
    title : string
        The header/title for the table.
    headings : list of strings
        These are the headings for each table column.

    """
    # for each column, we need to figure out how wide it should be
    # 1) Check headings
    col_len = [len(col) + 2 for col in headings]
    # 2) Check if some of the columns are wider
    for row in table:
        for i, col in enumerate(row):
            if len(col) >= col_len[i] - 2:
                col_len[i] = len(col) + 2  # add some extra space

    width = len(col_len) + sum(col_len) - 1
    topline = '+' + width * '-' + '+'
    # create the header
    str_header = f'|{{0:<{width}s}}|'
    str_header = str_header.format(title)
    # make format for header:
    head_fmt = [f'{{0:^{col}s}}' for col in col_len]
    # create first row, which is the header on columns:
    row_line = [fmt.format(col) for fmt, col in zip(head_fmt, headings)]
    row_line = [''] + row_line + ['']
    row_line = '|'.join(row_line)
    # also set-up the horizontal line:
    hline = [''] + [col * '-' for col in col_len] + ['']
    hline = '+'.join(hline)
    # generate table
    str_table = [topline, str_header, hline, row_line, hline.replace('-', '=')]
    for row in table:
        row_line = [fmt.format(col) for fmt, col in zip(head_fmt, row)]
        row_line = [''] + row_line + ['']
        row_line = '|'.join(row_line)
        str_table.extend([row_line, hline])
    return str_table


def generate_latex_table(table, title, headings, fixnum=None):
    r"""Generate latex code for a table.

    This function will generate latex code for a table. The table is
    given with a title, headings for the columns and the contents of
    the table. For latex we might wish to make some numbers more pretty
    by removing exponential notation: i.e. ``1.e-10`` can be replaced
    by ``1.0 \times 10^{-10}`` (which should render like
    :math:`1.0 \times 10^{-10}`).

    Parameters
    ----------
    table : list of lists
        `table[i][j]` is the contents of column `j` of row `i` of the
        table.
    title : string
        The header/title for the table.
    headings : list of strings
        These are the headings for each table column.
    fixnum : list of integers, optional
        These integers identify the columns where `latexify_number` is
        to be applied.

    """
    str_table = [
        r'\begin{longtable}{' + len(headings) * '| c ' + '|}',
        r'\caption{' + title + r'}\\',
        r'\hline',
        r' & '.join(headings) + r'\\ \hline',
    ]
    for row in table:
        rowl = row
        if fixnum:
            rowl = [latexify_number(col) if i in fixnum else col for i, col
                    in enumerate(row)]
        str_table.append(' & '.join(rowl) + r'\\')
    str_table.append(r'\hline')
    str_table.append(r'\end{longtable}')
    return str_table


def latexify_number(str_float):
    r"""Change exponential notation into something nicer for latex.

    This will change exponential notation, e.g ``1.2e-03``, into
    ``1.2 \times 10^{-3}`` for latex output which should be rendered
    like :math:`1.2 \times 10^{-3}`.

    Parameters
    ----------
    str_float : string
        This is the string representation of a float.

    Returns
    -------
    out : string
        A formatted string for latex.

    """
    if 'e' in str_float:
        base, exp = str_float.split('e')
        return fr'${base} \times 10^{{{int(exp)}}}$'
    return fr'${str_float}$'

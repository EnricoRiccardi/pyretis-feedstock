# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Helper to create tables for PyRETIS keywords."""
import textwrap
import pathlib


def make_table(data, short_refs):
    """Make the actual table."""
    header1 = 'Keyword'
    maxlength = 79 - 3
    length = [len(i.split('|')[1]) + 2 for i in short_refs]
    length.append(len(header1))
    length1 = max(length)
    fmt1 = '| {{:<{}s}} |'.format(length1)
    border1 = '+' + '-' * (length1 + 2) + '+'
    eborder1 = '+' + '=' * (length1 + 2) + '+'
    rem = maxlength - len(border1)
    length2 = maxlength - length1 - 4 - 3
    fmt2 = ' {{:<{}s}} |'.format(length2)
    border2 = '-' * (rem - 1) + '+'
    eborder2 = '=' * (rem - 1) + '+'
    table = []
    table.append(border1 + border2)
    table.append(fmt1.format(header1) + fmt2.format('Description'))
    table.append(eborder1 + eborder2)

    for refi, datai in zip(short_refs, data):
        try:
            txt = (' '.join(datai)).split('>`:')[1].strip()
            txtl = textwrap.wrap(txt, width=length2)
        except IndexError:
            print('NO DATA FOUND!')
            txtl = ['']
        col1 = []
        col2 = []
        for i, txt in enumerate(txtl):
            if i == 0:
                col1.append(fmt1.format(refi.split()[1]))
            else:
                col1.append(fmt1.format(''))
            col2.append(fmt2.format(txt))
        for i, j in zip(col1, col2):
            table.append(i + j)
        table.append(border1 + border2)

    for ref in short_refs:
        print('\n{}'.format(ref))

    print()
    print('.. _table-NAMEOFTABLE:')
    print()
    print('.. table:: Supported XXX in |pyretis|')
    print('   :class: table-striped')
    print()

    for line in table:
        print('   {}'.format(line))


def shorten_ref(ref, add_to_ref):
    """Shorten references."""
    link = ':ref:{}'.format(ref.split(':ref:')[1].split(':')[0])
    linkr = link.split('<')[1].split('>')[0]
    short = '.. |{}_{}| replace:: {}'.format(add_to_ref,
                                             linkr.split('-')[-1],
                                             link)
    print(short.split('|')[1])
    return short


def read_keyword_section(filename, add_to_ref):
    """Read a keyword section in a .rst file."""
    read = False
    data = []
    sections = []
    section = []
    new_section = True
    filename = pathlib.Path(filename)
    with filename.open() as fileobj:
        for lines in fileobj:
            if lines.startswith('* :ref:'):
                if new_section:
                    if section:
                        sections.append(section)
                    section = []
                    new_section = False
                if data:
                    section.append(data)
                data = []
                read = True
            else:
                if lines.startswith('..'):
                    new_section = True
            if not lines.strip():
                read = False
            if read:
                data.append(lines.strip())
        if data:
            section.append(data)
        if section:
            sections.append(section)

    for section in sections:
        short = [shorten_ref(data[0], add_to_ref) for data in section]
        make_table(section, short)


if __name__ == '__main__':
    read_keyword_section('unitsystem.rst', 'units')

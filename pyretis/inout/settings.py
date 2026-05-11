# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This module handles parsing of input settings.

This module defines the file format for PyRETIS input files.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

parse_settings_file (:py:func:`.parse_settings_file`)
    Method for parsing settings from a given input file.

write_settings_file (:py:func:`.write_settings_file`)
    Method for writing settings from a simulation to a given file.
"""
import ast
import copy
import logging
import os
import pprint
import re
from pyretis.inout.restart import read_restart_file
from pyretis.inout.common import create_backup, create_empty_ensembles
from pyretis.inout.formats.cp2k import cp2k_settings
from pyretis.inout.formats.gromacs import gromacs_settings
from pyretis.info import PROGRAM_NAME, URL
from pyretis.inout.checker import check_interfaces, check_for_bullshitt
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())

__all__ = ['parse_settings_file', 'write_settings_file', 'SECTIONS',
           'RESTART_OVERRIDE_KEYWORDS', 'RESTART_FORBIDDEN_KEYWORDS',
           'fill_up_tis_and_retis_settings', 'add_default_settings',
           'add_specific_default_settings']


SECTIONS = {}
TITLE = f'{PROGRAM_NAME} input settings'
SECTIONS['heading'] = {
    'text': (f"{TITLE}\n{'=' * len(TITLE)}\n"
             f"For more info, please see: {URL}\n"
             f"Have Fun!")
}
HERE = os.path.abspath('.')

SECTIONS['simulation'] = {
    'endcycle': None,
    'exe_path': HERE,
    'flux': None,
    'interfaces': None,
    'restart': None,
    'rgen': 'rgen',
    'seed': None,
    'startcycle': None,
    'steps': None,
    'task': 'md',
    'zero_ensemble': None,
    'zero_left': None,
    'permeability': None,
    'swap_attributes': None,
    'priority_shooting': False,
    'umbrella': None,
    'overlap': None,
    'maxdx': None,
    'mincycle': None
}

SECTIONS['system'] = {
    'dimensions': 3,
    'input_file': None,
    'temperature': 1.0,
    'units': 'lj',
    # The random generator used for path-ensemble sampling is taken from
    # this section by pyretis.setup.createsimulation.create_ensemble. The
    # key is intentionally accepted here so that "rgen" in a [System]
    # block is honoured rather than silently dropped.
    'rgen': 'rgen',
}

SECTIONS['unit-system'] = {
    'charge': None,
    'energy': None,
    'length': None,
    'mass': None,
    'name': None,
}

SECTIONS['engine'] = {
    'class': None,
    'exe_path': HERE,
    'module': None,
    'rgen': 'rgen',
}

SECTIONS['box'] = {
    'cell': None,
    'high': None,
    'low': None,
    'periodic': None,
}

SECTIONS['particles'] = {
    'mass': None,
    'name': None,
    'npart': None,
    'position': None,
    'ptype': None,
    'type': 'internal',
    'velocity': None,
}

SECTIONS['forcefield'] = {
    'description': None
}

SECTIONS['potential'] = {
    'class': None,
    'parameter': None
}

SECTIONS['orderparameter'] = {
    'class': None,
    'module': None,
    'name': None,
}

SECTIONS['collective-variable'] = {
    'class': None,
    'module': None,
    'name': None
}

SECTIONS['output'] = {
    'backup': 'append',
    'cross-file': 1,
    'energy-file': 1,
    'pathensemble-file': 1,
    'prefix': None,
    'order-file': 1,
    'restart-file': 1,
    'screen': 10,
    'trajectory-file': 100,
}

SECTIONS['tis'] = {
    'allowmaxlength': False,
    'aimless': True,
    'ensemble_number': None,
    'detect': None,
    'freq': None,
    'maxlength': None,
    'nullmoves': None,
    'n_jumps': None,
    'high_accept': False,
    'interface_sour': None,
    'interface_cap': None,
    'relative_shoots': None,
    'rescale_energy': False,
    'rgen': 'rgen',
    'seed': None,
    'shooting_move': 'sh',
    'shooting_moves': [],
    'sigma_v': -1,
    'zero_momentum': False,
    'mirror_freq': 0,
    'target_freq': 0,
    'target_indices': [],
}

SECTIONS['initial-path'] = {
    'method': None
}

SECTIONS['retis'] = {
    'nullmoves': None,
    'relative_shoots': None,
    'rgen': None,
    'seed': None,
    'swapfreq': None,
    'swapsimul': None,
}

SECTIONS['repptis'] = {
    'memory': None
}

SECTIONS['ensemble'] = {
    'interface': None
}

SECTIONS['analysis'] = {
    'blockskip': 1,
    'bins': 100,
    'maxblock': 1000,
    'maxordermsd': -1,
    'ngrid': 1001,
    'plot': {'plotter': 'mpl', 'output': 'png',
             'style': 'pyretis'},
    'report': ['latex', 'rst', 'html'],
    'report-dir': None,
    'skipcross': 1000,
    'txt-output': 'txt.gz',
    'tau_ref_bin': [],
    'skip': 0
}


SPECIAL_KEY = {'parameter'}

# Section aliases: alternate section-header names that are silently mapped
# to a canonical section.  Listed here so they are not flagged as unknown.
# Keys are the alias (first word of a section header, lower-cased);
# values are the canonical SECTIONS key they merge into.
SECTION_ALIASES = {
    'md': 'simulation',
}

# This dictionary contains sections where the keywords
# can not be defined before we parse the input. The reason
# for this is that we support user-defined external modules
# and that the user should have the freedom to define keywords
# for these modules:
ALLOW_MULTIPLE = {
    'collective-variable',
    'engine',
    'ensemble',
    'initial-path',
    'orderparameter',
    'potential',
}

# This dictionary contains sections that can be defined
# multiple times. When parsing, these sections will be
# prefixed with a number to distinguish them.
SPECIAL_MULTIPLE = {
    'collective-variable',
    'ensemble',
    'potential',
}

# Keywords that the new input file is allowed to change on a restart.
# When a value in the new input differs from the value stored in the
# restart, a WARNING is emitted and the new value replaces the restarted
# one. A restart is meant to *extend* an existing simulation, so the keys
# here are limited to those that do not invalidate previously saved state
# (cycle counts, RNG seeds, log verbosity, … but not the interface set or
# any other topology-defining quantity — for those see
# RESTART_FORBIDDEN_KEYWORDS).
RESTART_OVERRIDE_KEYWORDS = {
    'simulation': [
        'steps',
        'seed',
        'priority_shooting',
    ],
    'output': [
        'screen',
        'backup',
        'trajectory-file',
        'energy-file',
        'order-file',
        'cross-file',
        'restart-file',
        'pathensemble-file',
    ],
    'tis': [
        # General shooting
        'freq',
        'sigma_v',
        'aimless',
        'allowmaxlength',
        'maxlength',
        'shooting_move',
        'shooting_moves',
        # wt / ss / wf specific
        'high_accept',
        'n_jumps',
        'mirror_freq',
        'target_freq',
        'target_indices',
        'nullmoves',
    ],
    'retis': [
        'swapfreq',
        'swapsimul',
        'nullmoves',
    ],
}


# Keywords that define the topology of the simulation (state space,
# ensemble layout, …). Changing any of these in a restart would
# invalidate the previously saved path data, so the restart must abort
# with a clear error if the new input disagrees with the stored value.
# Keys are section names; values are lists of keyword names.
RESTART_FORBIDDEN_KEYWORDS = {
    'simulation': [
        'interfaces',
        'task',
        'zero_left',
        'zero_ensemble',
        'flux',
        'permeability',
    ],
}


def parse_settings_file(filename, add_default=True):
    """Parse settings from a file name.

    Here, we read the file line-by-line and check if the current line
    contains a keyword, if so, we parse that keyword.

    Parameters
    ----------
    filename : string
        The file to parse.
    add_default : boolean
        If True, we will add default settings as well for keywords
        not found in the input.

    Returns
    -------
    settings : dict
        A dictionary with settings for PyRETIS.

    """
    with open(filename, 'r', encoding='utf-8') as fileh:
        raw_sections = _parse_sections(fileh)
    settings = _parse_all_raw_sections(raw_sections)
    if add_default:
        logger.debug('Adding default settings')
        add_default_settings(settings)
        add_specific_default_settings(settings)
    if settings['simulation']['task'] in {'retis', 'tis',
                                          'repptis', 'pptis',
                                          'explore', 'make-tis-files'}:
        fill_up_tis_and_retis_settings(settings)
        # Set up checks before to continue. This section shall GROW.
        check_interfaces(settings)
        check_for_bullshitt(settings)
    return _clean_settings(settings)


def parse_primitive(text):
    """Parse text to Python using the ast module.

    Parameters
    ----------
    text : string
        The text to parse.

    Returns
    -------
    out[0] : string, dict, list, boolean, or other type
        The parsed text.
    out[1] : boolean
        True if we managed to parse the text, False otherwise.

    """
    parsed = None
    success = False
    try:
        parsed = ast.literal_eval(text.strip())
        success = True
    except SyntaxError:
        parsed = text.strip()
        success = True
    except ValueError:
        parsed = text.strip()
        success = True
    return parsed, success


def look_for_keyword(line):
    """Search for a keyword in the given string.

    A string is assumed to define a keyword if the keyword appears as
    the first word in the string, ending with a `=`.

    Parameters
    ----------
    line : string
        A string to check for a keyword.

    Returns
    -------
    out[0] : string
        The matched keyword. It may contain spaces and it will also
        contain the matched `=` separator.
    out[1] : string
        A lower-case, stripped version of `out[0]`.
    out[2] : boolean
        `True` if we found a possible keyword.

    """
    # Match a word followed by a '=':
    key = re.match(r'(.*?)=', line)
    if key:
        keyword = ''.join([key.group(1), '='])
        keyword_low = key.group(1).strip().lower()
        for i in SPECIAL_KEY:
            if keyword_low.startswith(i):
                return keyword, i, True

        # Here we assume that keys with len One or Two are Atoms names
        if len(keyword_low) <= 2:
            keyword_low = key.group(1).strip()

        return keyword, keyword_low, True
    return None, None, False


def _parse_sections(inputtxt):
    """Find sections in the input file with raw data.

    This method will find sections in the input file and
    collect the corresponding raw data.

    Parameters
    ----------
    inputtxt : list of strings or iterable file object
        The raw data to parse.

    Returns
    -------
    raw_data : dict
        A dictionary with keys corresponding to the sections found
        in the input file. `raw_data[key]` contains the raw data
        for the section corresponding to `key`.

    """
    multiple = {key: 0 for key in SPECIAL_MULTIPLE}
    raw_data = {'heading': []}
    previous_line = None
    add_section = 'heading'
    data = []
    for lines in inputtxt:
        current_line, _, _ = lines.strip().partition('#')
        if not current_line:
            continue
        if current_line.startswith('---'):
            if previous_line is None:
                continue
            section_title = previous_line.split()[0].lower()
            if section_title in SPECIAL_MULTIPLE:
                new_section_title = f'{section_title}{multiple[section_title]}'
                multiple[section_title] += 1
                section_title = new_section_title
            if section_title not in raw_data:
                raw_data[section_title] = []
            raw_data[add_section].extend(data[:-1])
            data = []
            add_section = section_title
        else:
            data += [current_line]
        previous_line = current_line
    if add_section is not None:
        raw_data[add_section].extend(data)
    return raw_data


def _parse_section_heading(raw_section):
    """Parse the heading section.

    Parameters
    ----------
    raw_section : list of strings
        The text data for a given section which will be parsed.

    Returns
    -------
    setting : dict
        A dict with keys corresponding to the settings.

    """
    if not raw_section:
        return None
    return {'text': '\n'.join(raw_section)}


def _merge_section_text(raw_section):
    """Merge text for settings that are split across lines.

    This method supports keyword settings that are split across several
    lines. Here we merge these lines by assuming that keywords separate
    different settings.

    Parameters
    ----------
    raw_section : string
        The text we will merge.

    """
    merged = []
    for line in raw_section:
        _, _, found_keyword = look_for_keyword(line)
        if found_keyword or not merged:
            merged.append(line)
        else:
            merged[-1] = ''.join((merged[-1], line))
    return merged


def _parse_section_default(raw_section):
    """Parse a raw section.

    This is the default parser for sections.

    Parameters
    ----------
    raw_section : list of strings
        The text data for a given section which will be parsed.

    Returns
    -------
    setting : dict
        A dict with keys corresponding to the settings.

    """
    merged = _merge_section_text(raw_section)
    setting = {}
    for line in merged:
        match, keyword, found_keyword = look_for_keyword(line)
        if found_keyword:
            raw = line[len(match):].strip()
            parsed, success = parse_primitive(raw)
            if success:
                special = None
                for skey in SPECIAL_MULTIPLE:
                    # To avoid a false True for ensemble_number
                    if keyword.startswith(skey) and keyword[len(skey)] != '_':
                        special = skey

                if special is not None:
                    var = [''.join(line.split(keyword.split()[0])[1])]
                    new_setting = _parse_section_default(var)
                    var = line.split(special)[1].split()[0]
                    num = 0 if not var.isdigit() else int(var)

                    if special not in setting:
                        setting[special] = [{}]
                    while num >= len(setting[special]):
                        setting[special].append({})
                    setting[special][num].update(new_setting)

                elif keyword in SPECIAL_KEY:
                    if keyword not in setting:
                        setting[keyword] = {}
                    var = line.split(keyword)[1].split()[0]
                    # Yes, in some cases we really want an integer.
                    # Note: This will only work for positive numbers
                    # (which we are assuming here).
                    if var.isdigit():
                        setting[keyword][int(var)] = parsed
                    else:
                        setting[keyword][var] = parsed

                elif len(keyword.split()) > 1:
                    key_0 = match.split()[0]
                    var = [' '.join(line.split()[1:])]
                    new_setting = _parse_section_default(var)
                    if key_0 not in setting:
                        setting[key_0] = {}
                    for key, val in new_setting.items():
                        if key in setting[key_0]:
                            setting[key_0][key].update(val)
                        else:
                            setting[key_0][key] = val
                else:
                    setting[keyword] = parsed

            else:  # pragma: no cover
                msg = [f'Could read keyword {keyword}']
                msg += ['Keyword was skipped, please check your input!']
                msg += [f'Input setting: {raw}']
                msgtxt = '\n'.join(msg)
                logger.critical(msgtxt)
    return setting


def _parse_raw_section(raw_section, section):
    """Parse the raw data from a section.

    Parameters
    ----------
    raw_section : list of strings
        The text data for a given section which will be parsed.
    section : string
        A text identifying the section we are parsing for. This is
        used to get a list over valid keywords for the section.

    Returns
    -------
    out : dict
        A dict with keys corresponding to the settings.

    """
    if section in SECTION_ALIASES:
        # Parse against the canonical section's keyword list.
        canonical = SECTION_ALIASES[section]
        return _parse_raw_section(raw_section, canonical)
    if section not in SECTIONS:
        # Unknown section, just ignore it and give a warning.
        msgtxt = f'Ignoring unknown input section "{section}"'
        logger.warning(msgtxt)
        return None
    if section == 'heading':
        return _parse_section_heading(raw_section)
    return _parse_section_default(raw_section)


def _parse_all_raw_sections(raw_sections):
    """Parse all raw sections.

    This method is helpful for running tests etc.

    Parameters
    ----------
    raw_sections : dict
        The dictionary with the raw data in sections.

    Returns
    -------
    settings : dict
        The parsed settings, with one key for each section parsed.

    """
    settings = {}
    for key, val in raw_sections.items():
        special = None
        for i in SPECIAL_MULTIPLE:
            if key.startswith(i):
                special = i
        if special is not None:
            new_setting = _parse_raw_section(val, special)
            if special not in settings:
                settings[special] = []
            settings[special].append(new_setting)
        else:
            new_setting = _parse_raw_section(val, key)
            if new_setting is None:
                continue
            # Aliases (e.g. 'md') are merged into their canonical section.
            dest_key = SECTION_ALIASES.get(key, key)
            if dest_key not in settings:
                settings[dest_key] = {}
            for sub_key in new_setting:
                settings[dest_key][sub_key] = new_setting[sub_key]
    return settings


def fill_up_tis_and_retis_settings(settings):
    """Make the life of sloppy users easier.

    The full input set-up will be here completed.

    Parameters
    ----------
    settings : dict
        The current input settings.

    Returns
    -------
    None, but this method might add data to the input settings.

    """
    create_empty_ensembles(settings)
    ensemble_save = copy.deepcopy(settings['ensemble'])

    # The previously constructed dictionary is inserted in the settings.
    # This is done such that the specific input given per ensemble
    # OVERWRITES the general input.
    for i_ens, val in enumerate(ensemble_save):
        for key in settings:
            if key in val:
                if key not in SPECIAL_MULTIPLE:
                    val[key] = {**copy.deepcopy(settings[key]),
                                **copy.deepcopy(val[key])}
                else:
                    for i_sub, sub in enumerate(settings[key]):
                        while len(val[key]) < len(settings[key]):
                            val[key].append({})
                        val[key][i_sub] = {
                            **copy.deepcopy(sub),
                            **copy.deepcopy(val[key][i_sub])
                            }

        ensemble_save[i_ens] = {**copy.deepcopy(settings),
                                **copy.deepcopy(val)}
        del ensemble_save[i_ens]['ensemble']

    for i_ens, ens in enumerate(ensemble_save):
        add_default_settings(settings)
        add_specific_default_settings(settings)
        settings['ensemble'][i_ens] = copy.deepcopy(ens)
        if 'make-tis-files' in settings['simulation']['task']:
            settings['ensemble'][i_ens]['simulation']['task'] = 'tis'
    if settings['tis'].get('shooting_moves', False):
        for i_ens, ens_set in enumerate(settings['ensemble']):
            ens_set['tis']['shooting_move'] = \
                settings['tis']['shooting_moves'][i_ens]


def add_default_settings(settings):
    """Add default settings.

    Parameters
    ----------
    settings : dict
        The current input settings.

    Returns
    -------
    None, but this method might add data to the input settings.

    """
    if settings.get('initial-path', {}).get('method') == 'restart':
        if settings['simulation'].get('restart') is None:
            settings['simulation']['restart'] = 'pyretis.restart'

    for sec, sec_val in SECTIONS.items():
        if sec not in settings:
            settings[sec] = {}
        for key, val in sec_val.items():
            if val is not None and key not in settings[sec]:
                settings[sec][key] = val
    to_remove = [key for key in settings if len(settings[key]) == 0]
    for key in to_remove:
        settings.pop(key, None)


def add_specific_default_settings(settings):
    """Add specific default settings for each simulation task.

    Parameters
    ----------
    settings : dict
        The current input settings.

    Returns
    -------
    None, but this method might add data to the input settings.

    """
    task = settings['simulation'].get('task')
    if task not in settings and task in SECTIONS:
        settings[task] = {}

    if 'exp' in task:
        settings['tis']['shooting_move'] = 'exp'

    if task in {'pptis', 'tis', 'make-tis-files'}:
        if 'flux' not in settings['simulation']:
            settings['simulation']['flux'] = False
        if 'zero_ensemble' not in settings['simulation']:
            settings['simulation']['zero_ensemble'] = False

    if task in {'repptis', 'retis'}:
        if 'flux' not in settings['simulation']:
            settings['simulation']['flux'] = True
        if 'zero_ensemble' not in settings['simulation']:
            settings['simulation']['zero_ensemble'] = True

    eng_name = settings['engine'].get('class', 'NoneEngine')
    if eng_name[:7].lower() in {'gromacs', 'cp2k', 'lammps'}:
        settings['particles']['type'] = 'external'
        settings['engine']['type'] = 'external'
        input_path = os.path.join(settings['engine'].get('exe_path', '.'),
                                  settings['engine'].get('input_path', '.'))
        engine_checker = {'gromacs': gromacs_settings,
                          'cp2k': cp2k_settings}
        # Checks engine specific settings
        if engine_checker.get(eng_name[:7].lower()):
            engine_checker[eng_name[:7].lower()](settings, input_path)
    else:
        settings['particles']['type'] = 'internal'
        settings['engine']['type'] = settings['engine'].get('type', 'internal')


def _check_restart_forbidden(new_set, settings):
    """Raise if topology-defining keywords disagree between input and restart.

    Keywords listed in :data:`RESTART_FORBIDDEN_KEYWORDS` describe the
    state space of the simulation (interfaces, task, …). They must match
    the value stored in the restart, otherwise the already-saved paths
    would be reinterpreted against a different topology and the
    continuation would be silently wrong.

    Parameters
    ----------
    new_set : dict
        Settings loaded from the restart file.
    settings : dict
        Settings parsed from the new input file.

    Raises
    ------
    ValueError
        When the new input changes a forbidden keyword.

    """
    for section, keys in RESTART_FORBIDDEN_KEYWORDS.items():
        if section not in settings or section not in new_set:
            continue
        for key in keys:
            if key not in settings[section]:
                continue
            new_val = settings[section][key]
            old_val = new_set[section].get(key)
            if new_val != old_val:
                raise ValueError(
                    f'Restart refused: "{section}.{key}" differs between '
                    f'the new input ({new_val!r}) and the restart file '
                    f'({old_val!r}). This keyword is part of the '
                    f'simulation topology and cannot be changed by a '
                    f'restart — use "load" instead to start a new '
                    f'simulation with different settings.'
                )


def _apply_restart_overrides(new_set, settings):
    """Apply allowed keyword overrides from a new input onto restart settings.

    Topology-defining keywords (see :data:`RESTART_FORBIDDEN_KEYWORDS`)
    are first checked for equality; any disagreement aborts the restart.
    For every section and key listed in :data:`RESTART_OVERRIDE_KEYWORDS`,
    if the value supplied in the new input file differs from the one stored
    in the restart, a WARNING is logged and the new value replaces the old
    one.  Sections or keys absent from *settings* are silently skipped.

    Parameters
    ----------
    new_set : dict
        Settings loaded from the restart file (modified in-place).
    settings : dict
        Settings parsed from the new input file.

    """
    _check_restart_forbidden(new_set, settings)

    for section, keys in RESTART_OVERRIDE_KEYWORDS.items():
        if section not in settings or section not in new_set:
            continue
        for key in keys:
            if key not in settings[section]:
                continue
            new_val = settings[section][key]
            old_val = new_set[section].get(key)
            if new_val != old_val:
                logger.warning(
                    'Restart override: "%s.%s" changed from %r to %r.',
                    section, key, old_val, new_val,
                )
                new_set[section][key] = new_val

    for i, ens_new in enumerate(settings.get('ensemble', [])):
        if i >= len(new_set.get('ensemble', [])):
            break
        ens_restart = new_set['ensemble'][i]
        for section, keys in RESTART_OVERRIDE_KEYWORDS.items():
            if section not in ens_new or section not in ens_restart:
                continue
            for key in keys:
                if key not in ens_new[section]:
                    continue
                new_val = ens_new[section][key]
                old_val = ens_restart[section].get(key)
                if new_val != old_val:
                    logger.warning(
                        'Restart override (ensemble %i): "%s.%s" '
                        'changed from %r to %r.',
                        i, section, key, old_val, new_val,
                    )
                    ens_restart[section][key] = new_val


def settings_from_restart(settings):
    """Overwrite the settings with restart info.

    Here, we attempt to remove unwanted stuff from the input settings.

    NOTE: This structure doesn't allow modifications to a simulation
          with a restart. That is, restart ONLY extends one simulation.
          Load should be used for any other purpose.

    Keywords listed in :data:`RESTART_OVERRIDE_KEYWORDS` are an exception:
    when their value in the new input differs from the restarted value, a
    WARNING is emitted and the new value is applied.

    Parameters
    ----------
    settings : dict
        The current input settings that is going to be mostly overwritten.

    Returns
    -------
    new_set : dict
        The current input settings with the restart info.
    restart_info : dict
        The info to restart the various simulation objects.

    """
    cycle = settings['simulation']['steps']
    exe_path = settings['simulation']['exe_path']
    filename = os.path.join(settings['simulation']['exe_path'],
                            settings['simulation']['restart'])
    restart = read_restart_file(filename)

    if settings.get('initial-path', {}).get('flexible_restart', False):
        new_set = copy_settings(settings)
    else:
        new_set = restart.pop('settings')
        _apply_restart_overrides(new_set, settings)
    new_set['restart'] = restart
    new_set['simulation']['startcycle'] = new_set['simulation']['steps']
    new_set['simulation']['steps'] = cycle
    new_set['simulation']['restart'] = filename
    new_set['simulation']['exe_path'] = exe_path
    # This won't loop if it is an empty list
    for i in range(len(new_set.get('ensemble', []))):
        new_set['ensemble'][i]['simulation']['exe_path'] = exe_path
        new_set['ensemble'][i]['initial-path'] = {'method': 'restart'}
        new_set['ensemble'][i]['engine']['input_path'] = \
            settings['ensemble'][i]['engine'].get('input_path', '.')
    # Priority shooting aims is to recover a failed job, so it must be
    # allowed here to interfere to the initial settings.
    new_set['simulation']['priority_shooting'] = settings['simulation'].get(
        'priority_shooting', False)
    new_set['engine']['input_path'] = settings['engine'].get('input_path', '.')

    return new_set, restart


def look_for_input_files(input_path, required_files,
                         extra_files=None):
    """Check that required files for external engines are present.

    It will first search for the default files.
    If not present, it will search for the files with the
    same extension. In this search,
    if there are no files or multiple files for a required
    extension, the function will raise an Error.
    There might also be optional files which are not required, but
    might be passed in here. If these are not present we will
    not fail, but delete the reference to this file.

    Parameters
    ----------
    input_path : string
        The path to the folder where the input files are stored.
    required_files : dict of strings
        These are the file names types of the required files.
    extra_files : list of strings, optional
        These are the file names of the extra files.

    Returns
    -------
    out : dict
        The paths to the required and extra files we found.

    """
    if not os.path.isdir(input_path):
        msg = f'Input path folder {input_path} not existing'
        raise ValueError(msg)

    # Get the list of files in the input_path folder
    files_in_input_path = \
        [i.name for i in os.scandir(input_path) if i.is_file()]

    input_files = {}
    # Check if the required files are present
    for file_type, file_to_check in required_files.items():
        req_ext = os.path.splitext(file_to_check)[1][1:].lower()
        if file_to_check in files_in_input_path:
            input_files[file_type] = os.path.join(input_path, file_to_check)
            logger.debug('%s input: %s', file_type, input_files[file_type])
        else:
            # If not present, let's try to explore the folder by extension
            file_counter = 0
            selected_file = None
            for file_input in files_in_input_path:
                file_ext = os.path.splitext(file_input)[1][1:].lower()
                if req_ext == file_ext:
                    file_counter += 1
                    selected_file = file_input

            # Since we are guessing the correct files, give an error if
            # multiple entries are possible.
            if file_counter == 1:
                input_files[file_type] = os.path.join(input_path,
                                                      selected_file)
                logger.warning('using %s as "%s" file',
                               input_files[file_type], file_type)
            else:
                msg = f'Missing input file "{file_to_check}" '
                if file_counter > 1:
                    msg += f'and multiple files have extension ".{req_ext}"'
                raise ValueError(msg)

    # Check if the extra files are present. If so, add them to the input_files.
    # Gromacs engine takes a dictionary as extra_files, while cp2k takes a list
    # I'm not familiar with cp2k, so I'm assuming the list format is required.
    # Either way, both types have their merits, so I add a check for enginetype
    if extra_files:
        if isinstance(extra_files, dict):
            for file_type, file_to_check in extra_files.items():
                logger.debug('Checking for key %s and file %s', file_type,
                             file_to_check)
                if file_to_check in files_in_input_path:
                    logger.debug('Found %s', file_to_check)
                    input_files[file_type] = os.path.join(input_path,
                                                          file_to_check)
                else:
                    msg = f'Extra file {file_to_check} not present '
                    msg += f'in {input_path}'
                    logger.info(msg)
        elif isinstance(extra_files, list):
            input_files['extra_files'] = []
            for file_to_check in extra_files:
                if file_to_check in files_in_input_path:
                    input_files['extra_files'].append(file_to_check)
                else:
                    msg = f'Extra file {file_to_check} not present '
                    msg += f'in {input_path}'
                    logger.info(msg)
        else:
            msg = 'Extra files should be given in a dict or list format'
            msg += f', but got {type(extra_files)}'
            raise ValueError(msg)

    return input_files


def _clean_settings(settings):
    """Clean up input settings.

    Here, we attempt to remove unwanted stuff from the input settings.

    Parameters
    ----------
    settings : dict
        The current input settings.

    Returns
    -------
    settingsc : dict
        The cleaned input settings.

    """
    settingc = {}
    # Add other sections:
    for sec in settings:
        if sec not in SECTIONS:  # Well, ignore unknown ones:
            msgtxt = f'Ignoring unknown section "{sec}"'
            logger.warning(msgtxt)
            continue
        if sec in SPECIAL_MULTIPLE:
            settingc[sec] = list(settings[sec])
        else:
            settingc[sec] = {}
            if sec in ALLOW_MULTIPLE:  # Here, just add multiple sections:
                for key in settings[sec]:
                    settingc[sec][key] = settings[sec][key]
            else:
                for key in settings[sec]:
                    if key not in SECTIONS[sec]:  # Ignore junk:
                        msgtxt = f'Ignoring unknown "{key}" in "{sec}"'
                        logger.warning(msgtxt)
                    else:
                        settingc[sec][key] = settings[sec][key]
    to_remove = [key for key, val in settingc.items() if len(val) == 0]
    for key in to_remove:
        settingc.pop(key, None)
    return settingc


def settings_to_text(settings):
    """Turn settings into text usable for an output file.

    Parameters
    ----------
    settings : dict
        The dictionary to write

    Returns
    ------
    out : string
        Text representing the settings.

    """
    txt = []
    for section in SECTIONS:
        if section not in settings:
            continue
        if section in SPECIAL_MULTIPLE:
            for sec in settings[section]:
                title = section.capitalize()
                line = '-' * len(title)
                if section == 'ensemble':
                    _, raw_data = multiple_section_to_text(sec,
                                                           prefix=None,
                                                           pure=True)
                else:
                    raw_data = section_to_text(sec)
                txt.append(f'{title}\n{line}\n{raw_data}\n\n')
        elif section == 'heading':
            txt.append(f'{settings[section]["text"]}\n\n')
        else:
            if section in ('tis', 'retis'):
                title = f'{section.upper()} settings'
            else:
                title = f'{section.capitalize()} settings'
            line = '-' * len(title)
            raw_data = section_to_text(settings[section])
            txt.append(f'{title}\n{line}\n{raw_data}\n\n')
    return ''.join(txt)


def multiple_section_to_text(settings, prefix=None, pure=False):
    """Turn settings for the ensemble into text for output.

    Parameters
    ----------
    settings : dict
        A dictionary with settings to transform.
    prefix : string, optional
        If this string is given, it will be prepended to
        the setting we are writing.
    pure: boolean, optional
        The flag is used to track if subroutine works on a
        main section (True) or in a sub section (False).
        In the first case, prefix has to be re-set.

    Returns
    -------
    out[0] : string
        Formatted text representing the prefix to use in
        a recursive key-word search.
    out[1] : string
        Formatted text representing the settings.

    """
    data = []
    for key in settings:
        prefix = None if pure else prefix
        if key in SPECIAL_MULTIPLE:
            for i, entry in enumerate(settings[key]):
                temp_prefix = f'{key}{i:d}'
                _, txt = multiple_section_to_text(entry,
                                                  prefix=temp_prefix)
                data.append(txt)

        elif key == 'interface':
            pretty = pprint.pformat(settings[key], width=67)
            pretty = pretty.replace('\n', '\n' + ' ' * 67)
            txt = f'{key} = {pretty}'
            data.append(txt)

        elif key == 'heading':
            txt = f'{key} = {settings[key]}'
            data.append(txt)

        elif isinstance(settings[key], dict):
            base = prefix
            if prefix is None:
                prefix = key
            else:
                base = prefix
                prefix += ' ' + key
            _, txt = multiple_section_to_text(settings[key], prefix=prefix)
            prefix = base
            data.append(txt)

        else:
            txt = f'{prefix} {key} = {settings[key]}'
            data.append(txt)

    return prefix, '\n'.join(data)


def section_to_text(settings, prefix=None):
    """Turn settings for a section into text for output.

    Parameters
    ----------
    settings : dict
        A dictionary with settings to transform.
    prefix : string, optional
        If this string is given, it will be prepended to
        the setting we are writing.

    Returns
    -------
    out : string
        Formatted text representing the settings.

    """
    data = []
    for key in settings:
        if key == 'parameter':
            txt = section_to_text(settings[key], prefix='parameter')
        else:
            if prefix is not None:
                leng = len(str(key)) + 3 + len(prefix) + 1
            else:
                leng = len(str(key)) + 3
            pretty = pprint.pformat(settings[key], width=79-leng)
            pretty = pretty.replace('\n', '\n' + ' ' * leng)
            if prefix is not None:
                txt = f'{prefix} {key} = {pretty}'
            else:
                txt = f'{key} = {pretty}'
        if len(txt) >= 5:  # Shortest text, e.g: "a = 1".
            data.append(txt)
    return '\n'.join(data)


def write_settings_file(settings, outfile, backup=True):
    """Write simulation settings to an output file.

    This will write a dictionary to an output file in the PyRETIS input
    file format.

    Parameters
    ----------
    settings : dict
        The dictionary to write.
    outfile : string
        The file to create.
    backup : boolean, optional
        If True, we will backup existing files with the same file
        name as the provided file name.

    Note
    ----
    This will currently fail if objects have made it into the supplied
    ``settings``.

    """
    if backup:
        msg = create_backup(outfile)
        if msg:
            logger.info(msg)
    with open(outfile, 'w', encoding='utf-8') as fileh:
        txt = settings_to_text(settings)
        fileh.write(txt.strip())


def copy_settings(settings):
    """Return a copy of the given settings.

    Parameters
    ----------
    settings : dict of dicts
        A dictionary which we will return a copy of.

    Returns
    -------
    lsetting : dict of dicts
        A copy of the settings.

    """
    lsetting = {}
    for sec in settings:  # this is common for all simulations:
        lsetting[sec] = {}
        if sec in SPECIAL_MULTIPLE:
            lsetting[sec] = [copy.deepcopy(j) for j in settings[sec]]
        else:
            for key in settings[sec]:
                lsetting[sec][key] = settings[sec][key]
    return lsetting

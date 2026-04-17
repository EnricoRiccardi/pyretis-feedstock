# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Some method for interacting with input/output from CP2K.

Important methods defined here
------------------------------

update_cp2k_input (:py:func:`.update_cp2k_input`)
    A method for updating a CP2K input file and creating a new one.

read_cp2k_input (:py:func:`.read_cp2k_input`)
    A method to read a CP2K input file.
"""
import logging
import os
import numpy as np
from pyretis.core.box import box_matrix_to_list, box_vector_angles
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())

np.set_printoptions(legacy='1.21')


class SectionNode:
    """A class representing a section in the CP2K input.

    Attributes
    ----------
    title : string
        The title of the section
    parent : string
        The parent section if this node represents a
        sub-section.
    settings : list of strings
        The setting(s) for this particular node.
    data : string
        A section of settings if the node defines several
        settings.
    children : set of objects like :py:class:`.SectionNode`
        A set with the sub-sections of this section.
    level : integer
        An integer to remember how far down this node is.
        E.g. if the level is 2, this node is a sub-sub-section.
        This is used for printing.
    parents : list of strings or None
        A list representing the path from the node to the top
        section.

    """

    def __init__(self, title, parent, settings, data=None):
        """Initialise a node.

        Parameters
        ----------
        title : string
            The title of the section.
        parent : object like :py:class:`.SectionNode`
            The parent if this section is a sub-section.
        settings : list of strings
            The settings defined in this section.
        data : list of strings, optional
            A section of settings.

        """
        self.title = title
        self.parent = parent
        self.settings = settings
        if data:
            self.data = list(data)
        else:
            self.data = []
        self.children = set()
        self.level = 0
        self.parents = None

    def add_child(self, child):
        """Add a sub-section to the current section."""
        self.children.add(child)

    def get_all_parents(self):
        """Find the path to the top of the tree."""
        parents = [self.title]
        prev = self.parent
        while prev is not None:
            parents.append(prev.title)
            prev = prev.parent
        self.parents = parents[::-1]


def dfs_print(node, visited):
    """Walk through the nodes and print out text.

    Parameters
    ----------
    node : object like :py:class:`.SectionNode`
        The object representing a C2PK section.
    visited : set of objects like :py:class:`.SectionNode`
        The set contains the nodes we have already visited.

    Returns
    -------
    out : list of strings
        These strings represent the CP2K input file.

    """
    out = []
    pre = ' ' * (2 * node.level)
    if not node.settings:
        out.append(f'{pre}&{node.title}')
    else:
        out.append(f'{pre}&{node.title} {" ".join(node.settings)}')
    for lines in node.data:
        out.append(f'{pre}  {lines}')
    visited.add(node)
    for child in node.children:
        if child not in visited:
            for lines in dfs_print(child, visited):
                out.append(lines)
    out.append(f'{pre}&END {node.title}')
    return out


def set_parents(listofnodes):
    """Set parents for all nodes."""
    node_ref = {}

    def dfs_set(node, vis):
        """DFS traverse the nodes."""
        if node.parents is None:
            node.get_all_parents()
            par = '->'.join(node.parents)
            if par in node_ref:
                prev = node_ref.pop(par)
                par1 = f'{par}->{" ".join(prev.settings)}'
                par2 = f'{par}->{" ".join(node.settings)}'
                node_ref[par1] = prev
                node_ref[par2] = node
            else:
                node_ref[par] = node
        vis.add(node)
        for child in node.children:
            if child not in visited:
                dfs_set(child, vis)

    for nodes in listofnodes:
        visited = set()
        dfs_set(nodes, visited)
    return node_ref


def read_cp2k_input(filename):
    """Read a CP2K input file.

    Parameters
    ----------
    filename : string
        The file to open and read.

    Returns
    -------
    nodes : list of objects like :py:class:`.SectionNode`
        The root section nodes found in the file.

    """
    nodes = []
    current_node = None
    with open(filename, 'r', encoding="utf-8") as infile:
        for lines in infile:
            lstrip = lines.strip()
            if not lstrip:
                # skip empty lines
                continue
            if lstrip.startswith('&'):
                strip = lstrip[1:].split()
                if lstrip[1:].lower().startswith('end'):
                    current_node = current_node.parent
                else:
                    if len(strip) > 1:
                        setts = strip[1:]
                    else:
                        setts = []
                    new_node = SectionNode(strip[0].upper(),
                                           current_node, setts)
                    if current_node is None:
                        nodes.append(new_node)
                    else:
                        new_node.level = current_node.level + 1
                        current_node.add_child(new_node)
                    current_node = new_node
            else:
                if current_node is not None:
                    current_node.data.append(lstrip)
    return nodes


def _add_node(target, settings, data, nodes, node_ref):
    """Just add a new node."""
    # check if this is a root node:
    root = target.find('->') == -1
    if root:
        new_node = SectionNode(target, None, settings, data=data)
        nodes.append(new_node)
    else:
        parents = target.split('->')
        title = parents[-1]
        par = '->'.join(parents[:-1])
        if par not in node_ref:
            _add_node(par, None, None, nodes, node_ref)
        parent = node_ref['->'.join(parents[:-1])]
        new_node = SectionNode(title, parent, settings, data=data)
        new_node.level = parent.level + 1
        parent.add_child(new_node)
    node_ref[target] = new_node


def update_node(target, settings, data, node_ref, nodes,
                replace=False):
    """Update the given target node.

    If the node does not exist, it will be created.

    Parameters
    ----------
    target : string
        The target node, on form root->section->subsection
    settings : list of strings
        The settings for the node.
    data : list of strings or dict
        The data for the node.
    node_ref : dict of :py:class:`.SectionNode`
        A dict of all nodes in the tree.
    nodes : list of :py:class:`.SectionNode`
        The root nodes.
    replace : boolean, optional
        If this is True and if the nodes have some data, the already
        existing data will be ignored. We also assume that the data
        is already formatted.

    """
    if target not in node_ref:  # add node
        # TODO: remove decommented try-except construction later
        # try:
        _add_node(target, settings, data, nodes, node_ref)
        # except KeyError:
        #    pass
        return None
    node = node_ref[target]
    new_data = []
    done = set()
    if not replace:
        for line in node.data:
            key = line.split()[0]
            if key in data:
                new_data.append(f'{key} {data[key]}')
                done.add(key)
            else:
                new_data.append(line)
        for key in data:
            if key in done:
                continue
            if data[key] is None:
                new_data.append(str(key))
            else:
                new_data.append(f'{key} {data[key]}')
        node.data = list(new_data)
    else:
        node.data = list(data)
    if settings is not None:
        if replace:
            node.settings = list(settings)
        else:
            node.settings += settings
    return node


def remove_node(target, node_ref, root_nodes):
    """Remove a node (and it's children) from the tree.

    Parameters
    ----------
    target : string
        The target node, on form root->section->subcection.
    node_ref : dict
        A dict with all the nodes.
    root_nodes : list of objects like :py:class:`.SectionNode`
        The root nodes.

    """
    to_del = node_ref.pop(target, None)
    if to_del is None:
        pass
    else:
        # remove all it's children:
        visited = set()
        nodes = [to_del]
        while nodes:
            node = nodes.pop()
            if node not in visited:
                visited.add(node)
                for i in node.children:
                    nodes.append(i)
        # remove the reference to this node from the parent
        parent = to_del.parent
        if parent is None:
            # This is a root node.
            root_nodes.remove(to_del)
        else:
            parent.children.remove(to_del)
        del to_del
        for key in visited:
            _ = node_ref.pop(key, None)


def update_cp2k_input(template, output, update=None, remove=None):
    """Read a template input and create a new CP2K input.

    Parameters
    ----------
    template : string
        The CP2K input file we use as a template.
    output : string
        The CP2K input file we will create.
    update : dict, optional
        The settings we will update.
    remove : list of strings, optional
        The nodes we will remove.

    """
    nodes = read_cp2k_input(template)
    node_ref = set_parents(nodes)
    if update is not None:
        for target in update:
            value = update[target]
            settings = value.get('settings', None)
            replace = value.get('replace', False)
            data = value.get('data', [])
            update_node(target, settings, data, node_ref, nodes,
                        replace=replace)
    if remove is not None:
        for nodei in remove:
            remove_node(nodei, node_ref, nodes)
    with open(output, 'w', encoding='utf-8') as outf:
        for i, nodei in enumerate(nodes):
            vis = set()
            if i > 0:
                outf.write('\n')
            outf.write('\n'.join(dfs_print(nodei, vis)))
            outf.write('\n')


def read_box_data(box_data):
    """Read the box data.

    Parameters
    ----------
    box_data : list of strings
        The settings for the SUBSYS->CELL section.

    Returns
    -------
    out[0] : numpy.array, 1D
        The box vectors, in the correct order.
    out[1] : list of booleans
        The periodic boundary setting for each dimension.

    """
    to_read = {'A': 'vec', 'B': 'vec', 'C': 'vec', 'PERIODIC': 'string',
               'ABC': 'vec', 'ALPHA_BETA_GAMMA': 'vec'}
    data = {}
    for lines in box_data:
        for key, val in to_read.items():
            keyword = f'{key} '
            if lines.startswith(keyword):
                if val == 'vec':
                    data[key] = [float(i) for i in lines.split()[1:]]
                elif val == 'string':
                    data[key] = ' '.join(lines.split()[1:])
    if all(('A' in data, 'B' in data, 'C' in data)):
        box_matrix = np.zeros((3, 3))
        box_matrix[:, 0] = data['A']
        box_matrix[:, 1] = data['B']
        box_matrix[:, 2] = data['C']
        box = box_matrix_to_list(box_matrix)
    elif 'ABC' in data:
        if 'ALPHA_BETA_GAMMA' in data:
            box_matrix = box_vector_angles(
                data['ABC'],
                data['ALPHA_BETA_GAMMA'][0],
                data['ALPHA_BETA_GAMMA'][1],
                data['ALPHA_BETA_GAMMA'][2],
            )
            box = box_matrix_to_list(box_matrix)
        else:
            box = np.array(data['ABC'])
    else:
        box = None
    periodic = []
    periodic_setting = data.get('PERIODIC', 'XYZ')
    for val in ('X', 'Y', 'Z'):
        periodic.append(val in periodic_setting.upper())
    return box, periodic


def read_cp2k_energy(energy_file):
    """Read and return CP2K energies.

    Parameters
    ----------
    energy_file : string
        The input file to read.

    Returns
    -------
    out : dict
        This dict contains the energy terms read from the CP2K energy file.

    """
    data = np.loadtxt(energy_file)
    energy = {}
    for i, key in ((1, 'time'), (2, 'ekin'), (3, 'temp'), (4, 'vpot')):
        try:
            energy[key] = data[:, i]
        except IndexError:
            logger.warning('Could not read energy term %s from CP2kfile %s',
                           key, energy_file)
    if 'ekin' in energy and 'vpot' in energy:
        energy['etot'] = energy['ekin'] + energy['vpot']
    return energy


def read_cp2k_restart(restart_file):
    """Read some info from a CP2K restart file.

    Parameters
    ----------
    restart_file : string
        The file to read.

    Returns
    -------
    pos : numpy.array
        The positions.
    vel : numpy.array
        The velocities.
    box_size : numpy.array
        The box vectors.
    periodic : list of booleans
        For each dimension, the list entry is True if periodic
        boundaries should be applied.

    """
    nodes = read_cp2k_input(restart_file)
    node_ref = set_parents(nodes)
    velocity = 'FORCE_EVAL->SUBSYS->VELOCITY'
    coord = 'FORCE_EVAL->SUBSYS->COORD'
    cell = 'FORCE_EVAL->SUBSYS->CELL'

    atoms, pos, vel = [], [], []

    for posi, veli in zip(node_ref[coord].data, node_ref[velocity].data):
        pos_split = posi.split()
        atoms.append(pos_split[0])
        pos.append([float(i) for i in pos_split[1:4]])
        vel.append([float(i) for i in veli.split()])
    pos = np.array(pos)
    vel = np.array(vel)
    box, periodic = read_box_data(node_ref[cell].data)
    return atoms, pos, vel, box, periodic


def read_cp2k_box(inputfile):
    """Read the box from a CP2K file.

    Parameters
    ----------
    inputfile : string
        The file we will read from.

    Returns
    -------
    out[0] : numpy.array
        The box vectors.
    out[1] : list of booleans
        For each dimension, the list entry is True if periodic
        boundaries should be applied.

    """
    nodes = read_cp2k_input(inputfile)
    node_ref = set_parents(nodes)
    try:
        box, periodic = read_box_data(
            node_ref['FORCE_EVAL->SUBSYS->CELL'].data
        )
    except KeyError:
        logger.warning('No CELL found in CP2K file "%s"', inputfile)
        box = np.array([[100., 0., 0.], [0., 100., 0.], [0., 0., 100.]])
        periodic = [True, True, True]
    return box, periodic


def cp2k_settings(settings, input_path):
    """Read and processes cp2k settings.

    Parameters
    ----------
    settings : dict
        The current input settings..
    input_path : string
        The CP2K input path

    """
    ext = settings['engine'].get('cp2k_format', 'xyz')
    default_files = {'conf': f'initial.{ext}',
                     'template': 'cp2k.inp'}
    settings['engine']['input_files'] = {}
    for key in ('conf', 'template'):
        # Add input path and the input files if input is no given:
        settings['engine']['input_files'][key] = \
            settings['engine'].get(key,
                                   os.path.join(input_path,
                                                default_files[key]))
    nodes = read_cp2k_input(
        settings['engine']['input_files']['template'])
    md_data = set_parents(nodes)['MOTION->MD'].data

    # Checks temperature
    cp2k_temp_inp = False
    for lines in md_data:
        if lines.startswith('TEMPERATURE'):
            cp2k_temp_inp = True
            cp2k_temp = float(lines.split()[1])
            if 'temperature' not in settings['system']:
                settings['system']['temperature'] = cp2k_temp
            elif abs(cp2k_temp - settings['system']['temperature']) < 10e-8:
                pass
            else:
                msg = 'Inequal Pyretis and CP2K input temperatures!' +\
                        ' Temperature is set to CP2K input temperature.'
                logger.warning(msg)
                settings['system']['temperature'] = cp2k_temp
            break
    if not cp2k_temp_inp:
        if 'temperature' not in settings['system']:
            msg = 'Temperature not set in CP2K file!' +\
                ' Temperature is set to CP2K default temperature.'
            logger.warning(msg)
        elif abs(300.0 - settings['system']['temperature']) > 10e-8:
            msg = 'Temperature not set in CP2K file!' +\
                ' And temperature in input rst file does not equal 300.0 K!' +\
                ' Temperature is set to CP2K default temperature.'
            logger.warning(msg)
        # CP2K defaults to 300K when temperature is unspecified.
        settings['system']['temperature'] = 300.0

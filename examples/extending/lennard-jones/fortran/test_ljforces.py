# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the FORTRAN implementation of the Lennard Jones potential.

This test is comparing the three versions of the Lennard Jones
potential:

1) The pure python implementation

2) The numpy python implementation

3) The FORTRAN implementation.
"""
# pylint: disable=invalid-name
import numpy as np
from pyretis.core import System, create_box, Particles
from pyretis.core.units import create_conversion_factors
from pyretis.forcefield import ForceField
from pyretis.forcefield.potentials import PairLennardJonesCut
from pyretis.forcefield.potentials import PairLennardJonesCutnp
from pyretis.tools import generate_lattice
from ljpotentialf import PairLennardJonesCutF


def set_up_initial_state():
    """Create particles for the test."""
    create_conversion_factors('lj')
    lattice, size = generate_lattice('fcc', [5, 5, 5], density=0.9)
    npart = len(lattice)
    lattice += np.random.randn(npart, 3) * 0.05
    size = np.array(size)
    box = create_box(low=size[:, 0], high=size[:, 1],
                     periodic=[True, True, True])
    system = System(temperature=1.0, units='lj', box=box)
    system.particles = Particles(dim=3)
    for pos in lattice:
        system.add_particle(name='Ar', pos=pos, mass=1.0, ptype=0)
    msg = 'Created lattice with {} atoms.'
    print(msg.format(system.particles.npart))
    return system


def run_calculations(system, parameters):
    """Evaluate the LJ potential."""
    # Calculate with FORTRAN:
    potential_ext = PairLennardJonesCutF(dim=3, shift=True,
                                         mixing='geometric')
    forcefield_ext = ForceField('Python with external FORTRAN force field',
                                potential=[potential_ext],
                                params=[parameters])
    system.forcefield = forcefield_ext
    print(f'Evaluating with: {forcefield_ext.print_potentials()}')
    vpot_ext, forces_ext, virial_ext = system.potential_and_force()
    vpot_ext /= float(system.particles.npart)
    # Calculate with pure python implementation:
    potential = PairLennardJonesCut(dim=3, shift=True)
    forcefield = ForceField('Pure Python force field',
                            potential=[potential],
                            params=[parameters])
    system.forcefield = forcefield
    print(f'Evaluating with: {forcefield.print_potentials()}')
    vpot, forces, virial = system.potential_and_force()
    vpot /= float(system.particles.npart)
    # Calculate with numpy python implementation:
    potentialnp = PairLennardJonesCutnp(dim=3, shift=True)
    forcefieldnp = ForceField('Python force field with numpy',
                              potential=[potentialnp],
                              params=[parameters])
    system.forcefield = forcefieldnp
    print(f'Evaluating with: {forcefieldnp.print_potentials()}')
    vpotnp, forcesnp, virialnp = system.potential_and_force()
    vpotnp /= float(system.particles.npart)
    return ((vpot, forces, virial),
            (vpotnp, forcesnp, virialnp),
            (vpot_ext, forces_ext, virial_ext))


def _assert_results_close(result, keys):
    """Assert that all implementations give the same result."""
    for i, keyi in enumerate(keys[:-1]):
        for j, key2 in enumerate(keys[i+1:]):
            print(f'\nCompare {keyi} and {key2}')
            force = np.allclose(result[i][1], result[i+j+1][1])
            print(f' -> Forces close: {force}')
            assert force
            virial = np.allclose(result[i][2], result[i+j+1][2])
            print(f' -> Virial close: {virial}')
            assert virial
            assert abs(result[i][0] - result[i+j+1][0]) < 1e-7
            vdiff = np.abs(result[i][0] - result[i+j+1][0])
            print(f' -> Difference in pot. energy: {vdiff:.15e}')


class TestLennardJonesFortran:
    """Run the tests for the FORTRAN potential class."""

    def test_lj(self):
        """Test one-component system."""
        print('\nTesting for a one-component system')
        system = set_up_initial_state()
        param = {0: {'sigma': 1.0, 'epsilon': 1.0, 'rcut': 2.5}}
        maxcut = 0.5 * min(system.box.length)
        for key in param:
            if 'rcut' in param[key]:
                assert maxcut >= param[key]['rcut']
        result = run_calculations(system, param)
        _assert_results_close(result, ['python', 'python-numpy', 'fortran'])

    def test_lj_mix(self):
        """Test for mixture."""
        print('\nTesting for a two-component mixture')
        system = set_up_initial_state()
        param = {0: {'sigma': 1.0, 'epsilon': 1.0, 'rcut': 2.5},
                 1: {'sigma': 2.0, 'epsilon': 1.2, 'rcut': 3.5}}
        idx2 = np.random.choice(list(range(system.particles.npart)),
                                size=int(system.particles.npart * 0.5),
                                replace=False)
        maxcut = 0.5 * min(system.box.length)
        for key in param:
            if 'rcut' in param[key]:
                assert maxcut >= param[key]['rcut']
        print(f'Mutating {len(idx2)} particles')
        for i in idx2:
            system.particles.ptype[i] = 1
        result = run_calculations(system, param)
        _assert_results_close(result, ['python', 'python-numpy', 'fortran'])

    def test_lj_multi_mix(self):
        """Test for multi-mixture."""
        ncomp = np.random.randint(3, 11)
        print(f'\nTesting for a {ncomp}-component mixture')
        system = set_up_initial_state()
        param = {0: {'sigma': 1.0, 'epsilon': 1.0, 'rcut': 2.5}}
        maxcut = 0.5 * min(system.box.length)
        assert maxcut >= param[0]['rcut']
        idx = np.array(list(range(system.particles.npart)), dtype='i4')
        np.random.shuffle(idx)
        for i, idx2 in enumerate(np.array_split(idx, ncomp)):
            system.particles.ptype[idx2] = i
            if i not in param:
                param[i] = {'sigma': np.random.uniform(low=0.5, high=1.5),
                            'epsilon': np.random.uniform(low=0.5, high=1.5),
                            'rcut': np.random.uniform(low=2.0, high=maxcut)}
        natoms = {}
        for i in range(system.particles.npart):
            ptype = system.particles.ptype[i]
            if ptype not in natoms:
                natoms[ptype] = 0
            natoms[ptype] += 1
        for atom in natoms:
            print(f'{natoms[atom]} atoms of type {atom}')
        result = run_calculations(system, param)
        _assert_results_close(result, ['python', 'python-numpy', 'fortran'])

# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Some common methods for the engine tests."""
import os
from pyretis.core.system import System
from pyretis.core.particles import ParticlesExt
from pyretis.engines.openmm import HAS_OPENMM
if HAS_OPENMM:
    try:
        import openmm as mm
        from openmm import app
        from openmm import unit
        from simtk import unit
    except ImportError:  # Openmm < 7.6
        from simtk.openmm import app
        import simtk.openmm as mm
        from simtk import unit
else:
    mm, app, unit = None, None, None

PDB_STRING = (
    """CRYST1   68.478   68.478   68.472  60.00  60.00  90.00 P 1           1
MODEL        0
ATOM      1  O   HOH A   1      38.510  67.082  29.996  1.00  0.00           O
ATOM      2  H1  HOH A   1      39.260  66.728  29.518  1.00  0.00           H
ATOM      3  H2  HOH A   1      38.877  67.385  30.827  1.00  0.00           H
ATOM      4  O   HOH A   2      19.361  51.482   4.207  1.00  0.00           O
ATOM      5  H1  HOH A   2      19.004  51.112   5.014  1.00  0.00           H
ATOM      6  H2  HOH A   2      20.299  51.294   4.252  1.00  0.00           H
ENDMDL
END
""")


def remove_dir(dirname):
    """Remove a directory silently."""
    try:
        os.removedirs(dirname)
    except OSError:
        pass


def make_test_system(conf):
    """Make a system with particles for testing."""
    system = System()
    system.particles = ParticlesExt(dim=3)
    system.particles.set_pos(conf)
    system.particles.set_vel(False)
    return system


def create_openmm_simulation(pdb='test_openmm_pyretis.pdb'):
    """
    This creates a default OpenMM simulation object for testing purposses.
    Adapted from builder.openmm.org
    """
    pdb = app.PDBFile(pdb)
    forcefield = app.ForceField('amber99sbildn.xml', 'tip3p.xml')

    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME,
                                     nonbondedCutoff=1.0*unit.nanometers,
                                     constraints=app.HBonds, rigidWater=True,
                                     ewaldErrorTolerance=0.0005)
    integrator = mm.LangevinIntegrator(300*unit.kelvin,
                                       1.0/unit.picoseconds,
                                       2.0*unit.femtoseconds)
    integrator.setConstraintTolerance(0.00001)

    simulation = app.Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    return simulation


def write_openmm_simulation(name, pdb='test_openmm_pyretis.pdb'):
    CODE = (
        f"""
from openmm import app
import openmm as mm
from openmm import unit

pdb = app.PDBFile("{pdb}")
forcefield = app.ForceField('amber99sbildn.xml', 'tip3p.xml')

system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME,
                                 nonbondedCutoff=1.0*unit.nanometers,
                                 constraints=app.HBonds, rigidWater=True,
                                 ewaldErrorTolerance=0.0005)
integrator = mm.LangevinIntegrator(300*unit.kelvin,
                                   1.0/unit.picoseconds,
                                   2.0*unit.femtoseconds)
integrator.setConstraintTolerance(0.00001)

simulation = app.Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
""")
    with open(name, 'w', encoding='utf-8') as output:
        output.write(CODE)


def write_test_pdb(name):
    """
    This writes test.pdb with 2 waters into the current running directory.
    """
    with open(name, 'w', encoding='utf-8') as output:
        output.write(PDB_STRING)


class FakeOp:
    """Order parameter for testing the OpenMMEngine."""
    def __init__(self, values=None):
        self.values = values
        self.n = 0

    def calculate(self, _):
        """Calculate the fake order parameter."""
        value = self.values[self.n]
        self.n += 1
        return [value]

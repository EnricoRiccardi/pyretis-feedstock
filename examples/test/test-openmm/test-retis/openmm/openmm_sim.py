# Minimal 2-Argon-atom Lennard-Jones simulation for PyRETIS OpenMM testing.
#
# This module is loaded by the PyRETIS OpenMM engine via the
# openmm_module / openmm_simulation settings in retis.rst.
# It defines `simulation`, an app.Simulation object that PyRETIS uses to
# propagate the system.
#
# The integrator seed is fixed so that results are reproducible across runs.
try:
    import openmm as mm
    from openmm import app, unit
except ImportError:
    import simtk.openmm as mm  # legacy namespace
    from simtk.openmm import app
    from simtk import unit

import numpy as np

# ── System: two Argon atoms ──────────────────────────────────────────────────
system = mm.System()
for _ in range(2):
    system.addParticle(39.948 * unit.dalton)

# Lennard-Jones parameters for Ar (no electrostatics)
lj = mm.NonbondedForce()
lj.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
for _ in range(2):
    lj.addParticle(0.0,
                   0.34 * unit.nanometer,
                   0.997 * unit.kilojoule_per_mole)
system.addForce(lj)

# ── Topology ─────────────────────────────────────────────────────────────────
topology = app.Topology()
chain = topology.addChain()
residue = topology.addResidue('AR', chain)
argon = app.Element.getByAtomicNumber(18)
for atom_name in ('Ar1', 'Ar2'):
    topology.addAtom(atom_name, argon, residue)

# ── Integrator & Simulation ──────────────────────────────────────────────────
integrator = mm.LangevinIntegrator(
    300 * unit.kelvin,
    1.0 / unit.picoseconds,
    2.0 * unit.femtoseconds,
)
# Fixed seed ensures deterministic thermal noise during propagation steps.
integrator.setRandomNumberSeed(42)

simulation = app.Simulation(topology, system, integrator)

# Initial positions: atoms placed near the LJ minimum (~0.38 nm)
positions = np.array([[0.0, 0.0, 0.0],
                      [0.38, 0.0, 0.0]]) * unit.nanometer
simulation.context.setPositions(positions)
simulation.context.setVelocitiesToTemperature(300 * unit.kelvin, 42)

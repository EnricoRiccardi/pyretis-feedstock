Testing the LAMMPS engine by running a RETIS simulation
========================================================

This test will run a RETIS simulation with the LAMMPS engine
and compare the results to previously obtained reference results.

The system is a 3D Lennard-Jones fluid (512 atoms, units lj) with
the x-coordinate of atom 1 used as the order parameter.
Interfaces are set at [0.5, 1.5, 3.5] in reduced LJ units.

The folder is:

* ``lammps1``: Which contains a test using the PyRETIS ``LAMMPSEngine``.

Instructions
------------

The test is run using the ``run.sh`` script. You can provide the
LAMMPS executable as an argument, e.g. ``./run.sh lmp_mpi``. The
default is ``lmp_serial``.

Before comparing, reference results must be present in
``lammps1/results.tgz``. To generate them, run the simulation once
and then archive the output::

    cd lammps1
    pyretisrun -i retis.rst -p -l DEBUG
    tar czf results.tgz results/
    cd ..

Commit ``lammps1/results.tgz`` to version control. Subsequent calls
to ``./run.sh`` will compare against these reference results.

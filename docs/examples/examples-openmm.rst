.. _examples-openmm:

Running a |pyretis| simulation with OpenMM
==========================================

In this example, we show the interface between ``OpenMM`` and |pyretis|.

First, we need to generate a file with an OpenMM ``Simulation`` object. We
use the online tool `OpenMM Script Builder <http://builder.openmm.org/>`_ to
generate the example with this
:download:`PDB file </_static/engine-examples/input.pdb>` containing two
water molecules.

.. literalinclude:: /_static/engine-examples/openmm_sim.py


With this ``simulation`` object, we can now construct the required ``.rst``
options in the engine section:

.. literalinclude:: /_static/engine-examples/openmm_retis.rst    
   :language: rst
   :lines: 14-20

Here we say that:

  * The engine ``type`` is ``openmm``, which tells |pyretis| that it has a
    special system representation.

  * The engine ``class`` is also ``openmm`` to tell |pyretis| we want to use the
    ``OpenMM`` engine.

  * The ``openmm_module`` is the module from which we want to load the
    simulation. Here it is ``openmm_sim.py``.

  * The ``openmm_simulation`` is the name of the object in that file. In this
    case it is ``simulation``.

Finally, we give the number of ``subcycles``, which indicates how many MD
steps ``OpenMM`` performs before |pyretis| requests another frame. This value
should be relatively high when using GPUs, as it determines how often |pyretis|
communicates with the GPU. However, keep in mind that this also lowers the
time resolution of your |pyretis| simulation.

After this, the setup can be run as a regular |pyretis| simulation using
``pyretisrun``:

.. code-block:: pyretis

   pyretisrun -i openmm_retis.rst -p

.. _examples-permeability:

Studying permeability with |pyretis|
====================================

This example shows how to set up a permeability simulation with |pyretis|.
For further details on the derivation of the formulas and description of the
Monte Carlo moves, please read and cite the permeability from (RE)TIS paper [1]_.

This example uses three non-interacting particles on a
:download:`flat potential </_static/examples/permeability/flat_potential.py>`.
It walks through this
:download:`retis_perm.rst </_static/examples/permeability/retis_perm.rst>`
input file, together with this
:download:`initial.xyz </_static/examples/permeability/initial.xyz>`.

New simulation settings
-----------------------
The ``Simulation`` section has a few extra options:

.. literalinclude:: /_static/examples/permeability/retis_perm.rst
   :language: rst
   :lines: 9-12

Here we have:

  * ``zero_left``: Tells |pyretis| that the ``[0^-]`` ensemble has a left
    boundary that is not located at ``-inf``.

  * ``permeability``: If ``True``, any path in the ``[0^-]`` ensemble that
    starts and ends at one of the interfaces is accepted. If ``False``, any
    path that hits the ``zero_left`` interface will be rejected, leading to
    incorrect flux calculations. This option also triggers ``pyretisanalyse``
    to calculate
    :math:`\xi`, :math:`\frac{\tau}{dz}` and the permeability.

  * ``swap_attributes``: A list of ``path_ensemble`` attributes that need to
    be swapped with every REPEX swap. The newly implemented Monte Carlo moves
    alter the order parameter in ``[0^-]``. These altered order parameters
    need to be swapped with the path whenever the path is exchanged between
    ensembles.


New TIS settings
----------------
The new ``mirror`` and ``target-swap`` moves add a few options to the
``TIS`` section:

.. literalinclude:: /_static/examples/permeability/retis_perm.rst
   :language: rst
   :lines: 43-45

Here we have:

  * ``mirror_freq``: The probability of attempting the ``mirror`` move in the
    ``[0^-]`` ensemble.

  * ``target_freq``: The probability of attempting the ``target-swap`` move in
    the ``[0^-]`` ensemble.

  * ``target_indices``: A list of atom indices. The target swap is only
    attempted between these atoms. Make sure that the original
    ``orderparameter.index`` is included in this list.

New order parameter class
-------------------------
This simulation can be run using the new order parameter class
``Permeability``. This class is a subclass of ``Position``, but alters the
output depending on
``mirror_pos``, ``relative`` and ``offset``.

.. literalinclude:: /_static/examples/permeability/retis_perm.rst
   :language: rst
   :lines: 73-80

Here we have:

  * ``dim``: The same as for the class ``Position``.

  * ``index``: The index of the particle that will be tracked at the
    start of the simulation. This attribute will be changed by the ``target-swap`` move.

  * ``offset``: This order parameter adds an ``offset`` to the value of
    ``compute_s()`` before wrapping it into the periodic box. This alters
    the resulting OP value, but all values will fall within the box vectors.
    If you want to alter the box vectors instead and do not have access
    to them, you can use ``PermeabilityMinusOffset``, which subtracts the offset
    after wrapping, before returning the value.

  * ``relative``: If ``True``, the output is mapped relative to the box vector
    (between 0 and 1). Both ``offset`` and ``mirror_pos`` should be defined
    relative to this box vector as well.

  * ``mirror_pos``: The position of the mirror plane, before applying the
    offset. For the current implementation, this **must** be set halfway
    between the ``0-R`` and ``0-L`` interfaces.

The ``Permeability`` classes call the function ``compute_s()`` before applying
the offset and mirror. For the base class this calls the ``compute`` function of
``Position``. To use this with your own OP, you can make a subclass
of ``Permeability`` and override the ``self.compute_s()`` function to return
your own custom OP before applying the offset and mirroring.

Output of the new moves
-----------------------

The new moves also lead to some new possible responses in pathensemble.txt.
For the ``mirror`` move, which is always accepted subject to the constraint on
``mirror_pos``, this is just a new move type called ``mr``.

For the ``target_swap`` move:

  * A new generated label ``ts`` to indicate target-swap.

  * A new rejection reason: ``TSS``, which means there are no valid indices to
    swap to.

  * Another new rejection reason ``TSA``, which is a rejection based on the
    Monte Carlo acceptance.

Another changed response is ``BTS`` (backward too short), which is a more
common rejection for the ``[0^-]`` <-> ``[0^+]`` swap. This indicates that
the attempted trajectory in ``[0^-]`` ended at the ``L`` interface, so we do
not attempt to extend it into the ``[0^+]`` ensemble.

New analysis options
--------------------
Adding ``permeability = True`` to the ``Simulation`` settings triggers
``pyretisanalyse`` to also calculate and plot :math:`\xi`,
:math:`\frac{\tau}{dz}` and the permeability. This follows the formulas
described in the paper [1]_.
For the calculation of :math:`\tau` to work, a reference region has to be
chosen. This is done by adding ``tau_ref_bin`` to the ``Analysis`` section of
the ``.rst`` file.

.. literalinclude:: /_static/examples/permeability/retis_perm.rst
   :language: rst
   :lines: 148-150


This value can be changed in ``out.rst`` and the analysis can then be rerun
with another reference region. The analysis also plots a 10-bin histogram of
the ``[0^-]`` region to help the user select a flat histogram region in this
space.



References
----------

.. [1] Permeability from (RE)TIS, A. Ghysels et al. (manuscript in preparation)

.. _user-guide-intro:

Introduction to |pyretis| and rare event methods
================================================

|pyretis| is a computational library for performing molecular
simulations of rare events with a focus on transition
interface sampling (TIS) [1]_
and replica exchange transition interface sampling (RETIS) [2]_.

Rare events are called rare because they happen at time or length scales
much longer than we are able to simulate through brute-force simulations.
Raindrops form every day, but compared to the motion of water molecules,
raindrop formation is rare, and simulating it through brute-force
simulations would take forever and a day.

Rare event methods aim to sample these rare events inaccessible
through brute-force simulation.


.. _user-guide-tis-theory:

Transition Interface Sampling
-----------------------------

To explain the basic concepts, we
consider an illustrative example: The transition
between two stable states, from the reactant (labeled A) to the product
(labeled B) as illustrated in :numref:`fig_intro_pot`.

.. _fig_intro_pot:

.. figure:: /_static/img/intro/potential.png
    :class: img-responsive center-block
    :alt: Potential example
    :width: 50%
    :align: center

    A potential energy barrier separating two stable states A and B.
    The order parameter measures the extent of the reaction in this
    particular energy landscape and a possible trajectory for a
    particle making the transition is illustrated as a black arrow.


These two states are defined by a progress
coordinate (or order parameter), :math:`\lambda`,
where state A is for :math:`\lambda \leq \lambda_\text{A}`
and state B for :math:`\lambda \geq \lambda_\text{B}`.
The central quantity calculated in TIS/RETIS simulations is the
rate constant :math:`k_\text{AB}` for the transition
:math:`\text{A} \to \text{B}` which can be expressed as:

.. math::

   k_{\text{AB}} = f_\text{A} \mathcal{P}_{\text{A}} (\lambda_\text{B} \vert \lambda_\text{A}),

and we see that it contains two parts:

* The initial flux :math:`f_\text{A}` which measures how often
  trajectories start off at the foot of the reaction barrier
  from the reaction side :math:`\lambda_\text{A}`. In TIS this
  is obtained from a molecular dynamics (MD) simulation which
  in |pyretis| is requested by running a
  :ref:`md-flux <user-section-simulation-mdflux>` simulation. Note
  that this is not needed for the :ref:`RETIS method <user-guide-retis-theory>` as explained.

* The crossing probability
  :math:`\mathcal{P}_{\text{A}} (\lambda_\text{B} \vert \lambda_\text{A})`
  of reaching :math:`\lambda_\text{B}` before :math:`\lambda_\text{A}`
  given that :math:`\lambda_\text{A}` has just been crossed.
  In |pyretis| this can be accomplished by running a
  :ref:`tis <user-section-simulation-tis>` or
  :ref:`retis <user-section-simulation-retis>` simulation.


As this is a rare event, the crossing probability is
extremely small and nearly impossible to compute in brute-force simulations.
Transition interface sampling divides the region between A and B into
sub-regions using interfaces, denoted :math:`\lambda_i` (see :numref:`fig_intro_tis`).
The first interface, :math:`\lambda_0` is positioned at the interface defining
state A (:math:`\lambda_A`), and
the next interface, :math:`\lambda_1` is put
at :math:`\lambda_1 > \lambda_0` such that the probability of reaching
:math:`\lambda_1` from :math:`\lambda_0` is no longer extremely small. We can
continue in this fashion and place :math:`N` more interfaces until we
reach :math:`\lambda_{N} = \lambda_\text{B}`. What we effectively have
done with these :math:`N+1` interfaces is to split up the
computation of the extremely small crossing probability into the
computation of many, not-so small, crossing probabilities:

.. math::

   k_{\text{AB}} = f_{\text{A}}  \prod_{i=0}^{N} P_{\text{A}} (\lambda_{i+1}|\lambda_{i}).


Here, the intermediate crossing probabilities,
:math:`P_{\text{A}} (\lambda_{i+1}|\lambda_{i})`,
are  formally defined as the
probability of a path crossing
:math:`\lambda_{i+1}` given
that it originated from :math:`\lambda_\text{A}`,
ended in :math:`\lambda_\text{A}` or :math:`\lambda_\text{B}`,  and had
at least one crossing with :math:`\lambda_i` in the past.

The interfaces we have placed define the so-called **path ensembles**.
A path ensemble comprises all possible
trajectories that start at the foot of reaction
barrier from the reactant side (:math:`\lambda_\text{A}`),
end in it or at the product region (:math:`\lambda_\text{B}`) and
having reached a certain threshold value (:math:`\lambda_i`) between the start point
and the final point. This path ensemble is labelled as :math:`\left[i^{+} \right]`.
The probabilities,
:math:`\mathcal{P}_{\text{A}} (\lambda_{i+1} \vert \lambda_{i})` can then
be obtained as the fraction of paths in the :math:`\left[i^+\right]` ensemble
that also cross :math:`\lambda_{i+1}`.


.. _fig_intro_tis:

.. figure:: /_static/img/intro/tis.png
    :class: img-responsive center-block
    :alt: TIS interfaces
    :width: 50%
    :align: center

    Illustration of TIS interfaces placed along the order parameter
    in a system where a potential energy barrier separates two stable
    states A and B. The interfaces define different ensembles and here,
    two trajectories are shown. One (black) is reactive, reaching the final state,
    while the other (orange) just reaches the intermediate :math:`\lambda_2`
    interface.


What is left now, is to have an efficient way of generating trajectories
for the various path ensembles.
This is in fact done by making use of a selection of Monte Carlo (MC) moves.
For TIS [1]_ we choose between two moves:

* The **shooting move** which is adapted from the
  transition path sampling (TPS) shooting algorithm [3]_ [4]_
  to allow variable trajectory length. In this move, we generate
  a new trajectory from an existing trajectory one by:

  1. Picking randomly one of the discrete MD steps in the present
     trajectory.

  2. Modifying the velocities of this phase point (e.g. by randomly
     drawing new velocities from a Maxwellian distribution).

  3. Generating a new trajectory from this new phase point by integrating
     (i.e. running MD simulations) forward and backward in time until A or
     B is reached. The new trajectory is then obtained by merging the backward
     and forward trajectories.

  The new trajectory is accepted as part of :math:`[i^+]` only if all the following criteria are satisfied:

  1. A detailed balance condition for the energy and path length. [1]_

  2. It starts at :math:`\lambda_A`

  3. It has at least one crossing with :math:`\lambda_i` before ending in A or B.

  The shooting move gives a much higher
  chance to generate a valid trajectory at each
  trial compared to simply starting from a random phase point within the reactant well.

* The **time reversal move** which generates trajectories by simply changing the
  time direction of a path. [1]_

These two moves are illustrated in :numref:`fig_intro_retis_moves`.
See :ref:`user-guide-path-move-types` for a complete list of all move
type codes used in |pyretis| output and input files.



.. _user-guide-retis-theory:

Replica Exchange Transition Interface Sampling
----------------------------------------------

The RETIS method is similar to TIS and employs both the shooting
and time reversal moves.
In addition, RETIS makes use of the **swapping move** and defines a
new ensemble :math:`[0^-]`
which consist of trajectories that explore the reactant state
(see the illustration in :numref:`fig_intro_retis`).

.. _fig_intro_retis:

.. figure:: /_static/img/intro/retis.png
    :class: img-responsive center-block
    :alt: RETIS interfaces
    :width: 50%
    :align: center

    Illustration of RETIS interfaces placed along the order parameter
    in a system where a potential energy barrier separates two stable
    states A and B. The interfaces define different ensembles and here,
    two trajectories are shown. One (black) is reactive, reaching the final state,
    while the other (orange) just reaches the intermediate :math:`\lambda_2`
    interface. In RETIS a special path ensemble, :math:`[0^-]` is also considered
    as described in the text.


The swapping move acts between different path simulations. If two
simulations generate simultaneously two paths that are valid for each
other's path ensemble, these two paths can be swapped.
The swapping moves increase with negligible extra computational cost
the number of accepted paths in the ensembles and decrease significantly the
correlations between the consecutive paths within the same ensemble. All the
moves used for generating trajectories are illustrated in :numref:`fig_intro_retis_moves`.


.. _fig_intro_retis_moves:

.. figure:: /_static/img/intro/retis_moves.png
    :class: img-responsive center-block
    :alt: RETIS moves
    :width: 50%
    :align: center

    Illustration of the RETIS moves for generating trajectories.
    A contour plot of a hypothetical free energy
    surface along a progress coordinate and an arbitrary second coordinate
    is shown and 4 interfaces (:math:`\lambda_0`, :math:`\lambda_1`, :math:`\lambda_2`
    and :math:`\lambda_3`) have been positioned along the progress coordinate.
    Three different RETIS moves (shooting, time reversal and swapping) are shown for
    the :math:`[i^+] = [2^+]` path ensemble. The old paths are in blue and the new paths
    after (a successful) completion of the MC moves are shown in red. The orange line
    show the interface that needs to be crossed for a valid path in the current ensemble.

There is one notable exception where the swapping move is more computationally demanding: the
swap between the :math:`[0^+]` and :math:`[0^-]` ensembles requires MD simulations. [2]_

In a TIS simulation, as explained above, we have to perform an extra simulation in order
to calculate the initial flux. In RETIS, this initial flux can directly be obtained
by:

.. math::

   f_{\text{A}} = \frac{1}{\left \langle t_{\rm path}^{[0^+]} \right \rangle + \left \langle t_{\rm path}^{[0^-]}\right \rangle }

where :math:`\left \langle t_{\rm path}^{[0^+]} \right \rangle` is the average path length
in the :math:`[0^+]` ensemble and
:math:`\left \langle t_{\rm path}^{[0^-]}\right \rangle` the average path length
in the :math:`[0^-]` ensemble.

In |pyretis|, RETIS simulations are requested by setting the
:ref:`simulation task <user-section-simulation>` to :ref:`RETIS <user-section-simulation-retis>`.


References
~~~~~~~~~~

.. [1] T. S. van Erp, D. Moroni and P. G. Bolhuis,
       J. Chem. Phys. 118, 7762 (2003)
       https://dx.doi.org/10.1063%2F1.1562614

.. [2] T. S. van Erp,
       Phys. Rev. Lett. 98, 26830 (2007)
       https://dx.doi.org/10.1103/PhysRevLett.98.268301

.. [3] C. Dellago, P. G. Bolhuis, F. S. Csajka, and D. Chandler,
       J. Chem. Phys. 108, 1964 (1998)
       https://dx.doi.org/10.1063/1.475562

.. [4] P. G. Bolhuis, D. Chandler, C. Dellago and P. L. Geissler,
       Annu. Rev. Phys. Chem. 53, 291 (2002)
       https://dx.doi.org/10.1146/annurev.physchem.53.082301.113146

.. _examples-pyvisa-analysis:

Post-processing and visualization with |pyvisa|
===============================================

In this example, we perform post-processing on the methane
hydrate system from one of the previous examples called "Using GROMACS".
All post-processing is done with the |pyvisa| analysis tab, shown in
:numref:`fig_analysis_tab`, and the recalculation tool.
We will continue the study of cage-to-cage diffusion within a S1 hydrate
and add new collective variables.
Then we perform clustering, principal component analysis (PCA), and
visualization of the new results. These analysis methods require the
``scikit-learn`` and ``scipy`` packages.

.. _fig_analysis_tab:

.. figure:: /_static/examples/pyvisa-analysis/fig_analysis_tab.png
    :alt: PyVisA analysis-tab
    :width: 50%
    :align: center

    Analysis tab of |pyvisa| with options for interactivity, animation,
    and post-processing.


Before we proceed with the analysis, we need to have some simulation
results, so please run the example :ref:`Using GROMACS <examples-gromacs-hydrate>` first.

The example is structured as follows:

.. contents::
   :local:


Visualization and compression
-----------------------------

When the simulation has finished, the results can be visualized and
post-processed. It is always important to save the original data, because
the recalculation creates new |pyvisa| data from the trajectory files of
the old simulation. This loses the frames between stored trajectory files.
The following steps are suggested:

1. Compress the original data into an HDF5 file using:

   .. code-block:: pyretis

      pyvisa -i out.rst -cmp

2. Compress only the order parameter data into an ``.hdf5.zip`` file using:

   .. code-block:: pyretis

      pyvisa -i out.rst -cmp -oo

3. Rename these files to ``standard_simulation.hdf5.zip``.

The results can be visualized from the ``.rst`` files or from the
compressed files.
If you want to visualize the data, the following commands can be run:

.. code-block:: pyretis

  pyvisa -i <input-file>

.. code-block:: pyretis

  pyvisa -i <rst-file> -data <data>

where ``<input-file>`` can be either an ``.rst`` file or an
``.hdf5.zip`` file. If you are using the ``-data`` command, then
``<data>`` can be ``all`` to load all data, a trajectory file, or a
simulation directory. For example, the first path ensemble can be loaded
by using ``000`` as ``<data>``.

With |pyvisa|, you can visualize trajectories from all path ensembles and
sort them by status and Monte Carlo move. Try finding and visualizing a
reactive pathway by clicking on the plot and pressing "Show trj" in the
Analysis tab. This highlights the points that belong to the selected trajectory.
You can also customize how the visualization looks. In :numref:`fig_dens_reactive`
a reactive pathway has been shown.

.. _fig_dens_reactive:

.. figure:: /_static/examples/pyvisa-analysis/fig_dens_reactive.png
    :alt: A reactive pathway
    :width: 70%
    :align: center

    Density plot of the order parameter and the kinetic energy from the methane
    hydrate example with a reactive pathway shown in green.


Recalculation of new collective variables
-----------------------------------------

Now that the old data has been stored, we can add two more collective
variables to the simulation using |pyvisa|.
The collective variables we introduce are the area of the six-membered ring
that the methane jumps through, and the volume of the starting cage. These
descriptors will try to capture the breathing of the starting cage.
The following lines must be added to the retis.rst file:

.. literalinclude:: /_static/examples/pyvisa-analysis/order_rec.rst
   :language: rst

With the collective variable added to ``retis.rst``, we also need to add
the Python script for the calculation in the ``orderp.py`` file. The
``scipy`` package is needed to run this script, so make sure that it is
installed.
The full script is given here:

.. pyretis-collapse-block::
   :heading: Show/hide the Python script for the new collective variable.

   .. literalinclude:: /_static/examples/pyvisa-analysis/orderp.py
      :language: python
      :lines: 20-59

As we also have the z-coordinate of the methane molecule as a collective
variable, we add this to ``orderp.py`` so that the recalculation tool can
use it. The ``Position`` descriptor is shown here:

.. pyretis-collapse-block::
   :heading: Show/hide the Python script for the Position collective variable.

   .. literalinclude:: /_static/examples/pyvisa-analysis/orderp.py
      :language: python
      :lines: 62-125


To begin the recalculation, start |pyvisa| by loading all data with the
modified ``retis.rst`` file:

.. code-block:: pyretis

  pyvisa -i retis.rst

It is important to load |pyvisa| with the modified ``retis.rst`` file so
the program knows what has been added. The recalculation can then be
started from the |pyvisa| file menu by selecting "Recalculate data". The
warning shown in :numref:`fig_recalc_note` will appear. Press yes if you
have added the post-processing requirements to both ``retis.rst`` and
``orderp.py``.

.. _fig_recalc_note:

.. figure:: /_static/examples/pyvisa-analysis/recalc_note.png
    :alt: Warning issued before recalculation of new collective variables.
    :width: 50%
    :align: center

    The warning issued by |pyvisa| before the recalculations.

The recalculation will now start. When the procedure is finished, |pyvisa|
loads the new data into the GUI and displays a message that the new data
can be visualized.

Post-processing and clustering
------------------------------
Now that we have the new data, we can use |pyvisa| to perform clustering
and PCA. The following steps can be done in the data exploration:

1. Go to the Analysis tab of |pyvisa| and press the button "Show correlations".
   This will produce the correlation matrix. From these results, plot
   the order parameter and one of the collective variables.

2. Select a number of components to use in clustering and from the
   Analysis-selection, pick an algorithm for clustering and press the
   Analysis-button. This will produce a cluster plot of the
   chosen variables.
   Try to start with k-means, and then try Gaussian mixture and spectral
   clustering to see if there is a difference between the methods.

3. Try to perform a principal component analysis of the results. Begin by
   selecting 3 components, and PCA as the method, and press the analysis-
   button. This will produce the loading matrix, the scores plot from the
   first two components, and the cumulative explained variance. How much
   variance was retained? Were three components enough? Are there any
   strong correlations between the principal components and the original
   descriptors?

If you want to study the results from the principal component analysis
further, the data will be stored in an HDF5 file in your simulation
directory containing all the simulation data, and the data from the
analysis.

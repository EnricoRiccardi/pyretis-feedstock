# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This package contains functions for setting up plotters.

Specifically, it defines colors, color schemes and a function
for selecting a plotter.

Here we also hard-code some color schemes which may be useful for
the plotting:

* the colorblind 10 scheme [cb10]_,

* the deep color scheme from the seaborn project [deep]_ and

* the husl color scheme [husl]_.

Package structure
-----------------

Modules
~~~~~~~

__init__.py
    This file. Handles imports for PyRETIS and defines some colors
    and the method for creating a plotter.

mpl_plotting.py (:py:mod:`pyretis.inout.plotting.mpl_plotting`)
    Methods for generating plots using matplotlib.

plotting.py (:py:mod:`pyretis.inout.plotting.plotting`)
    Definition of a generic base class for the plotter(s).

txt_plotting.py (:py:mod:`pyretis.inout.plotting.txt_plotting`)
    Methods for generating text output.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

create_plotter (:py:func:`.create_plotter`)
    Method to create a plotter.

Folders
~~~~~~~

styles
    This folder contains style files for matplotlib which are used by
    PyRETIS.

References
~~~~~~~~~~

.. [cb10] The colorblind_10 color scheme,
          https://jiffyclub.github.io/palettable/tableau/
.. [deep] The deep color scheme, from the seaborn project
          http://stanford.edu/~mwaskom/software/seaborn/index.html
.. [husl] The husl color scheme, http://www.husl-colors.org/

"""
from pyretis.inout.settings import SECTIONS
from .mpl_plotting import MplPlotter, mpl_set_style
from .txt_plotting import TxtPlotter


__all__ = ['create_plotter']


# Custom named colors:
COLORS = {'almost_black': '#262626'}
# Custom color-schemes. The default will be defined by the style file.
# The husl schemes are suited when many different colors are needed. They
# are hard-coded with the different number of colors.
COLOR_SCHEME = {'colorblind_10': ['#006BA4', '#FF800E', '#ABABAB', '#595959',
                                  '#5F9ED1', '#C85200', '#898989', '#A2C8EC',
                                  '#FFBC79', '#CFCFCF'],
                'deep': ['#4C72B0', '#55A868', '#C44E52', '#8172B2',
                         '#CCB974', '#64B5CD'],
                'husl_10': ['#f67088', '#db8831', '#ad9c31', '#77aa31',
                            '#33b07a', '#35aca4', '#38a8c5', '#6e9af4',
                            '#cc79f4', '#f565cc'],
                'husl_15': ['#f67088', '#f37932', '#ca9131', '#ad9c31',
                            '#8ea531', '#4fb031', '#33b07a', '#34ad99',
                            '#36abae', '#38a8c5', '#3ba3ec', '#9491f4',
                            '#cc79f4', '#f45fe3', '#f569b7'],
                'husl_20': ['#f67088', '#f7754f', '#db8831', '#c29431',
                            '#ad9c31', '#96a331', '#77aa31', '#31b23e',
                            '#33b07a', '#34ae92', '#35aca4', '#36abb3',
                            '#38a8c5', '#3aa5de', '#6e9af4', '#a38cf4',
                            '#cc79f4', '#f45bf1', '#f565cc', '#f66bad']}


def create_plotter(plot_settings, out_dir=None):
    """Create a plotter from given settings.

    The input plot settings is assumed to be a dictionary which we use
    for creating the plotter. In case the plot settings are not given, we
    just return None. We are here assuming that they are not given
    simply because we do not want to create a plotter.

    Parameters
    ----------
    plot_settings : dict
        These are the settings to create the plotter from. Here, we
        look for the keys `plotter`, `output` and `style` which defines
        the plotter to use, the output format and the style to use.
    out_dir : string, optional
        This string selects if and where the output should be written
        to. It's mainly used for internal purposes to save the report
        files to a specific directory.

    Returns
    -------
    out : object like :py:class:`.MplPlotter`
        This is an object which can be used for plotting.

    """
    default = SECTIONS['analysis']['plot']
    plotter = plot_settings.get('plotter', default['plotter'])
    out_fmt = plot_settings.get('output', default['output'])
    style = plot_settings.get('style', default['style'])
    backup = plot_settings.get('backup', SECTIONS['output']['backup'])
    if plotter.lower() in ['mpl', 'matplotlib']:
        return MplPlotter(out_fmt, backup=backup, style=style,
                          out_dir=out_dir)
    msg = f'Unknown plotter: {plotter}'
    raise ValueError(msg)

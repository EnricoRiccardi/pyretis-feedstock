# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Functions for generating plots using matplotlib.

This module defines a plotter class for matplotlib and it also defines
some standard plots that are used in the analysis.

Important classes defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MplPlotter (:py:class:`.MplPlotter`)
    A class for plotting with matplotlib.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mpl_set_style (:py:func:`.mpl_set_style`)
    Method for setting the style for the plots, typically used here to
    load the *PyRETIS style*.
"""
import os
import logging
import numpy as np
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.collections import LineCollection
import matplotlib.style
from pyretis.inout.plotting.plotting import Plotter
from pyretis.inout.common import create_backup, name_file
from pyretis.inout.common import (ENERFILES, ENERTITLE, FLUXFILES,
                                  ORDERFILES, PATHFILES, PATH_MATCH)

np.set_printoptions(legacy='1.21')

logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


__all__ = ['MplPlotter']


# Define default style file:
_STYLEFILE = 'pyretis.mplstyle'
_MPL_STYLE_FILE = os.sep.join([os.path.dirname(__file__), 'styles',
                               _STYLEFILE])
_TITLE_SETTINGS = {'loc': 'right'}


class MplPlotter(Plotter):
    """A plotter using matplotlib.

    This class defines a plotter which makes use of matplotlib and
    it can be used as a starting point to create other plotters
    for PyRETIS bases on other tools (e.g. gnuplot etc.).

    Attributes
    ----------
    style : string
        Defines what style to use for the plotting.
    out_fmt : string
        Selects format for output plots.

    """

    def __init__(self, out_fmt, backup=False, style=None, out_dir=None):
        """Initiate, set style and check output format.

        Parameters
        ----------
        out_fmt : string
            This is the format to use for the images.
        backup : boolean, optional
            Determines if we should overwrite or backup old files.
        style : string, optional
            This selects the style to use, it can be a file path or the
            string with the style name.
        out_dir : string, optional
            Determines if we should write the files to a particular
            directory.

        """
        super().__init__(backup=backup, plotter_type='matplotlib',
                         out_dir=out_dir)
        self.style = style
        mpl_set_style(self.style)
        # Check if the requested file format is something we can do:
        fig = Figure()
        canvas = FigureCanvas(fig)
        supported = canvas.get_supported_filetypes()
        if out_fmt not in supported:
            msg = [f'Output format "{out_fmt}" is not supported.',
                   'Please try:']
            for key in supported:
                msg.append(key)
            raise ValueError(' '.join(msg))
        self.out_fmt = out_fmt
        del fig
        del canvas

    def _print_figures_to_file(self, canvas):
        """Save figure(s) to file(s).

        This function will save figures to files. It will add a file
        format extension to the files when writing.

        Parameters
        ----------
        canvas : dict of :class:`matplotlib.backend_bases.FigureCanvasBase`
            ``canvas[key]`` is assumed to define a figure which we will
            save to a file with a file name defined by ``key`` and the
            extension defined by ``self.out_fmt``.

        Returns
        -------
        out : dict
            The files written.

        """
        outputfiles = {}
        for key, val in canvas.items():
            local_file = name_file(key, self.out_fmt, path=None)
            full_path = name_file(key, self.out_fmt, self.out_dir)
            outputfiles[key] = local_file
            mpl_savefig(val, full_path, self.backup)
        return outputfiles

    def output_flux(self, results):
        """Plot flux results using py:func:`.mpl_plot_flux`.

        The parameters for this method is described in
        :py:func:`.mpl_plot_flux`.

        Returns
        -------
        outputfile : list of dicts
            A list containing the files created for the flux and for
            the error in the flux.

        Note
        ----
        The returned list is used to plot the figures in *pairs*.

        """
        canvas_run, canvas_err = mpl_plot_flux(results)
        # Restructure output files for reporting
        outputfiles = []
        for run, err in zip(canvas_run, canvas_err):
            local_run = name_file(run['name'], self.out_fmt, path=None)
            local_err = name_file(err['name'], self.out_fmt, path=None)
            full_run = name_file(run['name'], self.out_fmt, path=self.out_dir)
            full_err = name_file(err['name'], self.out_fmt, path=self.out_dir)
            mpl_savefig(run['canvas'], full_run, self.backup)
            mpl_savefig(err['canvas'], full_err, self.backup)
            outputfiles.append({'runflux': local_run,
                                'errflux': local_err})
        return outputfiles

    def output_energy(self, results, energies):
        """Plot energy results using :py:func:`.mpl_plot_energy`.

        The parameters for this method is described in
        :py:func:`.mpl_plot_energy`.

        Returns
        -------
        out : dict
            This dict contains the files created by the plotting.

        """
        canvas = mpl_plot_energy(results, energies)
        return self._print_figures_to_file(canvas)

    def output_orderp(self, results, orderdata):
        """Plot order parameter data using :py:func:`.mpl_plot_orderp`.

        The parameters for this method is described in
        :py:func:`.mpl_plot_orderp`.

        Returns
        -------
        out : dict
            This dict contains the files created by the plotting.

        """
        canvas = mpl_plot_orderp(results, orderdata)
        return self._print_figures_to_file(canvas)

    def output_path(self, results, path_ensemble):
        """Plot path results using :py:func:`.mpl_plot_path`.

        The parameters for this method is described in
        :py:func:`.mpl_plot_path`.

        Returns
        -------
        out : dict
            This dict contains the files created by the plotting.

        """
        canvas = mpl_plot_path(results, path_ensemble)
        return self._print_figures_to_file(canvas)

    def output_pppath(self, results, path_ensemble):
        """Plot repptis path results using :py:func:`.mpl_plot_pppath`.

        The parameters for this method is described in
        :py:func:`.mpl_plot_pppath`.

        Returns
        -------
        out : dict
            This dict contains the files created by the plotting.

        """
        canvas = mpl_plot_pppath(results, path_ensemble)
        return self._print_figures_to_file(canvas)

    def output_pp_global_cross(self, pp):
        """Plot repptis global cross prob using :py:func:`.mpl_plot_path`.

        The parameters for this method is described in
        :py:func:`.mpl_plot_path`.

        Returns
        -------
        out : dict
            This dict contains the files created by the plotting.

        """
        series = [{'type': 'xy',
                   'x': list(range(len(pp))),
                   'y': pp}]
        figset = {'xlabel': r'$\lambda_i$',
                  'ylabel': r'$P_A^{(\lambda_i|\lambda_A)}$',
                  'title': 'Global crossing probability', 'yscale': 'log'}
        canvas = mpl_simple_plot(series, fig_settings=figset)
        return self._print_figures_to_file({'pp_Pc': canvas})

    def output_matched_probability(self, path_ensembles, detect,
                                   matched, reptis=False):
        """Plot matched probabilities using :py:func:`.mpl_plot_matched`.

        The parameters for this method is described in
        :py:func:`.mpl_plot_matched`.

        Returns
        -------
        out : dict
            This dict contains the files created by the plotting.

        """
        canvas = mpl_plot_matched(path_ensembles, detect, matched, reptis)
        return self._print_figures_to_file(canvas)

    def output_xi(self, xi, xierror):
        """Plot Xi and Xi error using :py:func:`.mpl_plot_xi`.

        The parameters for this method is described in
        :py:func:`.mpl_plot_xi`.

        Returns
        -------
        out : dict
            This dict contains the files created by the plotting.

        """
        canvas_run, canvas_err = mpl_plot_xi(xi, xierror)
        local_run = name_file(canvas_run[0]['name'], self.out_fmt, path=None)
        full_run_path = name_file(canvas_run[0]['name'],
                                  self.out_fmt,
                                  self.out_dir)
        mpl_savefig(canvas_run[0]['canvas'], full_run_path, self.backup)

        local_err = name_file(canvas_err[0]['name'], self.out_fmt, path=None)
        full_err_path = name_file(canvas_err[0]['name'],
                                  self.out_fmt,
                                  self.out_dir)
        mpl_savefig(canvas_err[0]['canvas'], full_err_path, self.backup)

        return {"xirun": local_run, "xierror": local_err}

    def output_tau(self, tau, tauerror, bins, tau_ref, tau_ref_bins):
        """Plot Tau, Tau error, and tau reference histogram.

        The parameters for this method is described in
        :py:func:`.mpl_plot_tau` and :py:func:`.mpl_plot_tau_ref`.

        Returns
        -------
        out : dict
            This dict contains the files created by the plotting.

        """
        canvas_run, canvas_err = mpl_plot_tau(tau, tauerror)
        canvas_ref = mpl_plot_tau_ref(tau_ref, tau_ref_bins, bins)
        local_run = name_file(canvas_run[0]['name'], self.out_fmt, path=None)
        full_run_path = name_file(canvas_run[0]['name'],
                                  self.out_fmt,
                                  self.out_dir)
        mpl_savefig(canvas_run[0]['canvas'], full_run_path, self.backup)

        local_err = name_file(canvas_err[0]['name'], self.out_fmt, path=None)
        full_err_path = name_file(canvas_err[0]['name'],
                                  self.out_fmt,
                                  self.out_dir)
        mpl_savefig(canvas_err[0]['canvas'], full_err_path, self.backup)
        local_ref = name_file(canvas_ref['name'], self.out_fmt, path=None)
        full_ref_path = name_file(canvas_ref['name'],
                                  self.out_fmt,
                                  self.out_dir)
        mpl_savefig(canvas_ref['canvas'], full_ref_path, self.backup)
        return {"taurun": local_run, "tauerror": local_err,
                "tauref": local_ref}


def _mpl_read_style_file(filename):
    """Read style files for old versions of matplotlib.

    This function is just intended for old versions of matplotlib
    where we have to read parameters and set them ourselves.

    Parameters
    ----------
    filename : string
        This is the matplotlib rc file to open.

    Returns
    -------
    out : None
        Returns None, but modifies `matplotlib.rcParams`.

    """
    with open(filename, 'r', encoding="utf8") as fileh:
        for lines in fileh:
            linesc = lines.strip().split('#')[0]
            loc = linesc.find(':')
            if loc != -1:
                key = linesc[:loc].strip()
                value = linesc[loc+1:].strip()
                if key.find('color') != -1:
                    value = f'#{value}'
                try:
                    matplotlib.rcParams[key] = value
                except KeyError:
                    logger.warning(('Unknown setting "%s". '
                                    'Please update matplotlib'), key)


def mpl_set_style(style='pyretis'):
    """Set the plotting style for matplotlib.

    This will set up the plotting according to some given style.
    Styles can be given as string, for instance 'ggplot', 'bmh',
    'grayscale' (i.e. one of the styles in `matplotlib.style.available`)
    or as a file (full path is needed). The default PyRETIS style
    is stored in `_MPL_STYLE_FILE` and can be selected with 'pyretis'.
    Style equal to None is just the default matplotlib style.

    Parameters
    ----------
    style : string, optional
        This selects the style to use, it can be a file path or the
        string with the style name.

    """
    if style == 'pyretis':
        style = _MPL_STYLE_FILE
    if style in matplotlib.style.available:
        logger.info('Loading matplotlib style: %s', style)
        matplotlib.style.use(style)
    else:  # assume this is just a file
        logger.info('Loading matplotlib style from file: %s', style)
        rcpar = matplotlib.rc_params_from_file(
            style,
            use_default_template=False
        )
        matplotlib.rcParams.update(rcpar)


def mpl_savefig(canvas, outputfile, backup=False):
    """Write/save matplotlib figures to files.

    Parameters
    ----------
    canvas : object like :class:`matplotlib.backend_bases.FigureCanvasBase`
        This is the figure to be written to the file by
        using ``canvas.print_figure()``.
    outputfile : string
        This is the name of the output file to create.
    backup : boolean, optional
        This determines if we should try to back-up old versions of the
        figures.

    """
    if backup:
        msg = create_backup(outputfile)
        if msg:
            logger.warning(msg)
    canvas.print_figure(outputfile)


def mpl_plot_in_chunks(axs, series, chunksize=20000):
    """Plot a series in chunks using matplotlib.

    When plotting 'large' datasets, matplotlib might give an
    ``OverflowError: Allocated too many blocks`` error.
    Here we avoid this error by plotting the data in chunks. We could
    also downsample the data, but this is perhaps something best left
    to the user.

    Parameters
    ----------
    axs : Axes object from matplotlib
        Where to do the plotting
    series : dict
        Represents the data to be plotted.
    chunksize : int
        This is the maximum size we will try to plot in one go.

    """
    color = None
    line = None
    leny = min(len(series['y']), len(series['x']))
    if leny > chunksize:
        nchunk, rest = divmod(leny, chunksize)
        for i in range(nchunk):
            low = i * chunksize
            high = low + chunksize
            line = _mpl_plot_xy_chunk(axs, series, low=low, high=high,
                                      color=color)
            color = line.get_color()
        if rest > 0:
            # We want 1 less to force line overlap,
            # We want positive low in case x and y length don't overlap
            line = _mpl_plot_xy_chunk(axs, series, low=nchunk*chunksize-1,
                                      high=None, color=color)
    else:
        line = _mpl_plot_xy_chunk(axs, series)
    return line


def _mpl_plot_xy_chunk(axs, series, low=0, high=None, color=None):
    """Do the actual plotting in chunks for ``x`` vs ``y`` data.

    Parameters
    ----------
    axs : Axes object from matplotlib
        Where to do the plotting.
    series : dict
        Represents the data to be plotted.
    low : int, optional
        Lower index to start plotting. `low` can be negative.
    high : int, optional
        Index where to end the plotting, this index is not plotted.
        `high` is assumed to always be > 0 or None.
    color : string, optional
        A string representing the color to use.

    Returns
    -------
    handle : object like :py:class:`matplotlib.lines.Line2D`
        A handle for the plotted line.

    """
    # pick out just a few keys - we want to limit what we change here:
    kwargs = {'linestyle': series.get('ls', '-'),
              'alpha': series.get('alpha', 1.0),
              'linewidth': series.get('lw', 2.0)}
    if color is not None:
        kwargs['color'] = color
    elif 'color' in series:
        kwargs['color'] = series['color']

    if high is None:
        high = min(len(series['x']), len(series['y']))
    else:
        high = min(len(series['x']), len(series['y']), high)
    handle, = axs.plot(series['x'][low:high], series['y'][low:high],
                       **kwargs)
    return handle


def mpl_simple_plot(series, fig_settings=None):
    """Plot simple time series data (i.e. ``x`` vs ``y`` data).

    Parameters
    ----------
    series : list of dicts
        `series[i]` is the dict which contains the data to be plotted.
    fig_settings : dict, optional
        This dict contains settings for the figure, keys are:

        * `xlabel`: string, the label to use for the x-axis.
        * `ylabel`: string, the label to use for the y-axis.
        * `title`: string, title to use for the figure.
        * `yscale`: string, to change the scale for the y-axis.

    Returns
    -------
    out : object like :class:`matplotlib.backend_bases.FigureCanvasBase`
        This is the figure we create here.

    """
    fig = Figure()
    canvas = FigureCanvas(fig)
    axs = fig.add_subplot(111)
    handles = []
    labels = []
    for seri in series:
        handle = None
        if seri['type'] == 'xy':
            handle = mpl_plot_in_chunks(axs, seri)
        elif seri['type'] == 'bar':
            handle = mpl_plot_bar(axs, seri)
        elif seri['type'] == 'vline':
            handle = axs.axvline(x=seri['x'], ls=seri.get('ls', '-'),
                                 alpha=seri.get('alpha', 1.0),
                                 lw=seri.get('lw', 2.0))
        elif seri['type'] == 'hline':
            handle = axs.axhline(y=seri['y'], ls=seri.get('ls', '-'),
                                 alpha=seri.get('alpha', 1.0),
                                 lw=seri.get('lw', 2.0))
        legend = seri.get('label', None)
        if legend is not None and handle is not None:
            handles.append(handle)
            labels.append(legend)
    if 'xlabel' in fig_settings:
        axs.set_xlabel(fig_settings['xlabel'])
    if 'ylabel' in fig_settings:
        axs.set_ylabel(fig_settings['ylabel'])
    if 'title' in fig_settings:
        axs.set_title(fig_settings['title'], **_TITLE_SETTINGS)
    if len(labels) == len(handles) and len(labels) >= 1:
        ncol, rest = divmod(len(labels), 10)
        if rest > 0:
            ncol += 1
        if ncol <= 2:
            axs.legend(handles, labels, ncol=ncol)
    if 'yscale' in fig_settings:
        axs.set_yscale(fig_settings['yscale'])
    return canvas


def mpl_linecollection_gradient(axs, series):
    """Plot a line with a color gradient along the line.

    The line is split into segments and each segment is colored
    according to it's "location" along the line. By "location" we here
    mean the order of the segments.

    Parameters
    ----------
    axs : Axes object from matplotlib
        Where to do the plotting.
    series : dict
        Represents the data to be plotted.

    Returns
    -------
    handle : object like :class:`matplotlib.collections.LineCollection`
        A handle for the plotted line.

    """
    # pick out just a few keys - we want to limit what we change here:
    kwargs = {'linestyle': series.get('ls', '-'),
              'alpha': series.get('alpha', 1.0),
              'linewidth': series.get('lw', 2.0)}
    points = np.array([series['x'], series['y']]).T.reshape((-1, 1, 2))
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    values = np.linspace(0, 1, len(series['x']))
    linec = LineCollection(segments, array=values,
                           norm=matplotlib.colors.Normalize(vmin=0, vmax=1),
                           **kwargs)
    handle = axs.add_collection(linec)
    return handle


def mpl_chunks_gradient(axs, series, chunksize=20000):
    """Plot a line gradient in chunks.

    Here we will plot a line in chunks and color each chunk with one
    color. This function can be used as an alternative to
    :py:func:`.mpl_linecollection_gradient` when the number of points
    to plot is very large. Typically the chunk size here will be small
    compared to the size of the data to be plotted, so that each chunk,
    if plotted with :py:func:`.mpl_linecollection_gradient`, would have
    approximately the same color anyway.

    Parameters
    ---------
    axs : Axes object from matplotlib
        Where to do the plotting.
    series : dict
        Represents the data to be plotted.
    chunksize : int, optional
        This is the maximum size (number of points) we will try to plot
        in one go.

    Returns
    -------
    handle : object like :class:`matplotlib.lines.Line2D`
        A handle for the plotted line.

    Note
    ----
    Color maps in matplotlib will typically have 256 colors. The number
    of different colors we can get is currently limited to 256.

    """
    kwargs = {'linestyle': series.get('ls', '-'),
              'alpha': series.get('alpha', 1.0),
              'linewidth': series.get('lw', 2.0)}
    line = None
    lenx = len(series['x'])
    nchunk, rest = divmod(lenx, chunksize)
    cnorm = matplotlib.colors.Normalize(vmin=0, vmax=nchunk+1)
    color_map = matplotlib.cm.ScalarMappable(norm=cnorm)
    for i in range(nchunk):
        low = i * chunksize
        high = low + chunksize
        line, = axs.plot(series['x'][low:high], series['y'][low:high],
                         **kwargs)
        line.set_color(color_map.to_rgba(i))
    if rest > 0:
        last = rest + 1
        line, = axs.plot(series['x'][-last:], series['y'][-last:],
                         **kwargs)
        line.set_color(color_map.to_rgba(nchunk))
    return line


def mpl_line_gradient(series, fig_settings):
    """Plot time series with a color gradient.

    This function will plot time series data and color the lines with
    a gradient according to 'time'

    Parameters
    ----------
    series : list of dicts
        `series[i]` is the dict which contains the data to be plotted.
    fig_settings : dict
        This dict contains settings for the figure, keys are:

        * `xlabel`: string, the label to use for the x-axis.
        * `ylabel`: string, the label to use for the y-axis.
        * `title`: string, title to use for the figure.
        * `yscale`: string, to change the scale for the y-axis.

    Returns
    -------
    out : object like :class:`matplotlib.backend_bases.FigureCanvasBase`
        This is the figure we create here.

    Notes
    -----
    This function is based on the matplotlib example from:
    http://matplotlib.org/examples/pylab_examples/multicolored_line.html

    """
    fig = Figure()
    canvas = FigureCanvas(fig)
    axs = fig.add_subplot(111)
    handles = []
    labels = []
    for seri in series:
        lenx = len(seri['x'])
        if lenx >= 10**6:  # plot in chunks
            logger.info('Line gradient: LARGE dataset - plotting in chunks')
            handle = mpl_chunks_gradient(axs, seri)
        else:  # just plot it all
            handle = mpl_linecollection_gradient(axs, seri)
        legend = seri.get('label', None)
        if legend is not None and handle is not None:
            handles.append(handle)
            labels.append(legend)
    axs.autoscale_view()
    if 'xlabel' in fig_settings:
        axs.set_xlabel(fig_settings['xlabel'])
    if 'ylabel' in fig_settings:
        axs.set_ylabel(fig_settings['ylabel'])
    if 'title' in fig_settings:
        axs.set_title(fig_settings['title'], **_TITLE_SETTINGS)
    if len(labels) == len(handles) and len(labels) >= 1:
        axs.legend(handles, labels)
    return canvas


def mpl_error_plot(series, fig_settings):
    """Plot series with errors.

    Plot a series with error values displayed as a filled region.

    Parameters
    ----------
    series : list of tuples
        `series[i]` is the tuple which will be plotted. It is assumed
        to be on the form (x-values, y-values, y-error, legend).
    fig_settings : dict
        This dict contains settings for the figure, keys are:

        * `xlabel`: string, the label to use for the x-axis.
        * `ylabel`: string, the label to use for the y-axis.
        * `title`: string, title to use for the figure.

    Returns
    -------
    out : object like :class:`matplotlib.backend_bases.FigureCanvasBase`
        This is the figure we create here.

    """
    fig = Figure()
    canvas = FigureCanvas(fig)
    axs = fig.add_subplot(111)
    handles = []
    labels = []
    for seri in series:
        try:
            add_legend = seri[3] is not None
        except IndexError:
            add_legend = False
        handle, = axs.plot(seri[0], seri[1])
        axs.fill_between(seri[0], seri[1] + seri[2],
                         seri[1] - seri[2],
                         facecolor=handle.get_color(), alpha=0.3)
        if add_legend:
            handles.append(handle)
            labels.append(seri[3])
    if 'xlabel' in fig_settings:
        axs.set_xlabel(fig_settings['xlabel'])
    if 'ylabel' in fig_settings:
        axs.set_ylabel(fig_settings['ylabel'])
    if 'title' in fig_settings:
        axs.set_title(fig_settings['title'], **_TITLE_SETTINGS)
    if len(labels) == len(handles) and len(labels) >= 1:
        axs.legend(handles, labels)
    return canvas


def _mpl_shoots_histogram(histograms, ensemble):
    """Plot the histograms from the shoots analysis.

    Parameters
    ----------
    histograms : dict
        These are the histograms obtained in the shoots analysis.
    ensemble : string
        This is the ensemble identifier, e.g. 001, 002, etc.

    Returns
    -------
    out : object like :class:`matplotlib.backend_bases.FigureCanvasBase`
        This is the unscaled histogram.

    """
    series = []
    for key in ['ACC']:
        try:
            mid = histograms[key][2]
            hist = histograms[key][0]
            series.append({'type': 'xy', 'x': mid, 'y': hist,
                           'alpha': 0.8})
        except KeyError:
            continue
    figset = {'xlabel': 'Order parameter', 'ylabel': 'Probability density',
              'title': f'Ensemble ${ensemble}$'}
    canvas = mpl_simple_plot(series, fig_settings=figset)
    return canvas


def mpl_plot_path(results, path_ensemble):
    """Plot all figures from the path analysis.

    Parameters
    ----------
    results : dict
        This dict contains the result from the analysis.
    path_ensemble : object like :py:class:`.PathEnsemble`
        This is the path ensemble we have analyzed.

    Returns
    -------
    canvas : dict
        This dictionary contains the different canvases we have
        created.

    """
    ens = path_ensemble.ensemble_name
    ens_simplified = path_ensemble.ensemble_name_simple
    detect = results['detect']
    canvas = {}
    out = {}
    for key, item in PATHFILES.items():
        out[key] = item.format(ens_simplified)
    if 'pcross' in results:
        # First plot `pcross` vs `lambda` with the `detect` surface:
        series = [
            {'type': 'xy',
             'x': results['pcross'][0],
             'y': results['pcross'][1]},
            {'type': 'vline',
             'x': detect,
             'ls': '--',
             'alpha': 0.8},
        ]
        figset = {'xlabel': r'Order parameter ($\lambda$)',
                  'ylabel': 'Crossing probability',
                  'title': f'Ensemble ${ens}$'}
        canvas[out['pcross']] = mpl_simple_plot(series, fig_settings=figset)
    if 'prun' in results:
        # Next plot running ` pcross`:
        series = [
            {'type': 'xy',
             'x': results['cycle'],
             'y': results['prun']},
            {'type': 'hline',
             'y': results['prun'][-1],
             'ls': '--',
             'alpha': 0.8},
        ]
        figset = {'xlabel': 'Cycle number',
                  'ylabel': 'Running estimate',
                  'title': f'Ensemble ${ens}$'}
        canvas[out['prun']] = mpl_simple_plot(series, fig_settings=figset)
    if 'blockerror' in results:
        # Plot results of block-error analysis:
        series = [
            {'type': 'xy',
             'x': results['blockerror'][0],
             'y': results['blockerror'][3]},
            {'type': 'hline',
             'y': results['blockerror'][4],
             'ls': '--',
             'alpha': 0.8},
        ]
        title = r'Ensemble ${0}$: Rel. err.: {1:9.6e}, Ncor: {2:9.6f}'
        figset = {'xlabel': 'Block length', 'ylabel': 'Estimated error',
                  'title': title.format(ens, results['blockerror'][4],
                                        results['blockerror'][6])}
        canvas[out['perror']] = mpl_simple_plot(series, fig_settings=figset)
    if 'lengtherror' in results:
        # Plot results of length-error analysis:
        series = [
            {'type': 'xy',
             'x': results['lengtherror'][0],
             'y': results['lengtherror'][3]},
            {'type': 'hline',
             'y': results['lengtherror'][4],
             'ls': '--',
             'alpha': 0.8},
        ]
        title = r'Ensemble ${0}$: Rel. err.: {1:9.6e}, Ncor: {2:9.6f}'
        figset = {'xlabel': 'Block length', 'ylabel': 'Estimated error',
                  'title': title.format(ens, results['lengtherror'][4],
                                        results['lengtherror'][6])}
        canvas[out['lengtherror']] = mpl_simple_plot(series,
                                                     fig_settings=figset)

    # Plot length-histogram:
    labfmt = r'Avg. {0:6.2f} st.dev. {1:6.2f}'
    series = [
        {'type': 'xy',
         'x': results['pathlength'][0][1],
         'y': results['pathlength'][0][0],
         'label': labfmt.format(results['pathlength'][0][2][0],
                                results['pathlength'][0][2][1])},
    ]
    figset = {'xlabel': 'No. of PyRETIS engine steps', 'ylabel': 'Probability',
              'title': f'Ensemble ${ens}$'}
    canvas[out['pathlength']] = mpl_simple_plot(series, fig_settings=figset)
    # Plot shoots-histogram
    can_tmp = _mpl_shoots_histogram(results['shoots'][0], ens)
    canvas[out['shoots']] = can_tmp
    return canvas


def mpl_plot_pppath(results, path_ensemble):
    """Plot all figures from the repptis path analysis.

    Parameters
    ----------
    results : dict
        This dict contains the result from the analysis.
    path_ensemble : object like :py:class:`.PathEnsemble`
        This is the path ensemble we have analyzed.

    Returns
    -------
    canvas : dict
        This dictionary contains the different canvases we have
        created.

    """
    ens = path_ensemble.ensemble_name
    ens_simplified = path_ensemble.ensemble_name_simple
    canvas = {}
    out = {}
    colors = get_color_map(10)
    for key, item in PATHFILES.items():
        out[key] = item.format(ens_simplified)
    if 'pp_Pc' in results:
        # First plot `pcross` vs `lambda` with the `detect` surface:
        series = [
            {'type': 'xy',
             'x': list(range(len(results['pp_Pc']))),
             'y': results['pp_Pc']},
        ]
        figset = {'xlabel': r'intf',
                  'ylabel': r'$P_a(\lambda_i/\lambda_A)$',
                  'title': r''}
        canvas[f'{ens_simplified}_pp_Pc'] = mpl_simple_plot(
            series, fig_settings=figset)
    if 'pcross' in results:
        intf_l, intf_m, intf_r, lmlpercs, lmllambs, rmrpercs, rmrlambs = (
            results['repptis'])
        # First plot `pcross` vs `lambda` with the `detect` surface:
        series = [
            {'type': 'xy',
             'x': lmllambs,
             'y': lmlpercs,
             'color': colors[0]},
            {'type': 'xy',
             'x': rmrlambs,
             'y': rmrpercs,
             'color': colors[1]},
            {'type': 'vline',
             'x': intf_l,
             'ls': '--',
             'alpha': 0.8},
            {'type': 'vline',
             'x': intf_m,
             'ls': '--',
             'alpha': 0.8},
            {'type': 'vline',
             'x': intf_r,
             'ls': '--',
             'alpha': 0.8},
        ]
        figset = {'xlabel': r'Order parameter ($\lambda$)',
                  'ylabel': 'Crossing probability',
                  'title': f'Ensemble ${ens}$'
                           + f', L = {intf_l}, M = {intf_m}, R = {intf_r}'}
        canvas[out['pcross']] = mpl_simple_plot(series, fig_settings=figset)
    if 'prun_sl' in results:
        # Next plot running ` pcross`:
        series = [
            {'type': 'xy',
             'x': results['cycle'],
             'y': results['prun_sl'],
             'color': colors[0]},
            {'type': 'hline',
             'y': results['prun_sl'][-1],
             'ls': '--',
             'alpha': 0.8},
        ]
        figset = {'xlabel': 'Cycle number',
                  'ylabel': 'Running estimate',
                  'title': f'Ensemble ${ens}$'
                           + f", {results['prun_sl'][-1]}"}
        canvas[out['prun_sl']] = mpl_simple_plot(series, fig_settings=figset)
    if 'prun_sr' in results:
        # Next plot running ` pcross`:
        series = [
            {'type': 'xy',
             'x': results['cycle'],
             'y': results['prun_sr'],
             'color': colors[1]},
            {'type': 'hline',
             'y': results['prun_sr'][-1],
             'ls': '--',
             'alpha': 0.8},
        ]
        figset = {'xlabel': 'Cycle number',
                  'ylabel': 'Running estimate',
                  'title': f'Ensemble ${ens}$'
                           + f", {results['prun_sr'][-1]}"}
        canvas[out['prun_sr']] = mpl_simple_plot(series, fig_settings=figset)
    if 'blockerror_sl' in results:
        # Plot results of block-error analysis:
        series = [
            {'type': 'xy',
             'x': results['blockerror_sl'][0],
             'y': results['blockerror_sl'][3],
             'color': colors[0]},
            {'type': 'hline',
             'y': results['blockerror_sl'][4],
             'ls': '--',
             'alpha': 0.8},
        ]
        title = r'Ensemble ${0}$: Rel. err.: {1:9.6e}, Ncor: {2:9.6f}'
        figset = {'xlabel': 'Block length', 'ylabel': 'Estimated error',
                  'title': title.format(ens, results['blockerror_sl'][4],
                                        results['blockerror_sl'][6])}
        canvas[out['perror_sl']] = mpl_simple_plot(series, fig_settings=figset)
    if 'blockerror_sr' in results:
        # Plot results of block-error analysis:
        series = [
            {'type': 'xy',
             'x': results['blockerror_sr'][0],
             'y': results['blockerror_sr'][3],
             'color': colors[1]},
            {'type': 'hline',
             'y': results['blockerror_sr'][4],
             'ls': '--',
             'alpha': 0.8},
        ]
        title = r'Ensemble ${0}$: Rel. err.: {1:9.6e}, Ncor: {2:9.6f}'
        figset = {'xlabel': 'Block length', 'ylabel': 'Estimated error',
                  'title': title.format(ens, results['blockerror_sr'][4],
                                        results['blockerror_sr'][6])}
        canvas[out['perror_sr']] = mpl_simple_plot(series, fig_settings=figset)
    if 'lengtherror' in results:
        # Plot results of length-error analysis:
        series = [
            {'type': 'xy',
             'x': results['lengtherror'][0],
             'y': results['lengtherror'][3]},
            {'type': 'hline',
             'y': results['lengtherror'][4],
             'ls': '--',
             'alpha': 0.8},
        ]
        title = r'Ensemble ${0}$: Rel. err.: {1:9.6e}, Ncor: {2:9.6f}'
        figset = {'xlabel': 'Block length', 'ylabel': 'Estimated error',
                  'title': title.format(ens, results['lengtherror'][4],
                                        results['lengtherror'][6])}
        canvas[out['lengtherror']] = mpl_simple_plot(series,
                                                     fig_settings=figset)

    # Plot length-histogram:
    labfmt = r'Avg8 {0:6.2f} st.dev. {1:6.2f}'
    series = [
        {'type': 'xy',
         'x': results['pathlength'][0][1],
         'y': results['pathlength'][0][0],
         'label': labfmt.format(results['pathlength'][0][2][0],
                                results['pathlength'][0][2][1])},
    ]
    figset = {'xlabel': 'No. of PyRETIS engine steps',
              'ylabel': 'Crossing probability', 'title': f'Ensemble ${ens}$'}
    canvas[out['pathlength']] = mpl_simple_plot(series, fig_settings=figset)
    # Plot shoots-histogram
    can_tmp = _mpl_shoots_histogram(results['shoots'][0], ens)
    canvas[out['shoots']] = can_tmp
    return canvas


def mpl_plot_orderp(results, orderdata):
    """Plot the output from the order parameter analysis using matplotlib.

    Parameters
    ----------
    results : dict
        Each item in `results` contains the results for the
        corresponding order parameter.
    orderdata : list of numpy.arrays
        This is the raw data for the order parameter analysis.

    Returns
    -------
    canvas : dict
        The different plots created by this function.

    Note
    ----
    We are here only outputting results for the first order parameter.
    I.e. other order parameters or velocities are not written here. This
    will be changed when the structure of the output order parameter
    file has been fixed. Also note that, if present, the first order
    parameter will be plotted against the second one - i.e. the second
    one will be assumed to represent the velocity here.

    """
    canvas = {}
    time = orderdata[:, 0]
    series = [{'type': 'xy', 'x': time, 'y': orderdata[:, 1]}]
    figset = {'xlabel': 'Time', 'ylabel': 'Order parameter'}
    canvas[ORDERFILES['order']] = mpl_simple_plot(series, fig_settings=figset)
    # make running average plot of the energies as function of time
    series = [{'type': 'xy', 'x': time, 'y': results[0]['running'],
               'label': 'Running average'}]
    canvas[ORDERFILES['run_order']] = mpl_simple_plot(series,
                                                      fig_settings=figset)
    # plot block-error results:
    block = results[0]['blockerror']
    series = [
        {'type': 'xy', 'x': block[0], 'y': block[3]},
        {'type': 'hline', 'y': block[4], 'ls': '--', 'alpha': 0.8}
    ]
    title = 'Order parameter. Rel. err.: {0:9.6e}, Ncor: {1:9.6f}'
    figset = {'xlabel': 'Block length',
              'ylabel': 'Estimated error',
              'title': title.format(block[4], block[6])}
    canvas[ORDERFILES['block']] = mpl_simple_plot(series, fig_settings=figset)
    # plot distributions
    dist = results[0]['distribution']
    series = [{'type': 'xy', 'x': dist[1], 'y': dist[0]}]
    title = '{0}. Average: {1:9.6e}, std: {2:9.6f}'
    title = title.format('Order parameter', dist[2][0], dist[2][1])
    figset = {'title': title}
    canvas[ORDERFILES['dist']] = mpl_simple_plot(series,
                                                 fig_settings=figset)
    # also try a orderp vs ordervel plot:
    _, col = orderdata.shape
    if col >= 3:
        series = [{'type': 'xyc', 'x': orderdata[:, 1], 'y': orderdata[:, 2]}]
        figset = {'xlabel': r'$\lambda$',
                  'ylabel': r'$\dot{\lambda}$',
                  'title': 'Order parameter vs velocity'}
        plot = ORDERFILES['ordervel']
        canvas[plot] = mpl_line_gradient(series, fig_settings=figset)
    # output msd if it was calculated:
    if 'msd' in results[0]:
        msd = results[0]['msd']
        series = [(np.arange(len(msd)), msd[:, 0], msd[:, 1])]
        figset = {'xlabel': 'Time', 'ylabel': 'MSD'}
        canvas[ORDERFILES['msd']] = mpl_error_plot(series,
                                                   fig_settings=figset)
    return canvas


def mpl_plot_energy(results, energies):
    """Plot the output from the energy analysis using matplotlib.

    Parameters
    ----------
    results : dict
        Each item in `results` contains the results for the
        corresponding energy. It is assumed to contains the keys
        'vpot', 'ekin', 'etot', 'ham', 'temp', 'elec'.
    energies : dict of numpy.arrays
        This is the raw data for the energy analysis.

    Returns
    -------
    canvas : dict
        The output figures created by this function.

    """
    canvas = {}
    time = energies['time']
    # make time series plot of the energies
    series = []
    for key in ['vpot', 'ekin', 'etot', 'ham']:
        if key not in energies:
            continue
        series.append({'type': 'xy', 'x': time, 'y': energies[key],
                       'label': ENERTITLE[key]})
    figset = {'xlabel': 'Time', 'ylabel': 'Energy'}
    canvas[ENERFILES['energies']] = mpl_simple_plot(series,
                                                    fig_settings=figset)
    # make running average plot of the energies as function of time
    series = []
    for key in ['vpot', 'ekin', 'etot', 'ham']:
        if key not in results:
            continue
        series.append({'type': 'xy', 'x': time,
                       'y': results[key]['running'],
                       'label': ENERTITLE[key]})
    canvas[ENERFILES['run_energies']] = mpl_simple_plot(series,
                                                        fig_settings=figset)
    # plot temperature
    series = [{'type': 'xy', 'x': time, 'y': energies['temp']}]
    figset = {'xlabel': 'Time', 'ylabel': 'Temperature'}
    canvas[ENERFILES['temperature']] = mpl_simple_plot(series,
                                                       fig_settings=figset)
    # and running average for temperature
    series = [{'type': 'xy', 'x': time, 'y': results['temp']['running']}]
    figset = {'xlabel': 'Time',
              'ylabel': 'Temperature',
              'title': 'Running average'}
    canvas[ENERFILES['run_temp']] = mpl_simple_plot(series,
                                                    fig_settings=figset)
    # plot block-error results:
    title = r'{0}: Rel. err.: {1:9.6e}, Ncor: {2:9.6f}'
    for key in ['vpot', 'ekin', 'etot', 'temp']:
        if key not in results:
            continue
        plot = ENERFILES['block'].format(key)
        block = results[key]['blockerror']
        series = [
            {'type': 'xy', 'x': block[0], 'y': block[3]},
            {'type': 'hline', 'y': block[4], 'ls': '--', 'alpha': 0.8}
        ]
        figset = {'xlabel': 'Block length',
                  'ylabel': 'Estimated error',
                  'title': title.format(ENERTITLE[key], block[4], block[6])}
        canvas[plot] = mpl_simple_plot(series, fig_settings=figset)
    # plot distributions
    for key in ['vpot', 'ekin', 'etot', 'temp']:
        if key not in results:
            continue
        dist = results[key]['distribution']
        series = [{'type': 'xy', 'x': dist[1], 'y': dist[0],
                   'label': ENERTITLE[key]}]
        title = '{0}. Average: {1:9.6e}, std: {2:9.6f}'
        title = title.format(ENERTITLE[key], dist[2][0], dist[2][1])
        if 'boltzmann-dist' in results[key]:
            series.append({'type': 'xy',
                           'x': results[key]['boltzmann-dist'][1],
                           'y': results[key]['boltzmann-dist'][0],
                           'label': 'Boltzmann distribution'})
        plot = ENERFILES['dist'].format(key)
        canvas[plot] = mpl_simple_plot(series, fig_settings={'title': title})
    return canvas


def mpl_plot_flux(results):
    """Plot the output from the flux analysis using matplotlib.

    Parameters
    ----------
    results : dict
        This is the dict with the results from the flux analysis.

    Returns
    -------
    out[0] : list of dicts
        The output figures created by this function for running
        averages. `out[0][i]['name']` is the name of the figure and
        `out[0][i]['canvas']` is the corresponding canvas object.
    out[1] : list of dicts
        The output figures created by this function for block errors.
        `out[0][i]['name']` is the name of the figure and
        `out[0][i]['canvas']` is the corresponding canvas object.

    """
    canvas_run = []
    canvas_err = []
    for i in range(len(results['flux'])):
        # Plot running average:
        flux = results['flux'][i]
        runflux = results['runflux'][i]
        series = [{'type': 'xy', 'x': flux[:, 0], 'y': runflux,
                   'label': 'Running average'}]
        title = f'Flux for interface no. {i + 1}'
        figset = {'xlabel': 'Time',
                  'ylabel': 'Flux / internal units',
                  'title': title}
        canvas = mpl_simple_plot(series, fig_settings=figset)
        canvas_run.append({'name': FLUXFILES['runflux'].format(i + 1),
                           'canvas': canvas})
        # Plot error results:
        errflux = results['errflux'][i]
        title = r'Block error: {0}: Rel. err.: {1:9.6e}, Ncor: {2:9.6f}'
        series = [{'type': 'xy', 'x': errflux[0], 'y': errflux[3]},
                  {'type': 'hline', 'y': errflux[4], 'ls': '--',
                   'alpha': 0.8}]
        figset = {'xlabel': 'Block length',
                  'ylabel': 'Estimated error',
                  'title': title.format(i + 1, errflux[4], errflux[6])}
        canvas = mpl_simple_plot(series, fig_settings=figset)
        canvas_err.append({'name': FLUXFILES['block'].format(i + 1),
                           'canvas': canvas})
    return canvas_run, canvas_err


def mpl_plot_xi(xi, xierror):
    """Plot the output from the Xi analysis using matplotlib.

    Parameters
    ----------
    xi: array-like
        The running average of xi
    xierror: tuple
        The results from the blockerror analysis of xi

    Returns
    -------
    out[0] : list of dicts
        The output figures created by this function for running
        averages. `out[0][i]['name']` is the name of the figure and
        `out[0][i]['canvas']` is the corresponding canvas object.
    out[1] : list of dicts
        The output figures created by this function for block errors.
        `out[0][i]['name']` is the name of the figure and
        `out[0][i]['canvas']` is the corresponding canvas object.

    """
    canvas_run = []
    canvas_err = []
    series = [{'type': 'xy', 'x': range(len(xi)), 'y': xi,
               'label': 'Running average'},
              {'type': 'hline', 'y': xi[-1], 'ls': '--', 'alpha': 0.8}]
    title = r'$\xi={}$:'
    figset = {'xlabel': 'cycles',
              'ylabel': r'$\xi$',
              'title': title.format(xi[-1])}
    canvas = mpl_simple_plot(series, fig_settings=figset)
    canvas_run.append({'name': "xirun",
                       'canvas': canvas})
    # Plot error results:
    title = r"relative $\xi$ error: Rel. err.: {0:9.6e}, Ncor: {1:9.6f}"
    series = [{'type': 'xy', 'x': xierror[0], 'y': xierror[3]},
              {'type': 'hline', 'y': xierror[4], 'ls': '--',
               'alpha': 0.8}]
    figset = {'xlabel': 'Block length',
              'ylabel': 'Estimated error',
              'title': title.format(xierror[4], xierror[6])}
    canvas = mpl_simple_plot(series, fig_settings=figset)
    canvas_err.append({'name': "xierr",
                       'canvas': canvas})

    return canvas_run, canvas_err


def mpl_plot_tau(tau, tauerror):
    """Plot the output from the Tau analysis using matplotlib.

    Parameters
    ----------
    tau: array-like
        The running average of tau.
    tauerror: tuple
        The results from the blockerror analysis of tau.

    Returns
    -------
    out[0] : list of dicts
        The output figures created by this function for running
        averages. `out[0][i]['name']` is the name of the figure and
        `out[0][i]['canvas']` is the corresponding canvas object.
    out[1] : list of dicts
        The output figures created by this function for block errors.
        `out[0][i]['name']` is the name of the figure and
        `out[0][i]['canvas']` is the corresponding canvas object.

    """
    canvas_run = []
    canvas_err = []
    if not tau.any():
        str_in = "no"
    else:
        str_in = "in"
    series = [{'type': 'xy', 'x': range(len(tau)), 'y': tau,
               'label': 'Running average'},
              {'type': 'hline', 'y': tau[-1], 'ls': '--', 'alpha': 0.8}]
    title = (r'$\frac{\tau}{dz}$ '
             f"({str_in} "
             r'selected $\tau$ region)$= $')
    figset = {'xlabel': 'cycles',
              'ylabel': r'$\frac{\tau}{dz}$',
              'title': title+f"${tau[-1]}$"}
    canvas = mpl_simple_plot(series, fig_settings=figset)
    canvas_run.append({'name': "taurun",
                       'canvas': canvas})
    # Plot error results:
    title = (r"relative error $\frac{\tau}{dz}$ "
             f"({str_in} "
             r"selected $\tau$ region):")
    title2 = r"Rel. err.: {0:9.6e}, Ncor: {1:9.6f}"
    series = [{'type': 'xy', 'x': tauerror[0], 'y': tauerror[3]},
              {'type': 'hline', 'y': tauerror[4], 'ls': '--',
               'alpha': 0.8,
               'label': title2.format(tauerror[4], tauerror[6])}]
    figset = {'xlabel': 'Block length',
              'ylabel': 'Estimated error',
              'title': title}
    canvas = mpl_simple_plot(series, fig_settings=figset)
    canvas_err.append({'name': "tauerr",
                       'canvas': canvas})

    return canvas_run, canvas_err


def mpl_plot_tau_ref(tau_ref, tau_ref_bins, bins):
    """Plot a histogram of the Tau reference bins using matplotlib.

    Parameters
    ----------
    tau_ref: array-like
        The final running average of all tau ref bins.
    tau_ref_bins: list of tuple
        The (L,R) bin edges of the reference bins.
    bins: tuple
        The (L,R) bin edges of the region used to calculate tau.

    Returns
    -------
    out : dict
        The output figure created by this function for histograming the
        tau reference bins. dict['name']` is the name of the figure and
        `dict['canvas']` is the corresponding canvas object.

    """
    x = [i[0] for i in tau_ref_bins]
    width = x[1]-x[0]
    series = [{'type': 'bar', 'x': x, 'height': tau_ref,
               'width': width, 'align': 'edge'}]
    if len(bins) == 2:
        series.append(
              {'type': 'vline', 'x': bins[0], 'ls': '--', 'alpha': 0.8,
               'label': r'selected $\tau$ region'}
        )
        series.append(
              {'type': 'vline', 'x': bins[1], 'ls': '--', 'alpha': 0.8}
        )
    title = r'$\frac{\tau}{dz}$ reference bins'
    figset = {'xlabel': 'Orderparameter',
              'ylabel': r'$\frac{\tau}{dz}$ count',
              'title': title}
    canvas = mpl_simple_plot(series, fig_settings=figset)
    return {"name": "tauref", "canvas": canvas}


def get_color_map(ncolors):
    """Return a color map with at least n colors."""
    if ncolors <= 10:
        name = 'tab10'
    elif 10 < ncolors <= 20:
        name = 'tab20'
    else:
        name = 'tab20'
    logger.info('Using color map %s', name)
    cmap = matplotlib.colormaps[name]
    return cmap(np.linspace(0, 1, ncolors))


def mpl_plot_matched(path_ensembles, detect, matched, reptis=False):
    """Plot matched probabilities using matplotlib.

    This function will plot the matched probabilities for the different
    ensembles and also make a plot with just the over-all matched
    probability.

    Parameters
    ----------
    path_ensembles : list of strings
        This is the name of the path ensembles we have calculated
        the probability for.
    detect : list of floats
        These are the detect interfaces used in the analysis.
    matched : dict
        This dict contains the results from the matching of the
        probabilities: `matched['overall-prob']`,
         and `matched['overall-prun']`
         or `matched['overall-rrun']` (if exists) are used here.

    Returns
    -------
    canvas : dict
        The output figures created by this function.

    """
    canvas = {}
    # Also make a plot with the overall matched probability:
    series = []
    for idetect in detect:
        series.append({'type': 'vline', 'x': idetect,
                       'ls': '--', 'alpha': 0.8, 'lw': 1})
    if reptis:
        pcross = matched['pcross']
        n = min(len(detect), len(pcross))
        matched['overall-prob'] = np.array([detect[:n], pcross[:n]]).T
    series.append({'type': 'xy',
                   'x': matched['overall-prob'][:, 0],
                   'y': matched['overall-prob'][:, 1],
                   'lw': 3})
    figset = {'xlabel': r'Order parameter ($\lambda$)',
              'ylabel': 'Crossing probability',
              'title': 'Matched probability',
              'yscale': 'log'}
    canvas[PATH_MATCH['match']] = mpl_simple_plot(series,
                                                  fig_settings=figset)
    # Also make a plot with the TIME evolution of the overall
    # probability or rate:
    if 'overall-rrun' in matched:
        title = 'Time evolution of the overall Rate'
        ylabel = 'Rate (running estimate)'
        data = matched['overall-rrun']
        label = 'overall-rrun'

    elif 'overall-prun' in matched:
        title = 'Time evolution of the overall crossing probability'
        ylabel = 'Crossing probability (running estimate)'
        data = matched['overall-prun']
        label = 'overall-prun'
    else:
        return canvas

    series = [{'type': 'xy',
               'x': matched['overall-cycle'],
               'y': data},
              {'type': 'hline',
               'y': data[-1],
               'ls': '--',
               'alpha': 0.8}]
    figset = {'xlabel': 'Cycle number',
              'ylabel': ylabel,
              'title': title}
    canvas[label] = mpl_simple_plot(series, fig_settings=figset)
    # and finish with the block error analysis.
    series = [{'type': 'xy',
               'x': matched['overall-error'][0],
               'y': matched['overall-error'][3]},
              {'type': 'hline',
               'y': matched['overall-error'][4],
               'ls': '--',
               'alpha': 0.8}]
    if reptis:
        title = 'Relative error of the running crossing probability'
    else:
        title = r'Overall: Rate Rel. err.: {0:9.6e}, Ncor: {1:9.6e}'
        title.format(matched['overall-error'][4], matched['overall-error'][6])
    figset = {'xlabel': 'Block length', 'ylabel': 'Estimated error',
              'title': title}
    canvas[PATH_MATCH['error']] = mpl_simple_plot(series,
                                                  fig_settings=figset)
    if reptis:
        # reptis does not deal with ..
        return canvas

    # First plot the matched probabilities for each ensemble:
    series = []
    for idetect in detect:
        series.append({'type': 'vline', 'x': idetect,
                       'ls': '--', 'alpha': 0.8, 'lw': 1})

    # color_cycle was deprecated in matplotlib 1.5, but to support old
    # versions:
    ckey = 'axes.prop_cycle'
    # Check if we need to have more colors:
    if len(matplotlib.rcParams[ckey]) < len(path_ensembles):
        logger.info('Overriding color cycle.')
        colors = get_color_map(len(path_ensembles))
    else:
        colors = None

    new_series = {'type': 'xy',
                  'x': matched['overall-prob'][:, 0],
                  'y': matched['overall-prob'][:, 1],
                  'alpha': 0.8,
                  'lw': 9, 'label': 'Overall', 'color': '#262626'}
    series.append(new_series)
    for i, (prob, path_e) in enumerate(zip(matched['matched-prob'],
                                           path_ensembles)):
        new_series = {'type': 'xy', 'x': prob[:, 0], 'y': prob[:, 1], 'lw': 3,
                      'label': f'${path_e}$'}
        if colors is not None:
            new_series['color'] = colors[i]
        series.append(new_series)

    figset = {'xlabel': r'Order parameter ($\lambda$)',
              'ylabel': 'Probability',
              'title': 'Overall probabilities',
              'yscale': 'log'}
    canvas[PATH_MATCH['total']] = mpl_simple_plot(series,
                                                  fig_settings=figset)
    return canvas


def mpl_plot_bar(axs, seri):
    """Draw a bar chart on the given axes using series parameters.

    Parameters
    ----------
    axs : matplotlib axes object
        The axes to draw on.
    seri : dict
        Series parameters passed to ``axs.bar`` (the 'type' key is removed).

    Returns
    -------
    out : matplotlib BarContainer
        The bar container returned by ``axs.bar``.
    """
    del seri['type']
    return axs.bar(**seri)

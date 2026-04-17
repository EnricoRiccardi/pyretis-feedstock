# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Compatibility tests for PyVisA's visualization helpers."""
import pytest
import pandas as pd
from matplotlib import colors
from matplotlib.figure import Figure

from pyretis.pyvisa import HAS_PYVISA_REQ


pytestmark = pytest.mark.skipif(not HAS_PYVISA_REQ,
                                reason="PyVisA reqs not installed")

if HAS_PYVISA_REQ:
    from pyretis.pyvisa.orderparam_density import Trajectory
    from pyretis.pyvisa.visualize import (
        DataSlave,
        _remove_last_line,
        _redraw_colorbar,
        _set_tick_fontsize,
    )


def _make_traj(cycle):
    """Create a compact trajectory with all variables needed by DataSlave."""
    frames = pd.DataFrame({
        'op1': [1.0 + 2 * cycle, 2.0 + 2 * cycle],
        'op2': [10.0 + cycle, 11.0 + cycle],
        'potE': [0.1 + cycle, 0.2 + cycle],
        'kinE': [0.3 + cycle, 0.4 + cycle],
        'totE': [0.4 + cycle, 0.6 + cycle],
    })
    info = {
        'stored': False,
        'reactive': True,
        'status': 'ACC',
        'length': len(frames.index),
        'cycle': cycle,
        'ensemble_name': '000',
    }
    return Trajectory(frames, info)


def _make_dataslave():
    """Create a DataSlave with in-memory trajectory data."""
    data = DataSlave()
    data.infos = {'op_labels': ['op1', 'op2']}
    data.traj_data = {'000': {0: _make_traj(0), 1: _make_traj(1)}}
    return data


def test_dataslave_dataframe_collection_without_append():
    """Collect analysis data on pandas versions without DataFrame.append."""
    assert not hasattr(pd.DataFrame, 'append')
    settings = {'stored': False, 'start_end': 'Reactive',
                'ACC': 'ACC', 'MC-move': 'All', 'fol': '000'}
    data = _make_dataslave()

    dataframe = data.load_to_df(settings)
    reactive = data.load_reactive_unreactive(settings)
    reduced, reactive_reduced = data.remove_values(settings, 1.5, 3.5)

    assert list(dataframe['op1']) == [1.0, 2.0, 3.0, 4.0]
    assert list(reactive[0]) == [True, True, True, True]
    assert list(reduced['op1']) == [2.0, 3.0]
    assert list(reactive_reduced) == [True, True]


def test_matplotlib_tick_label_compatibility():
    """Set tick fonts using current Matplotlib tick label attributes."""
    fig = Figure()
    ax = fig.add_subplot(111)
    _set_tick_fontsize(ax.xaxis, 17)

    for tick in ax.xaxis.get_major_ticks():
        assert tick.label1.get_fontsize() == 17


def test_matplotlib_tick_label_legacy_fallback():
    """Set tick fonts through the legacy tick.label fallback."""
    class Label:
        """Small label double for unit-testing fontsize updates."""

        def __init__(self):
            self.fontsize = None

        def set_fontsize(self, fontsize):
            self.fontsize = fontsize

    class Tick:
        """Small tick double exposing only the old label attribute."""

        def __init__(self):
            self.label = Label()

    class Axis:
        """Small axis double exposing major ticks."""

        def __init__(self):
            self.tick = Tick()

        def get_major_ticks(self):
            return [self.tick]

    axis = Axis()
    _set_tick_fontsize(axis, 13)

    assert axis.tick.label.fontsize == 13


def test_matplotlib_colorbar_redraw_compatibility():
    """Refresh colorbars on Matplotlib versions without draw_all."""
    fig = Figure()
    ax = fig.add_subplot(111)
    image = ax.imshow([[0.0, 1.0], [1.0, 0.0]])
    colorbar = fig.colorbar(image, ax=ax)

    image.set_norm(colors.Normalize(vmin=-1.0, vmax=2.0))
    _redraw_colorbar(colorbar)

    assert colorbar.mappable.norm.vmin == -1.0
    assert colorbar.mappable.norm.vmax == 2.0


def test_matplotlib_colorbar_redraw_draw_all_branch():
    """Use draw_all directly when an older Matplotlib colorbar exposes it."""
    class Colorbar:
        """Small colorbar double exposing the old draw_all API."""

        def __init__(self):
            self.called = False

        def draw_all(self):
            self.called = True

    colorbar = Colorbar()
    _redraw_colorbar(colorbar)

    assert colorbar.called


def test_matplotlib_artist_list_line_removal_compatibility():
    """Remove lines through the supported artist API."""
    fig = Figure()
    ax = fig.add_subplot(111)
    first = ax.plot([0.0, 1.0], [0.0, 1.0])[0]
    second = ax.plot([0.0, 1.0], [1.0, 0.0])[0]

    _remove_last_line(ax)

    assert list(ax.lines) == [first]
    assert second not in ax.lines


def test_matplotlib_artist_list_line_removal_empty_axes():
    """Removing from an empty axes is a no-op."""
    fig = Figure()
    ax = fig.add_subplot(111)

    _remove_last_line(ax)

    assert list(ax.lines) == []

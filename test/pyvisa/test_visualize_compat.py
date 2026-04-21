# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Compatibility tests for PyVisA's visualization helpers."""
import os

import pytest
import pandas as pd
from matplotlib import cm, colors
from matplotlib.figure import Figure

from pyretis.pyvisa import HAS_PYVISA_REQ


os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


pytestmark = pytest.mark.skipif(not HAS_PYVISA_REQ,
                                reason="PyVisA reqs not installed")

if HAS_PYVISA_REQ:
    import numpy as np
    from PyQt5 import QtWidgets
    from pyretis.pyvisa.orderparam_density import Trajectory
    from pyretis.pyvisa.visualize import (
        DataSlave,
        VisualApp,
        VisualObject,
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


def test_visualapp_option_handlers_before_plot_ready(monkeypatch):
    """Ignore GUI option changes until a figure has been initialized."""
    monkeypatch.setattr(VisualApp, "show", lambda self: None)
    monkeypatch.setattr(VisualApp, "center_on_screen", lambda self: None)

    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    window = VisualApp(folder=None, iofile=None,
                       sub_files={'rst': None, 'traj': None})

    window.toggle_titles()
    window.toggle_regl()
    window.toggle_intf()
    window._update_canvas_font()
    assert window.statusbar.currentMessage() == 'Figure not ready!'

    window._update_canvas_text()
    window._change_cmap()
    window._change_zoom()
    window.update_preview()
    assert window.statusbar.currentMessage() == 'Figure not ready!'

    window.toggle_only_stored_trjs()
    assert window.statusbar.currentMessage() == 'No data loaded'

    window.deleteLater()
    app.processEvents()


# ---------------------------------------------------------------------------
# DataSlave.find_cycle_ensemble
# ---------------------------------------------------------------------------

def test_dataslave_find_cycle_ensemble_no_attributes():
    """Return empty result when x/y attributes have been removed."""
    data = DataSlave()
    del data.x
    result, cycle, ensemble = data.find_cycle_ensemble(None)
    assert result == []
    assert cycle == ''
    assert ensemble == ''


def test_dataslave_find_cycle_ensemble_empty_lists():
    """Return empty result when no data points are stored."""
    data = DataSlave()
    data.x, data.y = [], []
    result, cycle, ensemble = data.find_cycle_ensemble(None)
    assert result == []
    assert cycle == ''
    assert ensemble == ''


def test_dataslave_find_cycle_ensemble_finds_nearest():
    """Return coordinates and metadata of the closest plotted point."""
    data = DataSlave()
    data.x = [0.0, 5.0, 10.0]
    data.y = [0.0, 5.0, 10.0]
    data.data_origin = [['000', 0], ['000', 1], ['000', 2]]

    class _Event:
        xdata, ydata = 4.9, 5.1

    closest, cycle, ensemble = data.find_cycle_ensemble(_Event())
    assert list(closest) == [5.0, 5.0]
    assert cycle == 1
    assert ensemble == '000'


# ---------------------------------------------------------------------------
# DataSlave.load_all
# ---------------------------------------------------------------------------

def _make_traj_for(cycle, ensemble_name):
    """Create a trajectory assigned to a specific ensemble."""
    frames = pd.DataFrame({
        'op1': [1.0 + 2 * cycle, 2.0 + 2 * cycle],
        'op2': [10.0 + cycle, 11.0 + cycle],
        'potE': [0.1 + cycle, 0.2 + cycle],
        'kinE': [0.3 + cycle, 0.4 + cycle],
        'totE': [0.4 + cycle, 0.6 + cycle],
    })
    info = {
        'stored': False, 'reactive': True, 'status': 'ACC',
        'length': len(frames.index), 'cycle': cycle,
        'ensemble_name': ensemble_name,
    }
    return Trajectory(frames, info)


_LOAD_ALL_SETTINGS = {
    'op1': 'op1', 'op2': 'op2', 'E': 'potE',
    'ACC': 'BOTH', 'MC-move': 'All', 'weight': False,
    'stored': False, 'start_end': 'All',
    'min_cycle': 0, 'max_cycle': 0,
    'fol': '000',
}


def test_dataslave_load_all_specific_ensemble():
    """load_all returns op coordinates from a named ensemble."""
    data = DataSlave()
    data.traj_data = {'000': {0: _make_traj_for(0, '000')}}
    data.infos = {'op_labels': ['op1', 'op2'], 'ensemble_names': ['000']}
    data.settings = dict(_LOAD_ALL_SETTINGS)

    x, y, z = data.load_all()
    assert list(x) == [1.0, 2.0]
    assert list(y) == [10.0, 11.0]
    assert list(z) == [0.1, 0.2]


def test_dataslave_load_all_reactive_filter():
    """load_all adds reactive criterion when start_end is not 'All'."""
    data = DataSlave()
    data.traj_data = {'000': {0: _make_traj_for(0, '000')}}
    data.infos = {'op_labels': ['op1', 'op2'], 'ensemble_names': ['000']}
    settings = dict(_LOAD_ALL_SETTINGS)
    settings['start_end'] = 'Reactive'
    data.settings = settings

    x, y, z = data.load_all()
    assert len(x) == 2


def test_dataslave_load_all_all_ensembles():
    """load_all aggregates data across all ensembles when fol='All'."""
    data = DataSlave()
    data.traj_data = {
        '000': {0: _make_traj_for(0, '000')},
        '001': {0: _make_traj_for(2, '001')},
    }
    data.infos = {'op_labels': ['op1', 'op2'],
                  'ensemble_names': ['000', '001']}
    settings = dict(_LOAD_ALL_SETTINGS)
    settings['fol'] = 'All'
    data.settings = settings

    x, y, z = data.load_all()
    assert len(x) == 4
    assert len(data.data_origin) == 4


def test_dataslave_load_all_uses_explicit_first_column_values():
    """Single-file mode can plot the stored first column, not only 0..N."""
    frames = pd.DataFrame({
        'time': [0.5, 1.5],
        'lambda': [2.0, 3.0],
    })
    info = {
        'stored': False, 'reactive': True, 'status': 'ACC',
        'length': len(frames.index), 'cycle': 0, 'ensemble_name': 'data',
    }
    data = DataSlave()
    data.traj_data = {'data': {0: Trajectory(frames, info)}}
    data.infos = {
        'op_labels': ['time', 'lambda'],
        'plot_labels': ['time', 'lambda'],
        'main_op_label': 'lambda',
        'ensemble_names': ['data'],
        'analysis_columns': 2,
    }
    settings = dict(_LOAD_ALL_SETTINGS)
    settings.update({'fol': 'data', 'op1': 'time', 'op2': 'lambda', 'E': 'None'})
    data.settings = settings

    x, y, z = data.load_all()

    assert list(x) == [0.5, 1.5]
    assert list(y) == [2.0, 3.0]
    assert list(z) == [1, 1]


# ---------------------------------------------------------------------------
# DataSlave.load_to_df — additional branches
# ---------------------------------------------------------------------------

def test_dataslave_load_to_df_all_ensembles():
    """load_to_df concatenates data from every ensemble when fol='All'."""
    data = DataSlave()
    data.infos = {'op_labels': ['op1', 'op2']}
    data.traj_data = {
        '000': {0: _make_traj(0)},
        '001': {0: _make_traj(1)},
    }
    settings = {'stored': False, 'start_end': 'Reactive',
                'ACC': 'ACC', 'MC-move': 'All', 'fol': 'All'}
    df = data.load_to_df(settings)
    assert len(df) == 4


def test_dataslave_load_to_df_acc_both():
    """load_to_df skips status filter when ACC='BOTH'."""
    data = _make_dataslave()
    settings = {'stored': False, 'start_end': 'Reactive',
                'ACC': 'BOTH', 'MC-move': 'All', 'fol': '000'}
    df = data.load_to_df(settings)
    assert len(df) == 4


def test_dataslave_load_to_df_mcmove_filter_no_match():
    """load_to_df returns empty df when MC-move filter finds no matches."""
    data = _make_dataslave()
    settings = {'stored': False, 'start_end': 'Reactive',
                'ACC': 'ACC', 'MC-move': 'sh', 'fol': '000'}
    df = data.load_to_df(settings)
    assert df.empty


def test_dataslave_load_to_df_column_count_mismatch():
    """load_to_df skips trajectories with unexpected column count."""
    frames = pd.DataFrame({'op1': [1.0, 2.0], 'op2': [3.0, 4.0]})
    info = {'stored': False, 'reactive': True, 'status': 'ACC',
            'length': 2, 'cycle': 0, 'ensemble_name': '000'}
    data = DataSlave()
    data.infos = {'op_labels': ['op1', 'op2']}
    data.traj_data = {'000': {0: Trajectory(frames, info)}}
    settings = {'stored': False, 'start_end': 'Reactive',
                'ACC': 'ACC', 'MC-move': 'All', 'fol': '000'}
    df = data.load_to_df(settings)
    assert df.empty


def test_dataslave_load_to_df_no_matches():
    """load_to_df returns empty df when selection criteria excludes all."""
    data = _make_dataslave()
    settings = {'stored': False, 'start_end': 'Unreactive',
                'ACC': 'REJ', 'MC-move': 'All', 'fol': '000'}
    df = data.load_to_df(settings)
    assert df.empty


# ---------------------------------------------------------------------------
# DataSlave.load_reactive_unreactive — additional branches
# ---------------------------------------------------------------------------

def test_dataslave_load_reactive_unreactive_mcmove_filter():
    """load_reactive_unreactive adds MC-move to criteria when not 'All'."""
    data = _make_dataslave()
    settings = {'stored': False, 'start_end': 'Reactive',
                'ACC': 'ACC', 'MC-move': 'sh', 'fol': '000'}
    reactive = data.load_reactive_unreactive(settings)
    assert reactive.empty


def test_dataslave_load_reactive_unreactive_status_mismatch():
    """load_reactive_unreactive skips mismatched-status trajectories."""
    data = _make_dataslave()
    settings = {'stored': False, 'start_end': 'Reactive',
                'ACC': 'REJ', 'MC-move': 'All', 'fol': '000'}
    reactive = data.load_reactive_unreactive(settings)
    assert reactive.empty


def test_dataslave_load_reactive_unreactive_all_ensembles():
    """load_reactive_unreactive collects flags from every ensemble."""
    data = DataSlave()
    data.infos = {'op_labels': ['op1', 'op2']}
    data.traj_data = {
        '000': {0: _make_traj(0)},
        '001': {0: _make_traj(1)},
    }
    settings = {'stored': False, 'start_end': 'Reactive',
                'ACC': 'ACC', 'MC-move': 'All', 'fol': 'All'}
    reactive = data.load_reactive_unreactive(settings)
    assert len(reactive) == 4
    assert list(reactive[0]) == [True, True, True, True]


# ---------------------------------------------------------------------------
# DataSlave.remove_values — edge cases
# ---------------------------------------------------------------------------

def test_dataslave_remove_values_all_filtered():
    """remove_values returns empty df when every point is out of range."""
    data = _make_dataslave()
    settings = {'stored': False, 'start_end': 'Reactive',
                'ACC': 'ACC', 'MC-move': 'All', 'fol': '000'}
    reduced, reactive_reduced = data.remove_values(settings, 0.0, 0.5)
    assert reduced.empty
    assert reactive_reduced.empty


# ---------------------------------------------------------------------------
# VisualApp._iffloat — static method
# ---------------------------------------------------------------------------

def test_visualapp_iffloat_valid_float():
    """_iffloat converts a numeric string to float."""
    assert VisualApp._iffloat('3.14', 0.0) == 3.14


def test_visualapp_iffloat_invalid_string():
    """_iffloat returns the default when conversion fails."""
    assert VisualApp._iffloat('not-a-number', -1.0) == -1.0


# ---------------------------------------------------------------------------
# VisualObject — construction
# ---------------------------------------------------------------------------

def test_visualobject_init_defaults():
    """VisualObject initializes with None pfile/trajfile and empty settings."""
    vo = VisualObject()
    assert vo.pfile is None
    assert vo.trajfile is None
    assert vo.settings == {}


def test_visualobject_init_with_paths():
    """VisualObject stores trajfile; pfile reset by PathVisualize.__init__."""
    vo = VisualObject(pfile='/data/run.hdf5', trajfile='/data/traj.xtc')
    assert vo.trajfile == '/data/traj.xtc'
    assert vo.settings == {}


# ---------------------------------------------------------------------------
# VisualApp.rev_trj — static method (no Qt required)
# ---------------------------------------------------------------------------

def test_visualapp_rev_trj():
    """rev_trj returns trajectory frames in reverse order (excluding first)."""
    import mdtraj as md

    top = md.Topology()
    chain = top.add_chain()
    res = top.add_residue('ALA', chain)
    top.add_atom('CA', md.element.carbon, res)
    xyz = np.zeros((3, 1, 3))
    xyz[:, 0, 0] = [1.0, 2.0, 3.0]
    traj = md.Trajectory(xyz, top)

    r = VisualApp.rev_trj(traj)
    assert list(r.xyz[:, 0, 0]) == [3.0, 2.0]


# ---------------------------------------------------------------------------
# VisualApp — shared window factory
# ---------------------------------------------------------------------------

def _make_window(monkeypatch):
    monkeypatch.setattr(VisualApp, 'show', lambda self: None)
    monkeypatch.setattr(VisualApp, 'center_on_screen', lambda self: None)
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    win = VisualApp(folder=None, iofile=None,
                    sub_files={'rst': None, 'traj': None})
    return app, win


def _make_settings(**overrides):
    """Build a compact settings dictionary for VisualApp tests."""
    settings = {
        'op1': 'op1',
        'op2': 'op2',
        'E': 'potE',
        'ACC': 'ACC',
        'fol': '000',
        'min_cycle': 0,
        'max_cycle': 1,
        'weight': False,
        'method': ['Density'],
        'dim': 2,
        'try_shift': False,
        'show_int': True,
        'reg_line': True,
        'res': 4,
        'colormap': 'viridis',
        'colorbar': True,
        'title font': 12,
        'axes font': 10,
        'show titles': True,
        'stored': False,
        'MC-move': 'All',
        'start_end': 'Reactive',
        'x-limits': (0.0, 1.0),
        'y-limits': (0.0, 1.0),
        'z-limits': (0.0, 1.0),
    }
    settings.update(overrides)
    return settings


class _DataObjectStub:
    """Small dataobject double for VisualApp helper tests."""

    def __init__(self, x, y, z, infos=None):
        self._xyz = (x, y, z)
        self.infos = infos or {
            'interfaces': [0.0, 0.5, 1.0],
            'op_labels': ['op1', 'op2'],
        }

    def load_all(self):
        return self._xyz


# ---------------------------------------------------------------------------
# VisualApp — save methods before data is loaded
# ---------------------------------------------------------------------------

def test_visualapp_save_methods_report_no_data(monkeypatch):
    """All save methods report 'No data selected' with no data."""
    app, win = _make_window(monkeypatch)
    try:
        for method in (win.save_hdf5, win.save_json, win.save_png,
                       win.save_textdata, win.save_script):
            method()
            assert win.statusbar.currentMessage() == 'No data selected'
    finally:
        win.deleteLater()
        app.processEvents()


# ---------------------------------------------------------------------------
# VisualApp — toggle_buttons
# ---------------------------------------------------------------------------

def test_visualapp_toggle_buttons(monkeypatch):
    """toggle_buttons enables/disables key UI controls."""
    app, win = _make_window(monkeypatch)
    try:
        win.toggle_buttons(True)
        assert win.updateBtn.isEnabled()
        assert win.saveFigBtn.isEnabled()
        assert win.refreshBtn.isEnabled()

        win.toggle_buttons(False)
        assert not win.updateBtn.isEnabled()
        assert not win.saveFigBtn.isEnabled()
        assert not win.refreshBtn.isEnabled()
    finally:
        win.deleteLater()
        app.processEvents()


# ---------------------------------------------------------------------------
# VisualApp — disable_zcombox / disable_reglinecheck
# ---------------------------------------------------------------------------

def test_visualapp_disable_zcombox(monkeypatch):
    """Energy combobox is disabled for Density plots, re-enabled otherwise."""
    app, win = _make_window(monkeypatch)
    try:
        win.plotTypeComBox.setCurrentText('2D Density')
        win.disable_zcombox()
        assert not win.energyComBox.isEnabled()

        win.plotTypeComBox.setCurrentText('2D Scatter')
        win.disable_zcombox()
        assert win.energyComBox.isEnabled()
    finally:
        win.deleteLater()
        app.processEvents()


def test_visualapp_disable_reglinecheck(monkeypatch):
    """Regression checkbox disabled for 3D/Contour/Surface plots."""
    app, win = _make_window(monkeypatch)
    try:
        for plot_type in ('3D Scatter', '2D Contour', '3D Surface'):
            win.plotTypeComBox.setCurrentText(plot_type)
            win.disable_reglinecheck()
            assert not win.regLineChkBtn.isEnabled(), plot_type

        win.plotTypeComBox.setCurrentText('2D Scatter')
        win.disable_reglinecheck()
        assert win.regLineChkBtn.isEnabled()
    finally:
        win.deleteLater()
        app.processEvents()


def test_visualapp_disable_reglinecheck_surface_branch(monkeypatch):
    """Surface-specific branch disables regression lines for 2D Surface."""
    app, win = _make_window(monkeypatch)
    try:
        win.plotTypeComBox.addItem('2D Surface')
        win.plotTypeComBox.setCurrentText('2D Surface')
        win.disable_reglinecheck()
        assert not win.regLineChkBtn.isEnabled()
    finally:
        win.deleteLater()
        app.processEvents()


# ---------------------------------------------------------------------------
# VisualApp — display_message_box
# ---------------------------------------------------------------------------

def test_visualapp_display_message_box(monkeypatch):
    """display_message_box completes without raising on OK click."""
    monkeypatch.setattr(
        QtWidgets.QMessageBox, 'question',
        staticmethod(lambda *args, **kwargs: QtWidgets.QMessageBox.Ok))
    app, win = _make_window(monkeypatch)
    try:
        win.display_message_box('Test title', 'Test body')
    finally:
        win.deleteLater()
        app.processEvents()


def test_visualapp_init_calls_load_file_when_iofile_given(monkeypatch):
    """Providing iofile triggers _load_file during construction."""
    monkeypatch.setattr(VisualApp, 'show', lambda self: None)
    monkeypatch.setattr(VisualApp, 'center_on_screen', lambda self: None)
    monkeypatch.setattr(
        VisualApp,
        '_load_file',
        lambda self: setattr(self, '_load_file_called', True),
    )
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    win = VisualApp(folder='.', iofile='dummy.rst',
                    sub_files={'rst': None, 'traj': None})
    try:
        assert win._load_file_called is True
    finally:
        win.deleteLater()
        app.processEvents()


def test_visualapp_update_canvas_text_and_font(monkeypatch):
    """Text and font helpers update labels, titles and tick sizes."""
    app, win = _make_window(monkeypatch)
    try:
        win.settings = _make_settings()
        win.titleLine.setText('Main title')
        win.xaxisLabel.setText('Order 1')
        win.yaxisLabel.setText('Order 2')
        win.zaxisLabel.setText('Density')

        image = win.myfig.ax.imshow([[0.0, 1.0], [1.0, 0.0]])
        win.myfig.cbar = win.myfig.fig.colorbar(image, cax=win.myfig.cbar_ax)

        win._update_canvas_text()
        assert win.myfig.title.get_text() == 'Main title'
        assert win.myfig.xaxis.get_text() == 'Order 1'
        assert win.myfig.yaxis.get_text() == 'Order 2'
        assert win.myfig.zaxis.get_text() == 'Density'

        win.titleSizeSpin.setValue(18)
        win.axesSizeSpin.setValue(14)
        win._update_canvas_font()
        assert win.settings['title font'] == 18
        assert win.settings['axes font'] == 14
        assert win.myfig.title.get_fontsize() == 18
        assert win.myfig.xaxis.get_fontsize() == 14
        assert win.myfig.yaxis.get_fontsize() == 14
        assert win.myfig.zaxis.get_fontsize() == 14
    finally:
        win.deleteLater()
        app.processEvents()


def test_visualapp_save_methods_with_selected_data(monkeypatch, tmp_path):
    """save_hdf5/save_json/save_png write files when data is selected."""
    app, win = _make_window(monkeypatch)
    recorded = {}
    original_to_hdf = pd.Series.to_hdf
    try:
        win.folder = str(tmp_path)
        win.settings = _make_settings()
        win.dataobject = _DataObjectStub([1.0, 2.0], [3.0, 4.0], [5.0, 6.0])

        def _fake_to_hdf(self, outfile, *args, **kwargs):
            recorded['outfile'] = outfile

        monkeypatch.setattr(pd.Series, 'to_hdf', _fake_to_hdf)
        win.save_hdf5()
        assert recorded['outfile'].endswith('.hdf5')
        assert win.statusbar.currentMessage().endswith('.hdf5')

        win.save_json()
        assert list(tmp_path.glob('*.json'))
        assert win.statusbar.currentMessage().endswith('.json')

        win.save_png()
        assert list(tmp_path.glob('*.png'))
        assert win.statusbar.currentMessage().endswith('.png')
    finally:
        pd.Series.to_hdf = original_to_hdf
        win.deleteLater()
        app.processEvents()


def test_visualapp_save_textdata_without_z(monkeypatch, tmp_path):
    """save_textdata writes two-column output when z data is absent."""
    app, win = _make_window(monkeypatch)
    try:
        win.folder = str(tmp_path)
        win.settings = _make_settings(E='None')
        win.dataobject = _DataObjectStub([1.0, 2.0], [3.0, 4.0], [])

        win.save_textdata()

        text_files = list(tmp_path.glob('*.txt'))
        assert len(text_files) == 1
        content = text_files[0].read_text(encoding='utf8')
        assert 'op1' in content
        assert 'op2' in content
        assert 'potE' not in content
    finally:
        win.deleteLater()
        app.processEvents()


def test_visualapp_save_textdata_with_z(monkeypatch, tmp_path):
    """save_textdata writes three-column output when z data exists."""
    app, win = _make_window(monkeypatch)
    try:
        win.folder = str(tmp_path)
        win.settings = _make_settings(method=['Contour'])
        win.dataobject = _DataObjectStub([1.0, 2.0], [3.0, 4.0], [5.0, 6.0])

        win.save_textdata()

        text_files = list(tmp_path.glob('*.txt'))
        assert len(text_files) == 1
        content = text_files[0].read_text(encoding='utf8')
        assert 'op1' in content
        assert 'op2' in content
        assert 'potE' in content
    finally:
        win.deleteLater()
        app.processEvents()


def test_visualapp_save_script_density_branch(monkeypatch, tmp_path):
    """save_script generates a histogram-based reconstruction script."""
    app, win = _make_window(monkeypatch)
    try:
        win.folder = str(tmp_path)
        win.settings = _make_settings(method=['Density'], E='None')
        win.dataobject = _DataObjectStub([1.0, 2.0], [3.0, 4.0], [])
        monkeypatch.setattr(win, '_save_sim_data_hdf5', lambda: None)
        monkeypatch.setattr(win, '_get_settings', lambda: None)

        win.save_script()

        scripts = list(tmp_path.glob('makefigure_*.py'))
        assert len(scripts) == 1
        content = scripts[0].read_text(encoding='utf8')
        assert 'hist2d' in content
        assert "cbar = fig.colorbar(surf[3], cax=cbar_ax)" in content
        assert 'Density' in scripts[0].name
    finally:
        win.deleteLater()
        app.processEvents()


def test_visualapp_save_script_surface_branch(monkeypatch, tmp_path):
    """save_script includes z-value reconstruction for non-density plots."""
    app, win = _make_window(monkeypatch)
    try:
        (tmp_path / '000').mkdir()
        (tmp_path / '001').mkdir()
        (tmp_path / 'notes').mkdir()
        win.folder = str(tmp_path)
        win.settings = _make_settings(method=['Contour'], fol='All')
        win.dataobject = _DataObjectStub([1.0, 2.0], [3.0, 4.0], [5.0, 6.0])
        monkeypatch.setattr(win, '_save_sim_data_hdf5', lambda: None)
        monkeypatch.setattr(win, '_get_settings', lambda: None)

        win.save_script()

        scripts = list(tmp_path.glob('makefigure_*.py'))
        assert len(scripts) == 1
        content = scripts[0].read_text(encoding='utf8')
        assert 'fol = [' in content
        assert "'000'" in content
        assert "'001'" in content
        assert "zl = 'potE'" in content
        assert 'contourf' in content
    finally:
        win.deleteLater()
        app.processEvents()


def test_visualapp_change_cmap_branches(monkeypatch):
    """_change_cmap handles density, scatter and failure branches."""
    app, win = _make_window(monkeypatch)
    try:
        win.cmapComBox.setEditText('plasma')

        win.settings = _make_settings(method=['Density'])
        win.myfig.surf = win.myfig.ax.hist2d(
            [0.0, 1.0], [1.0, 0.0], bins=(2, 2), cmap='viridis', density=True
        )
        win._change_cmap()
        assert win.myfig.surf[3].get_cmap().name == 'plasma'

        win.settings = _make_settings(method=['Scatter'])
        win._change_cmap()
        assert 're-drawn to update colors' in win.statusbar.currentMessage()

        class _Surface:
            def __init__(self):
                self.cmap = None

            def set_cmap(self, cmap):
                self.cmap = cmap

        surf = _Surface()
        win.settings = _make_settings(method=['Contour'])
        win.myfig.surf = surf
        win._change_cmap()
        assert surf.cmap == 'plasma'

        win.myfig.surf = object()
        win._change_cmap()
        assert 'Could not recognize colormap' in win.statusbar.currentMessage()

        win.settings = _make_settings(method=['Density'])
        win.myfig.surf = object()
        win._change_cmap()
        assert 'Figure not recognized' in win.statusbar.currentMessage()

        win.cmapComBox.setEditText('definitely-not-a-colormap')
        win.settings = _make_settings(method=['Contour'])
        win.myfig.surf = surf
        win._change_cmap()
        msg = win.statusbar.currentMessage()
        assert 'Chosen colormap not recognized' in msg
    finally:
        win.deleteLater()
        app.processEvents()


def test_visualapp_get_op_range_and_change_zoom(monkeypatch):
    """Range/zoom helpers update GUI fields and plotted limits."""
    app, win = _make_window(monkeypatch)
    try:
        win.dataobject = _DataObjectStub(
            [1.0], [2.0], [3.0], infos={'interfaces': [0.25, 0.5, 1.25]}
        )
        win._get_op_range()
        assert win.minDataRange.text() == '0.25'
        assert win.maxDataRange.text() == '1.25'

        image = win.myfig.ax.imshow([[0.0, 1.0], [1.0, 0.0]])
        win.myfig.cbar = win.myfig.fig.colorbar(image, cax=win.myfig.cbar_ax)
        win.settings = _make_settings(dim=2, **{
            'x-limits': (1.0, 2.0),
            'y-limits': (3.0, 4.0),
            'z-limits': (-1.0, 2.0),
        })
        win._change_zoom()
        xlim = tuple(round(v, 1) for v in win.myfig.ax.get_xlim())
        assert xlim == (1.0, 2.0)
        ylim = tuple(round(v, 1) for v in win.myfig.ax.get_ylim())
        assert ylim == (3.0, 4.0)

        win.myfig.set_up(dim=3)
        scalar_map = cm.ScalarMappable(cmap='viridis')
        scalar_map.set_clim(0.0, 1.0)
        cbar = win.myfig.fig.colorbar(scalar_map, cax=win.myfig.cbar_ax)
        win.myfig.cbar = cbar
        win.settings = _make_settings(dim=3, **{
            'x-limits': (0.0, 1.0),
            'y-limits': (2.0, 3.0),
            'z-limits': (4.0, 5.0),
        })
        win._change_zoom()
        xlim3d = tuple(round(v, 1) for v in win.myfig.ax.get_xlim3d())
        assert xlim3d == (0.0, 1.0)
        ylim3d = tuple(round(v, 1) for v in win.myfig.ax.get_ylim3d())
        assert ylim3d == (2.0, 3.0)
        zlim3d = tuple(round(v, 1) for v in win.myfig.ax.get_zlim3d())
        assert zlim3d == (4.0, 5.0)
    finally:
        win.deleteLater()
        app.processEvents()


def test_visualapp_update_preview_runs_refresh_steps(monkeypatch):
    """update_preview refreshes settings, text, zoom and colormap in order."""
    app, win = _make_window(monkeypatch)
    steps = []
    try:
        win.settings = _make_settings()
        win.myfig.surf = object()
        monkeypatch.setattr(
            win, '_get_settings', lambda: steps.append('settings'))
        monkeypatch.setattr(
            win, '_update_canvas_text', lambda: steps.append('text'))
        monkeypatch.setattr(
            win, '_change_zoom', lambda: steps.append('zoom'))
        monkeypatch.setattr(
            win, '_change_cmap', lambda: steps.append('cmap'))

        win.update_preview()

        assert steps == ['settings', 'text', 'zoom', 'cmap']
        assert win.statusbar.currentMessage() == 'Plot ready!'
    finally:
        win.deleteLater()
        app.processEvents()


def test_visualapp_toggle_helpers(monkeypatch):
    """Toggle helpers update interface lines, regression lines and titles."""
    app, win = _make_window(monkeypatch)
    emit_called = []
    try:
        interface_a = win.myfig.ax.axvline(0.1)
        interface_b = win.myfig.ax.axvline(0.2)
        win.myfig.intf = [interface_a, interface_b]
        win.intShowChkBtn.setChecked(False)
        win.toggle_intf()
        assert not interface_a.get_visible()
        assert not interface_b.get_visible()

        regline = win.myfig.ax.plot([0.0, 1.0], [1.0, 0.0])[0]
        win.settings = _make_settings(dim=2)
        win.myfig.regl = [regline]
        win.regLineChkBtn.setChecked(False)
        win.toggle_regl()
        assert not regline.get_visible()

        win.settings = _make_settings()
        win.titleLine.setText('Visible title')
        win.xaxisLabel.setText('x')
        win.yaxisLabel.setText('y')
        win.zaxisLabel.setText('z')
        win._update_canvas_text()
        win.showTitleChkBtn.setChecked(False)
        win.toggle_titles()
        assert not win.myfig.title.get_visible()
        assert not win.myfig.xaxis.get_visible()
        assert not win.myfig.yaxis.get_visible()
        assert not win.myfig.zaxis.get_visible()

        win.dataobject = object()
        monkeypatch.setattr(
            win, 'emit_settings', lambda: emit_called.append(True))
        win.toggle_only_stored_trjs()
        assert emit_called == [True]
    finally:
        win.deleteLater()
        app.processEvents()


def test_customfigcanvas_setup_without_colorbar_margin():
    """set_up(cbar=False) uses the wider subplot layout branch."""
    from pyretis.pyvisa.visualize import CustomFigCanvas

    canvas = CustomFigCanvas()
    canvas.set_up(cbar=False)
    right = round(canvas.fig.subplotpars.right, 2)
    assert right == 0.90

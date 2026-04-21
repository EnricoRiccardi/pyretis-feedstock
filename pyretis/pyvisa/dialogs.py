# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Dialog classes for the PyVisA GUI.

Contains all QDialog subclasses used by VisualApp for loading data.

Important classes defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

LoadHdf5Dialog (:py:class:`.LoadHdf5Dialog`)
    File picker for HDF5 / zip compressed simulation data.

LoadOrderParamDialog (:py:class:`.LoadOrderParamDialog`)
    File picker for a standalone order parameter text file.

RecalculateDialog (:py:class:`.RecalculateDialog`)
    Dialog to select either a PyRETIS .rst input file or a Python script
    for order parameter recalculation from snapshot trajectories.
"""
import os
from PyQt5 import QtWidgets


class LoadHdf5Dialog(QtWidgets.QDialog):  # pragma: no cover
    """Dialog for selecting an HDF5 or zip compressed data file."""

    def __init__(self, parent=None):
        """Initialise the dialog and build the UI."""
        super().__init__(parent)
        self.setWindowTitle('Load Data (HDF5 / zip)')
        self.setMinimumWidth(520)
        self._path = None
        self._setup_ui()

    def _setup_ui(self):
        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(QtWidgets.QLabel(
            'Select a PyVisA compressed data file (.hdf5 or .hdf5.zip).'
        ))
        row = QtWidgets.QHBoxLayout()
        self._line = QtWidgets.QLineEdit()
        self._line.setReadOnly(True)
        self._line.setPlaceholderText('No file selected…')
        btn = QtWidgets.QPushButton('Browse…')
        btn.clicked.connect(self._browse)
        row.addWidget(self._line)
        row.addWidget(btn)
        layout.addLayout(row)
        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel
        )
        buttons.accepted.connect(self._on_accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _browse(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Select HDF5 / zip file',
            filter='Compressed data (*.hdf5 *.zip);;All files (*.*)',
        )
        if path:
            self._path = path
            self._line.setText(path)

    def _on_accept(self):
        if not self._path:
            QtWidgets.QMessageBox.warning(
                self, 'No file selected',
                'Please select an HDF5 or zip file.')
            return
        self.accept()

    def get_path(self):
        """Return the selected file path."""
        return self._path


class LoadOrderParamDialog(QtWidgets.QDialog):  # pragma: no cover
    """Dialog for selecting a standalone order parameter file.

    The file is expected to follow the PyRETIS order.txt format::

        Recalculated data
        #     Time       Orderp
                 0     2.819435
                 1     2.184310
                 ...

    In this mode ensemble-specific controls are hidden and all PyVisA
    visualizations operate on the single time series.
    """

    def __init__(self, parent=None):
        """Initialise the dialog and build the UI."""
        super().__init__(parent)
        self.setWindowTitle('Load Order Parameter File')
        self.setMinimumWidth(560)
        self._path = None
        self._setup_ui()

    def _setup_ui(self):
        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(QtWidgets.QLabel(
            'Select a standalone order parameter file (e.g. order.txt).\n\n'
            'Expected format:\n'
            '  Recalculated data\n'
            '  #     Time       Orderp\n'
            '           0     2.819435\n'
            '           1     2.184310\n'
            '  ...\n\n'
            'Ensemble-specific controls will be hidden in this mode.\n'
            'All PyVisA visualisations and analyses remain available.'
        ))
        row = QtWidgets.QHBoxLayout()
        self._line = QtWidgets.QLineEdit()
        self._line.setReadOnly(True)
        self._line.setPlaceholderText('No file selected…')
        btn = QtWidgets.QPushButton('Browse…')
        btn.clicked.connect(self._browse)
        row.addWidget(self._line)
        row.addWidget(btn)
        layout.addLayout(row)
        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel
        )
        buttons.accepted.connect(self._on_accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _browse(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Select order parameter file',
            filter='Text files (*.txt);;All files (*.*)',
        )
        if path:
            self._path = path
            self._line.setText(path)

    def _on_accept(self):
        if not self._path:
            QtWidgets.QMessageBox.warning(
                self, 'No file selected',
                'Please select an order parameter file.')
            return
        self.accept()

    def get_path(self):
        """Return the selected file path."""
        return self._path


class RecalculateDialog(QtWidgets.QDialog):  # pragma: no cover
    """Dialog for recalculating order parameters from snapshot trajectories.

    Presents two mutually exclusive options:

    A – PyRETIS input file (.rst)
        The order parameter is recalculated using the definition stored in
        the simulation settings file.  PyVisA scans the directory tree for
        trajectory files, rewrites the order.txt files, and loads the result.

    B – Python script (.py)
        A user-supplied script is executed.  It must write order parameter
        data to *stdout* as a JSON list of lists (one inner list per frame,
        each element an order parameter value) or as a CSV table.
        PyVisA captures the output, writes an order.txt next to the script,
        and loads the result as a single order parameter file.
    """

    def __init__(self, parent=None):
        """Initialise the dialog and build the UI."""
        super().__init__(parent)
        self.setWindowTitle('Recalculate from Snapshot / Trajectory')
        self.setMinimumWidth(620)
        self._mode = None
        self._path = None
        self._setup_ui()

    def _setup_ui(self):
        layout = QtWidgets.QVBoxLayout(self)

        # --- Option A ---
        self._radio_rst = QtWidgets.QRadioButton(
            'A  –  PyRETIS input file (.rst)\n'
            '       Recalculate the order parameter using the PyRETIS\n'
            '       simulation settings. PyVisA will scan the directory\n'
            '       for trajectory files and rewrite order.txt.'
        )
        self._radio_rst.setChecked(True)
        layout.addWidget(self._radio_rst)

        row_rst = QtWidgets.QHBoxLayout()
        self._line_rst = QtWidgets.QLineEdit()
        self._line_rst.setReadOnly(True)
        self._line_rst.setPlaceholderText('Select a .rst file…')
        btn_rst = QtWidgets.QPushButton('Browse…')
        btn_rst.clicked.connect(self._browse_rst)
        row_rst.addWidget(self._line_rst)
        row_rst.addWidget(btn_rst)
        layout.addLayout(row_rst)

        layout.addSpacing(14)

        # --- Option B ---
        self._radio_script = QtWidgets.QRadioButton(
            'B  –  Python script (.py)\n'
            '       Run a user-supplied script that prints order parameter\n'
            '       data to stdout as a JSON list of lists  [[op1, …], …]\n'
            '       or as a CSV table (one row per frame).\n'
            '       The result is saved as order.txt next to the script.'
        )
        layout.addWidget(self._radio_script)

        row_script = QtWidgets.QHBoxLayout()
        self._line_script = QtWidgets.QLineEdit()
        self._line_script.setReadOnly(True)
        self._line_script.setPlaceholderText('Select a .py script…')
        btn_script = QtWidgets.QPushButton('Browse…')
        btn_script.clicked.connect(self._browse_script)
        row_script.addWidget(self._line_script)
        row_script.addWidget(btn_script)
        layout.addLayout(row_script)

        layout.addStretch()
        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel
        )
        buttons.accepted.connect(self._on_accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _browse_rst(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Select PyRETIS input file',
            filter='PyRETIS input (*.rst);;All files (*.*)',
        )
        if path:
            self._line_rst.setText(path)
            self._radio_rst.setChecked(True)

    def _browse_script(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Select Python script',
            filter='Python scripts (*.py);;All files (*.*)',
        )
        if path:
            self._line_script.setText(path)
            self._radio_script.setChecked(True)

    def _on_accept(self):
        if self._radio_rst.isChecked():
            path = self._line_rst.text().strip()
            if not path or not os.path.isfile(path):
                QtWidgets.QMessageBox.warning(
                    self, 'No file', 'Please select a valid .rst file.')
                return
            self._mode = 'rst'
            self._path = path
        else:
            path = self._line_script.text().strip()
            if not path or not os.path.isfile(path):
                QtWidgets.QMessageBox.warning(
                    self, 'No file', 'Please select a valid .py script.')
                return
            self._mode = 'script'
            self._path = path
        self.accept()

    def get_selection(self):
        """Return (mode, filepath) where mode is 'rst' or 'script'."""
        return self._mode, self._path

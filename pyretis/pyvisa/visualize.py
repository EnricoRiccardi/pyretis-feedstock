# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""GUI application for visualizing simulation data.

This is a PyQt5 file, using a custom made ui layout, with a window which
displays different plots created by the path density tool. Window can either
display data of a pre-compiled compressed file generated with the
orderparam_density.py module, or compile on in/out files by importing
orderparam_density and executing before displaying the results.


Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

visualize_main (:py:func: `.visualize_main`)
    Method to load the GUI.


Important classes defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

VisualApp (:py:class:`.VisualApp`)
    The PyQt5 GUI class that handles all the functionalities.

CustomFigCanvas (:py:class:`.VisualApp`)
    Class to handle the MatPlotLiv canvas.

DataSlave (:py:class:`.VisualApp`)
    The class to handle the data for the plots.

VisualObject (:py:class:`.VisualApp`)
    The class that handles the load from files.

DataObject (:py:class:`.VisualApp`)
    The class that loads the PyRETIS files.

"""
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
import os
import warnings
import subprocess
import sys
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas
)
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import mdtraj as md
import numpy as np
import pandas as pd
from PyQt5 import QtCore, QtGui
from PyQt5 import QtWidgets, uic
from scipy.spatial import distance  # pylint: disable=import-error
from pyretis.inout import settings
from pyretis.inout.common import TRJ_FORMATS
from pyretis.setup.common import create_orderparameter
from pyretis.tools.recalculate_order import recalculate_order
from pyretis.pyvisa.common import (read_traj_txt_file,
                                   recalculate_all,
                                   try_data_shift,
                                   find_rst_file,
                                   get_cv_names)
from pyretis.pyvisa.orderparam_density import (PathDensity,
                                               PathVisualize,
                                               Trajectory,
                                               pyvisa_zip,
                                               remove_nan)
from pyretis.pyvisa.plotting import (gen_surface,
                                     plot_regline,
                                     plot_int_plane)
from pyretis.pyvisa.statistical_methods import (pyvisa_pca,
                                                k_means,
                                                spectral,
                                                gaussian_mixture,
                                                hierarchical,
                                                correlation_matrix,
                                                random_forest,
                                                decision_tree)

warnings.filterwarnings('ignore', category=pd.errors.PerformanceWarning)

# Hard-coded labels for energies and time/cycle steps
ENERGYLABELS = ['time', 'cycle', 'potE', 'kinE', 'totE']


def _set_tick_fontsize(axis, fontsize):
    """Set major tick label font sizes across Matplotlib versions."""
    for tick in axis.get_major_ticks():
        labels = [getattr(tick, 'label1', None),
                  getattr(tick, 'label2', None)]
        if not any(labels):
            labels = [getattr(tick, 'label', None)]
        for label in labels:
            if label is not None:
                label.set_fontsize(fontsize)


def _redraw_colorbar(colorbar):
    """Refresh a Matplotlib colorbar across Matplotlib versions."""
    if hasattr(colorbar, 'draw_all'):
        colorbar.draw_all()
    else:
        colorbar.update_normal(colorbar.mappable)


def _remove_last_line(ax):
    """Remove the latest line from an axes across Matplotlib versions."""
    if ax.lines:
        ax.lines[-1].remove()


UI_VW, QT_BC = uic.loadUiType(
    os.path.dirname(os.path.realpath(__file__)) +
    '/' + 'pyvisaGUI.ui')


class VisualApp(QtWidgets.QMainWindow, UI_VW):
    """Class definition of the path visualization window GUI for PyRETIS.

    Application opens a QMainWindow object with a QFrame for inserting
    a matplotlib figure canvas into, and a QDropWidget with the relevant
    settings for plotting.

    Attributes
    ----------
    mainworker : A QThread object that stores (and executes) the path density
        data back-end, initialized by the VisualApp window. When
        called will update the data lists and return to VisualApp
        for plotting.

    Functions
    ---------
    close_event, action_reload, center_on_screen, toggle_buttons :
        Functions that relate to closing of window, reloading data,
        centering the window on screen and disable key buttons while
        pdating data/figure.
    _update_canvas_{title/text} :
        Functions that update the title and label text of plot settings
        and in the figure.
    _get_settings : Gets plot settings from VisualApp dropwidget.
    _get_savename : Generate generic savename based on current settings.
    save_png : Function that store the current figure shown in VisualApp
        window as a .png image, with name generated from data settings.
    save_hdf5 : Function that store the current figure shown in VisualApp
        window a hdf5 save, with name generated from data settings.
    save_json : Function that store the current figure shown in VisualApp
        window a json save, with name generated from data settings.
    save_script : Function that generates a script to reconstruct the current
        plot from the compressed file.
    _load_file : Opens QFileDialog for choosing file to load in VisualApp.
    _load_data : Loads data from hdf5 save.
    _load_data_output : Executes a PathDensity on input/output file before
        displaying results.
    _reload : Clears old data and initializes loading new data from file.
    _change_{cmap/zoom} :
        Updates figure with different colormap, and with xlims+ylims
    toggle_{intf/regl/cbar/titles} :
        Functions that shows/hides interfaces, regression line,
        colorbar, and titles/labes, respectively.
    emit_settings : Function that is called when pressing a QPushButton labeled
        'Update' after choosing desired settings for plot.
        Sends settings to mainworker, which updates the data of
        x, y (and z) lists before sending back to VisualApp for
        plotting.
    on_pick : Function that is called when clicking on the plot.
        Writes the x, y coordinate, cycle, ensemble name,
        status and mc move of the closest corresponding point to
        wherever the user clicks on the plot and fills in the slots
        in the Analysis tab in PyVisA and plots the closest point.
    draw_trajectory : Function that is called when pressing a QPushButton.
        When pressing the button labeled "show trj" in the "Analysis"
        tab in PyVisA. The function plots all points in the trajectory
        corresponding to the closest point defined by the
        users pick on the plot.
    play_cycle : Function that is called when pressing a QPushButton.
        Pressing the QPushButton labeled "play" button in the "Analysis"
        tab in PyVisA causes the function to check if the corresponding
        trajectory files exists, and reverses them through rev_trj if
        necessary before merging them to a pdb file, and loading them
        to PyMOL. NB: Requires the user to have PyMOL installed in
        order to work.
    create_full_traj: Function that creates a pdb-file for a trajectory.
    display_message_box : Function that displays a QMessageBox.
        Function displays a message box in use for error handling and
        messaging the user about operations and choices that are not
        availeable.
    toggle_only_stored_trjs: Only plot data with trajectory files.
        Function that loops over all trajectories in the chosen ensembles
        and cycles and only plots those trajectories which have
        existing trajectory files.
    recalculate : Recalculate OP and CV's in new order.txt files.
        Function activating from the QAction labeled
        "Recalculate data" in the file menu of PyVisA. Uses the tools
        available in the recalculate_order.py file to load the
        trajectory files of the simulation, and calculate new
        collective variables which need to be declared in the
        input files "retis.rst" and "orderp.py". The new data will be
        stored in the full order.txt in the ensemble folders, and
        also the order.txt in the cycle folders of each ensemble
        enabling PyVisA to load and visualize the new descriptors.
        The previous data files are renamed as order_old.
    statistic_analysis: Perform statistical analysis on plotted data.
        A function that is connected to the statistical function in the
        statistical_methods file, and that takes the users choice of analysis
        and run the said choice.
    show_correlation : Show the correlation between simulation variables.
        Shows the Pearson correlation in a matrix form for all the variables
        from the simulation.
    refresh_data: Signal a walk of the data.
        Performs a new walk of the data, to use during a running simulation
        to visualise its development.

    Signals
    -------
    start_cmd : PyQt Signal
        It tells worker to start using string, cmd.
    send_settings : PyQt Signal
        It sends settings for data_get to worker thread.

    Slots
    -----
    update_cycle : function bound to a pyqtSignal sent by the mainworker
        thread. Signal sets upper and lower bounds for max/min cycles used
        in plotting.
    update_fig : PyQt Slot
        Pyqt slot, worker thread calls function by sending 3 lists
        of floats for plotting in this function.

    """

    send_settings = QtCore.pyqtSignal(dict, name='Data settings')
    start_cmd = QtCore.pyqtSignal(str)

    def __init__(self, folder, iofile, sub_files):
        """Initialize the VisualApp QMainWindow.

        Sets up ui from file, centers on screen and sets up mainworker QThread
        for data loading.

        Parameters
        ----------
        folder : string
            The path of the target folder.
        iofile : string
            The name of the main input file.
        sub_files : dict
             Contains the name of multiple input files.

        """
        QtWidgets.QMainWindow.__init__(self)
        UI_VW.__init__(self)
        # Delete memory on close?
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setupUi(self)
        self.center_on_screen()
        self.setWindowTitle('PyVisA - Visualization and Analysis')
        self.folder = folder
        self.iofile = iofile

        # Link to the input file if PyVisA was loaded with a rst file
        self.rst_file = sub_files['rst']
        # Link to the trajectory file if -data was used
        self.trajfile = sub_files['traj']
        self.settings, self.dataobject, self.thread = None, None, None
        # Load data in separate thread
        if iofile is not None:
            self._load_file()
        self._init_widget()
        self.myfig.fig.canvas.mpl_connect(
            'button_press_event', self.on_pick)

        self.xmin, self.xmax = None, None
        self.ymin, self.ymax = None, None
        self.zmin, self.zmax = None, None

    def close_event(self):  # pragma: no cover
        """Event function interpreter to make the event to a boolean."""
        self.close()

    def closeEvent(self, event):  # pylint: disable=C0103  # pragma: no cover
        """Event function, activated when attempting to close VisualApp.

        Will create a QMessage prompt asking for confirmation of exit by
        user. If confirmed, will stop Qthread of DataObject and attempt to
        clear dataobject from memory.

        Parameters
        ----------
        event : QtWidget event
            It activates the closing of the VisualApp.

        """
        reply = QtWidgets.QMessageBox.question(
            self, 'Quit?', 'Are you sure you want to quit?',
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
            QtWidgets.QMessageBox.No)

        if reply == QtWidgets.QMessageBox.Yes:
            self.statusbar.showMessage('Closing window')
            try:
                self.thread.quit()
                self.thread.wait()
                del self.dataobject
                del self.thread
            except AttributeError:
                pass
            finally:
                event.accept()
        else:
            event.ignore()

    def action_reload(self):  # pragma: no cover
        """Display a QMessagebox to confirm reload action.

        Function that displays a QMessagebox for user, double confirming
        the reload data action.
        """
        if not self.dataobject:
            self._reload()
        else:
            load = QtWidgets.QMessageBox.question(
                self,
                "Reload data from file",
                "Are you sure to discard the current data and "
                "reload a new set ?",
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
            if load == QtWidgets.QMessageBox.Yes:
                self.statusbar.showMessage('Deleting old data')
                try:
                    self.thread.quit()
                    self.thread.wait()
                except AttributeError:
                    # No self.thread running.
                    pass
                self._reload()

    def center_on_screen(self):  # pragma: no cover
        """Centers widget/window on screen."""
        resolution = QtWidgets.QDesktopWidget().screenGeometry()
        self.move(int((resolution.width() - self.frameSize().width()) / 2),
                  int((resolution.height() - self.frameSize().height()) / 2))

    def toggle_buttons(self, onoff):
        """Toggles enable-state of key buttons of VisualApp.

        It avoids user-created conflicts while busy updating data and drawing
        figure.

        Parameters
        ----------
        onoff : boolean
            Determines the status of the button.

        """
        self.updateBtn.setEnabled(onoff)
        self.previewBtn.setEnabled(onoff)
        self.intShowChkBtn.setEnabled(onoff)
        self.regLineChkBtn.setEnabled(onoff)
        self.showTitleChkBtn.setEnabled(onoff)
        self.saveFigBtn.setEnabled(onoff)
        self.previewBtn.setEnabled(onoff)
        self.refreshBtn.setEnabled(onoff)

    def disable_zcombox(self):
        """Disabling combobox with z-values if plotting density."""
        if 'Density' in self.plotTypeComBox.currentText():
            self.energyComBox.setEnabled(False)
        else:
            self.energyComBox.setEnabled(True)

    def disable_reglinecheck(self):
        """Disabling regression line checkbox if 3D method."""
        if '3D' in self.plotTypeComBox.currentText():
            self.regLineChkBtn.setEnabled(False)
        elif 'Contour' in self.plotTypeComBox.currentText():
            self.regLineChkBtn.setEnabled(False)
        elif 'Surface' in self.plotTypeComBox.currentText():
            self.regLineChkBtn.setEnabled(False)
        else:
            self.regLineChkBtn.setEnabled(True)

    def _update_canvas_font(self):
        """Change the fontsize of title/labels of figure.

        Bound to event only.
        """
        if self.settings is None or any(
            item is None for item in (
                self.myfig.title,
                self.myfig.xaxis,
                self.myfig.yaxis,
                self.myfig.zaxis,
            )
        ):
            self.statusbar.showMessage('Figure not ready!')
            return
        titlefont = self.titleSizeSpin.value()
        self.settings['title font'] = titlefont
        axesfont = self.axesSizeSpin.value()
        self.settings['axes font'] = axesfont
        self.myfig.title.set_fontsize(titlefont)
        self.myfig.xaxis.set_fontsize(axesfont)
        self.myfig.yaxis.set_fontsize(axesfont)
        self.myfig.zaxis.set_fontsize(axesfont)
        if self.myfig.cbar is not None:
            self.myfig.cbar.ax.tick_params(labelsize=axesfont)
        _set_tick_fontsize(self.myfig.ax.xaxis, axesfont)
        _set_tick_fontsize(self.myfig.ax.yaxis, axesfont)
        # Default pyplot tick size: lenght=3.5, width=1.0
        tick_l = round(axesfont / 4, 2)
        tick_w = round(axesfont / 10, 2)
        self.myfig.ax.tick_params(length=tick_l, width=tick_w)
        self.myfig.fig.canvas.draw()

    def _update_canvas_text(self):
        """Create the titles and labels of figure.

        This is in VisualApp from QLineEdits current text.
        Called by update_fig() and bound to event.
        """
        if self.settings is None:
            self.statusbar.showMessage('Figure not ready!')
            return
        self.settings['show titles'] = self.showTitleChkBtn.isChecked()
        show = self.settings['show titles']
        self.myfig.title = self.myfig.ax.set_title(
            self.titleLine.text(),
            fontsize=self.settings['title font'],
            visible=show)
        self.myfig.xaxis = self.myfig.ax.set_xlabel(
            self.xaxisLabel.text(),
            fontsize=self.settings['axes font'],
            visible=show)
        self.myfig.yaxis = self.myfig.ax.set_ylabel(
            self.yaxisLabel.text(),
            fontsize=self.settings['axes font'],
            visible=show)
        self.myfig.zaxis = self.myfig.cbar_ax.set_title(
            self.zaxisLabel.text(),
            fontsize=self.settings['axes font'],
            visible=show)
        self.myfig.fig.canvas.draw()

    def _get_titles(self):  # pragma: no cover
        """Generate tiles and labels and update QLineEdits of VisualApp.

        Function that generates generic titles and labels for the figure
        axes, colorbar, etc. and updates QLineEdits of VisualApp with
        those titles/labels.
        """
        method = self.settings['method']
        name = method[0]
        acc = self.settings['ACC']
        shift = self.settings['try_shift']
        fol = self.settings['fol']
        op1 = self.settings['op1']
        op2 = self.settings['op2']
        energy = self.settings['E']

        xyz = op1 + '-' + op2
        if 'Density' not in method:
            xyz += '-' + energy
            self.myfig.cbar_ax.set_title(energy)
        if self.settings['weight']:
            name += '(W)'
        cmin, cmax = self.settings['min_cycle'], self.settings['max_cycle']
        title = (f'{name}, using {fol} folder(s), with {acc} paths, '
                 f'cycles: {cmin}-{cmax}')
        if shift:
            title += '\nData shifted'
        # Setting text in QLineEdits of VisualApp
        self.xaxisLabel.setText(op1)
        self.yaxisLabel.setText(op2)
        if self.settings['method'][0] == 'Density':
            self.zaxisLabel.setText('Density')
        else:
            self.zaxisLabel.setText(energy)
        self.titleLine.setText(title)
        # Updating text in plot from text in QLineEdits
        self._update_canvas_text()

    def _get_savename(self):  # pragma: no cover
        """Generate generic save-name.

        Function that generates a generic save-name for the current
        figure displayed in VisualApp based on settings.
        """
        method = self.settings['method']
        method_name = method[0]
        acc = self.settings['ACC']
        fol = self.settings['fol']
        op1 = self.settings['op1']
        op2 = self.settings['op2']
        xyz = op1 + '-' + op2
        cmin, cmax = self.settings['min_cycle'], self.settings['max_cycle']
        if 'Density' not in method and self.settings['E'] != 'None':
            xyz += '-' + self.settings['E']
        elif self.settings['weight']:
            method_name += '(W)'
        name = f'{method_name}_{xyz}_ACC={acc}_folder={fol}_({cmin}-{cmax})'
        return name

    @staticmethod
    def _iffloat(value, default_v):
        """Check if value is a float, and returns it."""
        try:
            return float(value)
        except ValueError:
            return default_v

    def _set_def_limits(self, x, y, z):  # pragma: no cover
        """Update x/y/z limits in GUI with default values.

        Parameters
        ----------
        x, y, z : list
            Floats, the coordinates used in plotting.

        """
        self.statusbar.showMessage('Find limits of data selected.')

        self.xmin = min(x, default=0)
        self.ymin = min(y, default=0)
        self.zmin = min(z, default=0)
        self.xmax = max(x, default=0)
        self.ymax = max(y, default=0)
        self.zmax = max(z, default=0)

        self.xminLine.setText(str(self.xmin))
        self.xmaxLine.setText(str(self.xmax))
        self.yminLine.setText(str(self.ymin))
        self.ymaxLine.setText(str(self.ymax))
        self.zminLine.setText(str(self.zmin))
        self.zmaxLine.setText(str(self.zmax))

    def _set_cur_limits(self):  # pragma: no cover
        """Update x/y/z limits in GUI with current values."""
        self.xminLine.setText(str(self.settings['x-limits'][0]))
        self.yminLine.setText(str(self.settings['y-limits'][0]))
        self.zminLine.setText(str(self.settings['z-limits'][0]))
        self.xmaxLine.setText(str(self.settings['x-limits'][1]))
        self.ymaxLine.setText(str(self.settings['y-limits'][1]))
        self.zmaxLine.setText(str(self.settings['z-limits'][1]))

    def _get_limits(self):  # pragma: no cover
        """Update x/y/z limits from GUI.

        Function that updates the min/max limits of data for plotting
        in the GUI. If invalid inputs, will return default values.
        """
        try:
            float(self.xmin)
        except (AttributeError, TypeError):
            # Not yet defined data, skip
            return

        # Check min/max of x
        nxi = self._iffloat(self.xminLine.text(), self.xmin)
        nxa = self._iffloat(self.xmaxLine.text(), self.xmax)
        if nxi >= nxa:
            self.settings['x-limits'] = (self.xmin, self.xmax)
        else:
            self.settings['x-limits'] = (nxi, nxa)

        # Check min/max of y
        nyi = self._iffloat(self.yminLine.text(), self.ymin)
        nya = self._iffloat(self.ymaxLine.text(), self.ymax)
        if nyi >= nya:
            self.settings['y-limits'] = (self.ymin, self.ymax)
        else:
            self.settings['y-limits'] = (nyi, nya)

        # Check min/max of z
        nzi = self._iffloat(self.zminLine.text(), self.zmin)
        nza = self._iffloat(self.zmaxLine.text(), self.zmax)
        if nzi >= nza:
            self.settings['z-limits'] = (self.zmin, self.zmax)
        else:
            self.settings['z-limits'] = (nzi, nza)
        self._set_cur_limits()

    def _get_settings(self):  # pragma: no cover
        """Update self.settings{} with the current option from the GUI.

        Function that updates the settings of dictionary
        self.settings{} with current options from GUI window.
        """
        full_method = self.plotTypeComBox.currentText()
        dim = int(full_method[0])
        method = full_method.split(' ')[1:]
        weight = bool(method[-1] == '(weight)')
        self.settings = {
            'op1': self.firstOpComBox.currentText(),
            'op2': self.secondOpComBox.currentText(),
            'E': self.energyComBox.currentText(),
            'ACC': self.accComBox.currentText(),
            'fol': self.folderComBox.currentText(),
            'min_cycle': self.minCyclSpin.value(),
            'max_cycle': self.maxCyclSpin.value(),
            'weight': weight,
            'method': method,
            'dim': dim,
            'try_shift': self.dataShiftChkBtn.isChecked(),
            'show_int': self.intShowChkBtn.isChecked(),
            'reg_line': self.regLineChkBtn.isChecked(),
            'res': self.resSpinBox.value(),
            'colormap': self.cmapComBox.currentText(),
            'colorbar': self.cbarChkBtn.isChecked(),
            'title font': self.titleSizeSpin.value(),
            'axes font': self.axesSizeSpin.value(),
            'show titles': self.showTitleChkBtn.isChecked(),
            'stored': self.storedTrjsChkBtn.isChecked(),
            'MC-move': self.mcComBox.currentText(),
            'start_end': self.start_end_box.currentText()
        }
        # Get limits of x/y/z data:
        self._get_limits()

    def save_hdf5(self):
        """Save in hdf5 format the current figure.

        Function that saves the current pyplot.figure object of
        VisualApp and CustomFigCanvas to a .hdf5 file. Function creates
        title and name for file by the current settings for the plot.
        """
        if self.settings is None:
            self.statusbar.showMessage('No data selected')
            return
        x, y, z = self.dataobject.load_all()
        dataframe = pd.Series({'x': x, 'y': y, 'z': z,
                               'settings': self.settings})
        name = self._get_savename() + '.hdf5'
        outfile = os.path.join(self.folder, name)
        dataframe.to_hdf(outfile, key='hdf5', mode='w')
        self.statusbar.showMessage(f'Figure saved as {outfile}')

    def save_json(self):
        """Save in json format the current figure.

        Function that saves the current pyplot.figure object of
        VisualApp and CustomFigCanvas to a .json file. Function creates
        title and name for file by the current settings for the plot.
        """
        if self.settings is None:
            self.statusbar.showMessage('No data selected')
            return
        name = self._get_savename() + '.json'
        outfile = os.path.join(self.folder, name)
        x, y, z = self.dataobject.load_all()
        dataframe = pd.Series({'x': x, 'y': y, 'z': z,
                               'settings': self.settings})
        dataframe.to_json(outfile)
        self.statusbar.showMessage(f'Figure saved as {outfile}')

    def save_png(self):
        """Save current figure as png.

        Function that saves the current plotted figure of VisualApp
        and CustomFigCanvas as a .png file.
        Function creates title and name for file by the current settings
        for the plot.
        """
        if self.settings is None:
            self.statusbar.showMessage('No data selected')
            return
        name = self._get_savename() + '.png'
        outfile = os.path.join(self.folder, name)
        self.myfig.fig.savefig(outfile)
        self.statusbar.showMessage(f'Figure saved: {outfile}')

    def save_textdata(self):
        """Save x, y, z data to file.

        Function that saves the current x, y, z data of the visualizer to
        a text file for use in other operations, plotting, etc.
        """
        if self.settings is None:
            self.statusbar.showMessage('No data selected')
            return
        self.statusbar.showMessage('Saving figure data to file...')
        name = self._get_savename() + '.txt'
        outfile = os.path.join(self.folder, name)
        x, y, z = self.dataobject.load_all()

        if not z:
            line = '\t{}\t\t{}\n'
            header = line.format(self.settings['op1'],
                                 self.settings['op2'])
        else:
            line = '\t{}\t\t{}\t\t{}\n'
            header = line.format(self.settings['op1'],
                                 self.settings['op2'],
                                 self.settings['E'])

        self.statusbar.showMessage(f'Writing file: {outfile}')
        with open(os.path.join(self.folder, outfile),
                  'w', encoding="utf8") as output:
            output.write(header)
            if not z:
                for i, j in zip(x, y):
                    output.write(line.format(i, j))
            else:
                for i, j, k in zip(x, y, z):
                    output.write(line.format(i, j, k))
        self.statusbar.showMessage(f'Figure data saved: {outfile}!')

    def save_script(self):
        """Make makefigure.py script.

        Function that generate and store the current plot from the GUI to
        a separate script file, that can re-generate a fairly similar
        figure to the GUI using stored settings.
        """
        if self.settings is None:
            self.statusbar.showMessage('No data selected')
            return
        datafile = "pyvisa_compressed_data.hdf5"
        self._save_sim_data_hdf5()
        scriptfileformat = 'makefigure_{}_{}_{}_{}.py'
        self._get_settings()
        # Text writes to makefigure.py
        txt = "# Makefigure script\n"
        txt += "import pandas as pd\n"
        txt += "import numpy as np\n"
        txt += "import os.path\n"
        txt += "import matplotlib.pyplot as plt\n"
        txt += "from scipy.interpolate import griddata as scgriddata\n"
        txt += "\n"
        txt += f"datafile = '{datafile}'\n"
        txt += "\n"
        txt += "traj = pd.read_hdf(datafile, key='data')\n"
        txt += "infos = pd.read_hdf(datafile, key='infos')\n"
        txt += "data = [traj, infos]\n"
        txt += "\n"
        txt += "# Dictionary with settings for data load and plotting:\n"
        txt += f"settings = {self.settings}\n"
        txt += "\n"

        if self.settings['fol'] == 'All':
            fol = []
            for folder in next(os.walk(self.folder))[1]:
                if folder[0] == '0':
                    fol.append(folder)
        else:
            fol = [self.settings['fol']]
        txt += f"fol = {fol}\n"
        x, y, z = self.dataobject.load_all()
        op1 = self.settings['op1']
        op2 = self.settings['op2']
        energy = self.settings['E']
        txt += f"xl = '{op1}'\n"
        txt += f"yl = '{op2}'\n"
        if self.settings['method'][0] == 'Density':
            index_data = 0
        else:
            index_data = 1
            txt += f"zl = '{energy}'\n"
        txt += "\n"
        txt += "fig = plt.figure()\n"
        txt += "ax = fig.add_subplot(111)\n"
        txt += "fig.subplots_adjust(left=0.1, right=0.85, "
        txt += "bottom=0.1, top=0.9)\n"
        txt += "cbar_ax = fig.add_axes([0.86, 0.1, 0.03, 0.8])\n"
        txt += "\n"
        txt += "# Get x, y, z data\n"
        txt += f"x, y, z = {x}, {y}, {z}\n"
        res = self.settings['res']
        cmap = self.settings['colormap']
        if index_data == 1:
            txt += f"xi = np.linspace(min(x), max(x), {res})\n"
            txt += f"yi = np.linspace(min(y), max(y), {res})\n"
            txt += "X, Y = np.meshgrid(xi, yi)\n"
            txt += "Z = scgriddata((x, y), np.array(z), (X, Y), "
            txt += "method='linear', fill_value=max(z))\n"
            txt += "\n"
            txt += f"surf = ax.contourf(X,Y,Z, cmap='{cmap}')\n"
            txt += "cbar = fig.colorbar(surf, cax=cbar_ax)\n"
        else:
            txt += f"surf = ax.hist2d(x, y, bins={res},"
            txt += f" cmap='{cmap}', density=True)\n"
            txt += "cbar = fig.colorbar(surf[3], cax=cbar_ax)\n"

        txt += "plt.show()"

        if index_data == 0:
            zlabel = 'Density'
        else:
            zlabel = energy
        scriptfile = scriptfileformat.format(
            op1, op2, zlabel, self.settings['fol'])
        with open(os.path.join(self.folder, scriptfile), 'w',
                  encoding="utf8") as output:
            output.write(txt)
        self.statusbar.showMessage(f'Script file saved: {scriptfile}')

    def _load_file(self):  # pragma: no cover
        """Load data file.

        Function that sets up QObject, moves to QThread, and begins the
        data extract. Or, loads pre-processed data from file and
        moves to a QThread.
        """
        if self.folder is None and self.iofile is None:
            # When VisualApp is called without a directory, opens a filedialog
            iofile = QtWidgets.QFileDialog.getOpenFileName(
                parent=self,
                caption="Select input/output file in simulation directory")
            self.iofile = iofile[0]
            folder = os.path.dirname(os.path.realpath(self.iofile))
            self.folder = folder
            os.chdir(self.folder)
        # Setting name of file to QLineEdit of VisualApp
        if self.rst_file:
            self.filenameLine.setText(self.rst_file)
            self._load_data_output()
            return

        self.filenameLine.setText(self.iofile)
        # If PyVisA-data has been loaded
        if self.iofile.endswith(tuple(TRJ_FORMATS)) or \
                self.trajfile and \
                self.trajfile.endswith(tuple(TRJ_FORMATS)):
            self._load_data(self.iofile)
        elif self.iofile.endswith('.rst'):
            self._load_data_output()
        elif self.iofile.endswith(('.hdf5', '.zip')):
            self._load_data(self.iofile)
        else:
            msg = f'Format Error, file {self.iofile} not recognized.'
            self.statusbar.showMessage(msg)

    def _save_sim_data_hdf5_z(self):  # pragma: no cover
        """Save the data to hdf5 zipped file."""
        self._save_sim_data_hdf5()
        file_o = 'pyvisa_compressed_data.hdf5'
        pyvisa_zip(file_o)

    def _save_sim_data_hdf5(self):  # pragma: no cover
        """Save the data to hdf5 file."""
        file_o = 'pyvisa_compressed_data.hdf5'
        self.statusbar.showMessage(f'Saving data to {file_o}')
        data = pd.Series(self.dataobject.traj_data)
        infos = pd.Series(self.dataobject.infos)
        data.to_hdf(file_o, key='data', mode='w')
        infos.to_hdf(file_o, key='infos')
        self.statusbar.showMessage(f'Figure saved as {file_o}')

    def _load_data_output(self):  # pragma: no cover
        """Load simulation data.

        Function that loads simulation data using in/out file of simulation
        and initializes PathDensity on the directory, before showing the
        results in the VisualApp window.
        """
        self.statusbar.showMessage('Loading data from output files')
        # Set-up data object using PathDensity class and moving to thread
        self.thread = QtCore.QThread()
        self.dataobject = DataObject(iofile=self.iofile)
        # Connecting cycle print to update_cycle func
        self.dataobject.cycle_printed.connect(self.update_cycle)
        self.dataobject.return_coords.connect(self.update_fig)
        self.send_settings.connect(self.dataobject.get_xyz_data)
        self.start_cmd.connect(self.dataobject.walk)
        self.dataobject.moveToThread(self.thread)
        # Starting thread
        self.thread.start()
        # Starting QObject's walk_Dirs in thread
        self.start_cmd.emit('')

    def _load_data(self, pfile):  # pragma: no cover
        """Load saved data.

        Function that loads simulation data from a pre-compiled
        file and initializes the VisualApp window for plotting said data.
        """
        if pfile.endswith(('.hdf5', '.zip')):
            self.statusbar.showMessage('Loading data from compressed file.')
        else:
            self.statusbar.showMessage('Loading data from trajectory.')
        # Set-up data object using PathDensity class and moving to thread
        self.thread = QtCore.QThread()
        self.dataobject = VisualObject(pfile=pfile,
                                       trajfile=self.trajfile)
        # Connecting cycle print to update_cycle func
        self.dataobject.cycle_printed.connect(self.update_cycle)
        self.dataobject.return_coords.connect(self.update_fig)
        # Function takes hdf5 formats
        self.start_cmd.connect(self.dataobject.get_data)
        self.send_settings.connect(self.dataobject.get_xyz_data)
        self.dataobject.moveToThread(self.thread)
        # Starting thread
        self.thread.start()
        self.start_cmd.emit(pfile)

    def _reload(self):  # pragma: no cover
        """Reload the data of VisualApp.

        Function that clears the old data of VisualApp and initializes
        the load of new.
        """
        self.folder, self.iofile = None, None
        self.trajfile, self.rst_file = None, None
        self.dataobject = None
        self.toggle_buttons(False)
        self._load_file()

    def _init_widget(self):
        """Initialize QDropWidget.

        Function initializing the settings of the QDropWidget of VisualApp
        window. Adds correct items to the ComboBoxes of op1,op2 and E, and
        add the CustomFigCanvas object to the QFrame.
        If no folder or file is given when initializing VisualApp, a
        QFileDialog is opened to get file/folder.
        """
        self.statusbar.showMessage('No data loaded')
        # Actions and menubar
        self.actionExit.triggered.connect(self.close_event)
        # Save figure as ...
        self.action_png.triggered.connect(self.save_png)
        self.action_makefig_script.triggered.connect(self.save_script)
        # Save figure data as ...
        self.action_hdf5.triggered.connect(self.save_hdf5)
        self.action_json.triggered.connect(self.save_json)
        self.action_datafile.triggered.connect(self.save_textdata)
        self.action_Load_data.triggered.connect(self.action_reload)
        # Save data as ...
        self.action_sim_hdf5.triggered.connect(self._save_sim_data_hdf5)
        self.action_sim_hdf5_zip.triggered.connect(self._save_sim_data_hdf5_z)
        # Connect show reg.line and interfaces, and method check
        self.intShowChkBtn.stateChanged.connect(self.toggle_intf)
        self.regLineChkBtn.stateChanged.connect(self.toggle_regl)
        self.cbarChkBtn.stateChanged.connect(self.toggle_cbar)
        self.showTitleChkBtn.stateChanged.connect(self.toggle_titles)
        self.plotTypeComBox.activated.connect(self.disable_zcombox)
        self.plotTypeComBox.activated.connect(self.disable_reglinecheck)
        # The Frame
        self.mainFrame.setStyleSheet(
            "QWidget { background-color: %s }" % QtGui.QColor(
                210, 210, 235, 255).name()
        )
        self.layout_a = QtWidgets.QGridLayout()
        self.mainFrame.setLayout(self.layout_a)
        # The Figure
        self.myfig = CustomFigCanvas()
        self.layout_a.addWidget(self.myfig)
        # Setting up the dropwidget of VisualApp
        self.show()
        # Refresh button
        self.refreshBtn.clicked.connect(self.refresh_data)
        # Bind update button to function
        self.updateBtn.clicked.connect(self.emit_settings)
        self.saveFigBtn.clicked.connect(self.save_png)
        self.previewBtn.clicked.connect(self.update_preview)
        # Bind <ENTER>-press in colormap line to function
        self.cmapComBox.lineEdit().returnPressed.connect(self._change_cmap)
        # Bind changes to title and labels to functions
        self.titleSizeSpin.valueChanged.connect(self._update_canvas_font)
        self.axesSizeSpin.valueChanged.connect(self._update_canvas_font)
        self.titleLine.returnPressed.connect(self._update_canvas_text)
        self.xaxisLabel.returnPressed.connect(self._update_canvas_text)
        self.yaxisLabel.returnPressed.connect(self._update_canvas_text)
        self.zaxisLabel.returnPressed.connect(self._update_canvas_text)
        # Show only stored trjs
        self.storedTrjsChkBtn.stateChanged.connect(
            self.toggle_only_stored_trjs)
        # Interactiv tab
        self.previewBtn_2.clicked.connect(self.emit_settings)
        self.drawButton.clicked.connect(self.draw_trajectory)
        self.playBtn.clicked.connect(self.play_cycle)
        self.store_trj_btn.clicked.connect(self.create_full_traj)
        # Recalculation
        self.actionRecalculate_data.triggered.connect(self.recalculate)
        # Analysis tab
        self.clusterButton.clicked.connect(self.statistic_analysis)
        self.correlationBtn.clicked.connect(self.show_correlation)
        self.RandForestBtn.clicked.connect(self.random_forest)
        self.DecisionTreeBtn.clicked.connect(self.decision_tree)
        self.PCABtn.clicked.connect(self.initiate_pca)

    def _change_cmap(self):
        """Set colormap of surface.

        Function that tries to set the colormap of object myfig.surf (plot)
        Bound to detection of returnPress in colormap combo box.
        """
        if self.settings is None or self.myfig.surf is None:
            self.statusbar.showMessage('Figure not ready!')
            return
        meth = self.settings['method']
        col = self.cmapComBox.currentText()

        try:
            # Test if entered colormap is valid
            _ = plt.get_cmap(col)
            if meth[0] == 'Density':
                # Method is 2dhist-plot/list where 4th item is the plot/surf.
                self.myfig.surf[3].set_cmap(col)
            elif meth[0] == 'Scatter':
                # If scatter, replot.
                self.statusbar.showMessage(
                    'Scatter plots have to be re-drawn to update colors!')
            else:
                # Should work for most surf-objects.
                self.myfig.surf.set_cmap(col)
        except AttributeError:
            # If surf has no .set_cmap()
            self.statusbar.showMessage(
                'Could not recognize colormap, try to "Update" figure...')
        except TypeError:
            # If no surf-object, or NoneType in general
            self.statusbar.showMessage(
                'Figure not recognized, try to "Update" figure...')
        except ValueError:
            # Did not recognize 'col' in cm
            self.statusbar.showMessage(
                'Chosen colormap not recognized!')
        self.myfig.fig.canvas.draw()

    def _get_op_range(self):
        """Display the range of the order parameter to the GUI."""
        interfaces = self.dataobject.infos['interfaces']
        self.minDataRange.setText(str(interfaces[0]))
        self.maxDataRange.setText(str(interfaces[-1]))

    def _change_zoom(self):
        """Set the zoom/x-&y-limits of the plot."""
        if self.settings is None or self.myfig.cbar is None:
            self.statusbar.showMessage('Figure not ready!')
            return
        self.statusbar.showMessage('Drawing figure...')
        if self.settings['dim'] == 3:
            self.myfig.ax.set_xlim3d(*self.settings['x-limits'], )
            self.myfig.ax.set_ylim3d(*self.settings['y-limits'], )
            self.myfig.ax.set_zlim3d(*self.settings['z-limits'], )
        elif self.settings['dim'] == 2:
            self.myfig.ax.set_xlim(*self.settings['x-limits'], )
            self.myfig.ax.set_ylim(*self.settings['y-limits'], )
        # Change range on colorbar
        # Colorbar of contour plots does not scale, but color is updated
        norm = mpl.colors.Normalize(vmin=self.settings['z-limits'][0],
                                    vmax=self.settings['z-limits'][1])
        self.myfig.cbar.mappable.set_norm(norm=norm)
        _redraw_colorbar(self.myfig.cbar)
        self.myfig.fig.canvas.draw()
        self.statusbar.showMessage('Plot ready!')

    def update_preview(self):
        """Update the visual of the plot.

        Button function, updates visual of plot based on all settings of
        'Plot'-tab in VisualApp.
        """
        if self.settings is None or self.myfig.surf is None:
            self.statusbar.showMessage('Figure not ready!')
            return
        self.statusbar.showMessage('Drawing figure...')
        self._get_settings()
        self._update_canvas_text()
        self._change_zoom()
        self._change_cmap()
        self.myfig.fig.canvas.draw()
        self.statusbar.showMessage('Plot ready!')

    # Defining button functions that toggle the visibility of
    # interfaces, regression line, colorbar and title/axis labels.
    def toggle_intf(self):
        """Toggle interface visibility."""
        if not self.myfig.intf:
            return
        select = self.intShowChkBtn.isChecked()
        for line in self.myfig.intf:
            line.set_visible(select)
        self.myfig.fig.canvas.draw()

    def toggle_regl(self):
        """Toggle regression line and legend visibility."""
        if self.settings is None or not self.myfig.regl:
            return
        if self.settings['dim'] == 2:
            select = self.regLineChkBtn.isChecked()
            self.myfig.regl[0].set_visible(select)
        self.myfig.fig.canvas.draw()

    def toggle_cbar(self):  # pragma: no cover
        """Toggle colorbar visibility."""
        select = self.cbarChkBtn.isChecked()
        self.myfig.cbar_ax.set_visible(select)
        self.myfig.fig.subplots_adjust(left=0.1, right=0.85,
                                       bottom=0.1, top=0.9)
        self.myfig.fig.canvas.draw()

    def toggle_titles(self):
        """Toggle labels and titles visibility."""
        if any(item is None for item in (
            self.myfig.title,
            self.myfig.xaxis,
            self.myfig.yaxis,
            self.myfig.zaxis,
        )):
            return
        select = self.showTitleChkBtn.isChecked()
        self.myfig.title.set_visible(select)
        self.myfig.xaxis.set_visible(select)
        self.myfig.yaxis.set_visible(select)
        self.myfig.zaxis.set_visible(select)
        self.myfig.fig.canvas.draw()

    def toggle_only_stored_trjs(self):
        """Toggle to show only data with trajectory files."""
        if self.dataobject is None:
            self.statusbar.showMessage('No data loaded')
            return
        self.emit_settings()

    def emit_settings(self):  # pragma: no cover
        """Update settings before sending to the data object.

        Function calls for an update of data/plot settings before sending
        dictionary to data object in main worker to update the data lists
        x, y, (z).

        """
        self.toggle_buttons(False)
        # Updating status bar of VisualApp window
        self.statusbar.showMessage('Loading data...')
        self._get_settings()
        if 'Contour' in self.settings['method'] and \
                self.settings['E'] == 'None':
            self.statusbar.showMessage('Invalid combination, '
                                       'contours need z-values')
            self.toggle_buttons(True)
        elif self.settings['method'][0] in ['Surface', 'Contour'] and \
                self.settings['op1'] == self.settings['op2']:
            self.statusbar.showMessage('Invalid combination, '
                                       'two axes are identical')
            self.toggle_buttons(True)
        else:
            self.send_settings.emit(self.settings)

    @QtCore.pyqtSlot(list)
    def update_cycle(self, cycles):  # pragma: no cover
        """Set upper and lower bound of cycles for data loading.

        Function bound to pyqtSignal from DataObject in QThread, sets upper
        and lower bound of cycles for data loading, and enables 'Update'
        QPushButton when finished.

        Parameters
        ----------
        cycles : list
            Integers defining the initial en final cycle number to plot.

        """
        # Adding items to ComboBoxes in VisualApp.Window
        self.firstOpComBox.clear()
        self.firstOpComBox.addItems(self.dataobject.infos['op_labels'])
        self.firstOpComBox.insertSeparator(
            len(self.dataobject.infos['op_labels']) + 1)
        self.firstOpComBox.addItems(ENERGYLABELS)
        self.secondOpComBox.clear()
        self.secondOpComBox.addItems(self.dataobject.infos['op_labels'])
        self.secondOpComBox.insertSeparator(
            len(self.dataobject.infos['op_labels']) + 1)
        labels = []
        for label in ENERGYLABELS:
            if label not in self.dataobject.infos['op_labels']:
                labels.append(label)
        self.secondOpComBox.addItems(labels)
        self.energyComBox.clear()
        self.energyComBox.addItems(self.dataobject.infos['op_labels'])
        self.energyComBox.insertSeparator(
            len(self.dataobject.infos['op_labels']) + 1)
        self.energyComBox.addItems(ENERGYLABELS)
        self.energyComBox.addItems(['None'])
        self.folderComBox.clear()
        self.folderComBox.addItems(self.dataobject.infos['ensemble_names'])
        self.folderComBox.insertSeparator(
            len(self.dataobject.infos['ensemble_names']) + 1)

        self.folderComBox.addItem('All')
        # Setting cycles in dropwidget
        minc, maxc = cycles[0], cycles[1]
        self.minCyclSpin.setRange(minc, maxc)
        self.maxCyclSpin.setRange(minc, maxc)
        self.minCyclSpin.setValue(minc)
        self.maxCyclSpin.setValue(maxc)
        # Enables the buttons and updates status
        self.statusbar.showMessage(
            'All %s cycles loaded from files, ready to plot!', maxc - minc)
        self.updateBtn.setEnabled(True)
        self.refreshBtn.setEnabled(True)
        self._get_op_range()

    @QtCore.pyqtSlot(list, list, list)  # noqa: C901
    def update_fig(self, x, y, z):  # pragma: no cover
        """Update figure by pyqtSignal.

        Function that updates figure canvas of CustomFigCanvas and the
        colorbar. Initiated by pyqtSignal from running QThread.
        Function checks method chosen in GUI and calls on built-in functions
        of orderparam_density module to generate a surface/contour/scatter/etc.
        If checked, function also plots interfaces and regression line.

        Parameters
        ----------
        x, y, z : list
            Floats, the coordinates used in plotting.

        Updates/Draws
        -------------
        myfig.fig : Re-draws canvas with method, data, and possible
            interfaces, reg.line and legend.
        myfig.ax : As above. (Axes containing plot, legend, lines, etc)
        myfig.cbar_ax : Axes with colorbar of plot.

        """
        self.statusbar.showMessage('re-Drawing figure...')

        remove_nan([x, y, z])
        order_parameter = self.dataobject.infos['op_labels'][0]
        l_x, l_y, l_z = len(x), len(y), len(z)
        fail = False
        if l_x == 0 or l_y == 0:
            if self.storedTrjsChkBtn.isChecked():
                self.display_message_box('Cannot plot...',
                                         'No stored trajectory files.')
                self.storedTrjsChkBtn.setChecked(False)
            else:
                self.statusbar.showMessage(
                    'One (or more) lists are empty, plotting halted!' +
                    f'x: {l_x}, y: {l_y}')
            fail = True
        if 'Density' not in self.settings['method'][0]:
            if l_x != l_z or l_y != l_z:
                self.statusbar.showMessage(
                    'Lists have different lengths! ' +
                    f'x: {l_x}, y: {l_y}, z: {l_z}')
                fail = True
        if self.settings['dim'] == 3 and l_z == 0:
            self.statusbar.showMessage('Z list empty, nothing to plot!')
            fail = True

        if fail:
            self.updateBtn.setEnabled(True)
            self.refreshBtn.setEnabled(True)
            return

        # Do a myfig.set_up(), in case dim/projection has changed
        self.myfig.set_up(dim=self.settings['dim'],
                          cbar=self.settings['colorbar'])
        # Setting default xlim and ylim
        self._set_def_limits(x, y, z)

        plot = self.settings['method'][0].lower()
        if self.settings['method'][-1] == '(filled)':
            plot = plot + 'f'

        if self.settings['dim'] == 3:
            self.myfig.surf, self.myfig.cbar = gen_surface(
                x, y, z,
                self.myfig.fig,
                self.myfig.ax,
                cbar_ax=self.myfig.cbar_ax,
                dim=self.settings['dim'],
                method=plot,
                res_x=self.settings['res'],
                res_y=self.settings['res'],
                colormap=self.settings['colormap'])
        elif self.settings['dim'] == 2:
            if self.settings['method'][0] == 'Density':
                if l_x != l_y:
                    min_length = min(l_x, l_y)
                    x, y = x[:min_length], y[:min_length]
                    self.statusbar.showMessage(
                        'Lists have different lengths! Plotting smallest.. ' +
                        f'x: {l_x}, y: {l_y}')

                self.myfig.surf = self.myfig.ax.hist2d(
                    x, y, (self.settings['res'], self.settings['res']),
                    cmap=self.settings['colormap'], density=True)
                self.myfig.cbar = self.myfig.fig.colorbar(
                    self.myfig.surf[3], cax=self.myfig.cbar_ax)
                # Set zmin and zmax equal to min and max value in histogram
                self.zminLine.setText(str(self.myfig.surf[0].min()))
                self.zmaxLine.setText(str(self.myfig.surf[0].max()))
            else:
                self.myfig.surf, self.myfig.cbar = gen_surface(
                    x, y, z,
                    self.myfig.fig,
                    self.myfig.ax,
                    cbar_ax=self.myfig.cbar_ax,
                    dim=self.settings['dim'],
                    method=plot,
                    res_x=self.settings['res'],
                    res_y=self.settings['res'],
                    colormap=self.settings['colormap'])
        # Showing interface planes(3D)/lines(2D)
        self.myfig.intf = []
        if self.settings['op1'] == order_parameter:
            for i in self.dataobject.infos['interfaces']:
                if i > self.xmax:
                    break
                if self.settings['dim'] == 3:
                    interf = plot_int_plane(
                        self.myfig.ax,
                        i, (self.ymin, self.ymax),
                        (self.zmin, self.zmax),
                        visible=self.settings['show_int'])
                else:
                    interf = self.myfig.ax.axvline(
                        i, linewidth=1,
                        color='grey', alpha=0.5,
                        visible=self.settings['show_int'])

                self.myfig.intf.append(interf)

        # Showing regression line(2D only)
        self.myfig.regl = []
        if self.settings['dim'] == 2:
            self.myfig.regl = plot_regline(self.myfig.ax, x, y)
            self.myfig.regl[0].set_visible(self.settings['reg_line'])
            self.myfig.legend = self.myfig.ax.legend()
            self.myfig.legend.set_visible(self.settings['reg_line'])
            self.myfig.cbar_ax.set_visible(self.settings['colorbar'])
        self._get_titles()
        # Update CustomFigure canvas and status
        self.myfig.fig.canvas.draw()
        self.statusbar.showMessage('Plot ready!')
        self.toggle_buttons(True)

    def on_pick(self, event):  # pragma: no cover
        """Write all information about the closest point to slots in PyVisA.

        Parameters
        ----------
        event : Matplotlib mouse event
            It activates the interactive function on the plot.

        """
        # If no data has been loaded
        if not self.dataobject:
            self.display_message_box('Cannot select!',
                                     'No data has been loaded.')
            return

        # Avoid picking with 3D plots
        self._get_settings()
        if self.settings['dim'] == 3:
            return

        # Remove old points
        _remove_last_line(self.myfig.ax)

        # Avoid MemoryError from scipy cdist
        try:
            closest, cycle, ensemble = self.dataobject. \
                find_cycle_ensemble(event)
        except MemoryError:
            self.display_message_box("Error...",
                                     "Dataset is too large")
            return

        # If no data is plotted
        if len(closest) == 0:
            self.display_message_box("Error...",
                                     "Please plot data before picking")
            return

        self.xLine.setText(str(round(closest[0], 2)))
        self.yLine.setText(str(round(closest[1], 2)))
        self.cycLine.setText(str(cycle))
        self.ensLine.setText(ensemble)
        traj = self.dataobject.traj_data[ensemble][cycle]
        if 'status' in traj.info:
            self.statusLine.setText(traj.info['status'])
        if 'MC-move' in traj.info:
            self.mcLine.setText(traj.info['MC-move'])

        self.myfig.ax.plot(closest[0], closest[1],
                           color=self.colorBox.currentText(),
                           marker='o',
                           markersize=self.sizeSpinBox.value(),
                           label='point')
        self.myfig.fig.canvas.draw()

    def draw_trajectory(self):  # pragma: no cover
        """Draw chosen trajectory on plot.

        Show a full chosen trajectory on the plot in a selected
        color with choices of marker size and line thickness.

        """
        colors = {'Green': 'g', 'Blue': 'b', 'Red': 'r',
                  'Cyan': 'c', 'Magenta': 'm', 'Yellow': 'y',
                  'Black': 'k', 'White': 'w'}

        # Remove old figures in the plot
        _remove_last_line(self.myfig.ax)

        if not self.xLine.text():
            self.display_message_box('Cannot draw...', 'Select a trajectory')
            return

        self.myfig.fig.canvas.draw()
        criteria = {'cycles': [int(self.cycLine.text()),
                               int(self.cycLine.text())],
                    'ensemble_name': self.ensLine.text(),
                    'x': self.settings['op1'],
                    'y': self.settings['op2'],
                    'status': self.settings['ACC']}
        x, y, _, _ = self.dataobject.load_traj(criteria)
        if x and y:
            self.myfig.ax.plot(x, y,
                               color=colors[self.colorBox.currentText()],
                               marker='o',
                               linestyle='-',
                               linewidth=self.widthSpinBox.value(),
                               markersize=self.sizeSpinBox.value())
            self.myfig.fig.canvas.draw()
            self.statusbar.showMessage('Plot ready!')
        else:
            self.statusbar.showMessage(
                'Different length of lists, drawing halted!' +
                ' x: %s, y: %s', self.settings['op1'], self.settings['op2'])

    def create_full_traj(self, cycle=None, ensemble_name=None,
                         rst_file=None):  # pragma: no cover
        """Create trajectory file.

        Parameters
        ----------
        cycle: string, optional
            Cycle number of trajectory.
        ensemble_name: string, optional
            Ensemble name of trajectory.
        rst_file: string, optional
            Path to the input/output file of the simulation.

        """
        # If the function has been called through the save_trj button
        if not (cycle or ensemble_name or rst_file):
            cycle = self.cycLine.text()
            ensemble_name = self.ensLine.text()
            if not cycle:
                self.display_message_box('Cannot store trajectory!',
                                         'Please select a trajectory first.')
                return None
            if not isinstance(self.iofile, PathDensity) and \
                    self.iofile.endswith('.rst'):
                rst_file = self.iofile
            elif self.rst_file:
                rst_file = self.rst_file
            else:
                self.display_message_box('Cannot store trajectory!',
                                         'Please load PyVisA with a .rst '
                                         'file in the simulation directory.')
                return None

        # Find trajectory in simulation directory
        traj = self.dataobject.traj_data[ensemble_name][int(cycle)]
        if 'status' not in traj.info.keys() or not traj.info['stored']:
            self.display_message_box('Trajectory files not present.',
                                     'Pick another trajectory')
            return None

        # Find path to folders
        trr_status = 'traj-acc' if traj.info['status'] == 'ACC'\
            else 'traj-rej'

        path_to_cycle = os.getcwd() +\
            f'/{ensemble_name}/traj/{trr_status}/{cycle}'
        path_to_traj = path_to_cycle + '/traj'
        # Find the sequence and velocity of files in the trajectory
        files = read_traj_txt_file(path_to_cycle + '/traj.txt')

        # Find the topology file from rst_file
        trj_list = []
        input_settings = \
            settings.parse_settings_file(rst_file)['engine']['input_files']
        if input_settings['conf'].endswith('.g96'):
            self.display_message_box('Cannot load file..',
                                     'Configuration format not supported.')
            return input_settings['conf']

        # Create the correct order of traj-files
        for _, value in files.items():
            if not value[0].endswith(tuple(TRJ_FORMATS)):
                continue
            if int(value[1]) == -1:
                trj_list.append(self.rev_trj(md.load(
                    path_to_traj + f'/{value[0]}', top=input_settings['conf']
                )))
            else:
                trj_list.append(md.load(
                    path_to_traj + f'/{value[0]}', top=input_settings['conf']
                ))

        trj = trj_list[0]
        for _ in trj_list[1:]:
            trj += _

        trj_name = f'trajectory_{ensemble_name}_{cycle}.pdb'
        trj.save_pdb(trj_name)
        return trj_name

    def play_cycle(self):  # pragma: no cover
        """Collect traj-files in order to animate the trajectory."""
        # Find cycle and ensemble name of chosen trj
        cycle = self.cycLine.text()
        ensemble_name = self.ensLine.text()

        if not cycle or not ensemble_name:
            self.display_message_box('Cannot play...',
                                     'No trajectory chosen.')
            return

        if not isinstance(self.iofile, PathDensity) and \
                self.iofile.endswith('.rst'):
            trj_name = self.create_full_traj(cycle, ensemble_name,
                                             self.iofile)
        elif self.rst_file:
            trj_name = self.create_full_traj(cycle, ensemble_name,
                                             self.rst_file)
        else:
            self.display_message_box('Invalid choice..',
                                     'Please load an rst file')
            return

        if trj_name:
            try:
                subprocess.run(  # nosec B603 B607
                    ["avogadro2", trj_name], check=True)
                os.remove(os.getcwd() + '/' + trj_name)
            except FileNotFoundError:
                self.display_message_box('Cannot play cycle',
                                         'Could not start avogadro2.' + '\n' +
                                         ' Make sure avogadro is installed' +
                                         ' and try again.')
                return

    @staticmethod
    def rev_trj(trj):
        """Reverse a trj via MD traj.

        Parameters
        ----------
        trj : mdtraj object
            Mdtraj MD trajectory.

        Returns
        -------
        r_trj : mdtraj object
            The reversed MD trajectory.

        """
        f_n = np.arange(len(trj) - 1, 0, -1)
        r_trj = trj.slice(f_n)
        return r_trj

    def display_message_box(self, line1, line2):
        """Display a message box with chosen text.

        Parameters
        ----------
        line1: string
            Title of message box.
        line2: string
            Main message for the box.

        """
        msg = QtWidgets.QMessageBox.question(self, line1,
                                             line2,
                                             QtWidgets.QMessageBox.Ok)
        if msg == QtWidgets.QMessageBox.Ok:
            return

    def recalculate(self):  # pragma: no cover
        """Recalculate OP and new CV's.

        Using the functions from recalculate_order.py this function
        will calculate new collective variables to the simulation while
        keeping the original OP in order to explore new descriptors of the
        system. The new files will replace the order files so PyVisA can
        load and visualize the new descriptors. The old files will
        be renamed as order_old.

        """
        msg = QtWidgets.QMessageBox.question(
            self, "Warning: Confirm descriptor recalculation...",
            "This operation might result in data-loss, " +
            "\n" + "are you sure to continue?",
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)

        if msg == QtWidgets.QMessageBox.Yes:

            if not self.rst_file or 'rst' not in self.rst_file:
                self.display_message_box('Cannot perform re-calculations...',
                                         'Please load PyVisA with a rst file.')
                return
            if not recalculate_all(self.folder, self.rst_file,
                                   self.dataobject.infos['ensemble_names']):
                self.display_message_box('Cannot perform re-calculations...',
                                         'Invalid input.')
                return
            self.statusbar.showMessage('Recalculating...')
            self.statusbar.showMessage('Plot ready!')
            self.refresh_data()
            self.display_message_box('Recalculation complete!',
                                     'The new data can now be visualized.')

    def statistic_analysis(self):  # pragma: no cover
        """Perform statistical analysis on plotted data.

        Function which performs a range of different clustering
        techniques like K-mean, Hierarchical, Gaussian Mixture,
        Spectral and PCA by choice of the user. The user specifies
        the type of analysis as well as the amount of cluster/components.

        """
        self.statusbar.showMessage('Performing analysis...')
        plt.ion()
        if not self.dataobject:
            self.display_message_box('Cannot analyse!',
                                     'No data has been loaded.')
            return

        self._get_settings()
        n_clusters = self.nrOfCluster.value()
        method = self.clusterTypeComBox.currentText()
        cmap = self.cmapComBox.currentText()

        if self.settings['dim'] == 3:
            self.display_message_box('Cannot perform statistics!',
                                     'Please choose a 2D plot.')
            return
        # If no data has been plotted:
        if not self.dataobject.settings:
            self.display_message_box('Cannot perform analysis...',
                                     'No data has been selected.')
            return
        data = np.column_stack([self.dataobject.x,
                                self.dataobject.y])

        methods = {'Kmeans': k_means,
                   'Hierarchical': hierarchical,
                   'Gaussian Mixture': gaussian_mixture,
                   'Spectral': spectral}
        try:
            methods[method](n_clusters, data, self.settings, cmap)
        except MemoryError:
            self.display_message_box('Dataset too large',
                                     'Cannot perform fitting.')
            return
        plt.ioff()
        self.statusbar.showMessage('Analysis completed!')

    def initiate_pca(self):  # pragma: no cover
        """Perform PCA."""
        self._get_settings()
        n_pca = self.nrOfpca.value()
        cmap = self.cmapComBox.currentText()
        if n_pca > len(self.dataobject.infos['op_labels']) + 2:
            self.display_message_box('Cannot perform PCA, ',
                                     'The amount of components exceeds' +
                                     '\n' + 'the amount of features in'
                                            'the dataset')
        else:
            # Collect data for analysis
            data = self.dataobject.load_to_df(self.settings)
            try:
                pyvisa_pca(n_pca, self.settings, data, cmap)
            except ValueError:
                self.display_message_box('Cannot perform analysis!',
                                         'No data could be used with the '
                                         'chosen criteria. Try to change '
                                         'the status type of MC-move or '
                                         'type of path.')

    def random_forest(self):  # pragma: no cover
        """Create a random forest from simulation results."""
        self.statusbar.showMessage('Performing analysis...')
        self._get_settings()
        min_op = float(self.minDataRange.text())
        max_op = float(self.maxDataRange.text())
        try:
            data, reactive = \
                self.dataobject.remove_values(self.settings, min_op, max_op)
        except AssertionError:
            self.display_message_box('Cannot perform analysis!',
                                     'No data could be used with the chosen'
                                     ' criteria. Try to change the status, '
                                     'type of MC-move or type of path.')
            return
        if len(data) == 0:
            self.display_message_box('Cannot perform analysis!',
                                     'No data was valid within the '
                                     'chosen data range.')
            return

        depth = self.LayerSpinBox.value()
        random_forest(data, reactive, depth=depth)
        self.statusbar.showMessage('Analysis completed!')

    def decision_tree(self):  # pragma: no cover
        """Create a decision tree from simulation results."""
        self.statusbar.showMessage('Performing analysis...')
        self._get_settings()
        min_op = float(self.minDataRange.text())
        max_op = float(self.maxDataRange.text())
        try:
            data, reactive = \
                self.dataobject.remove_values(self.settings, min_op, max_op)
        except AssertionError:
            self.display_message_box('Cannot perform analysis!',
                                     'No data could be used with the chosen'
                                     ' criteria. Try to change the status, '
                                     'type of MC-move or type of path.')
            return

        if len(data) == 0:
            self.display_message_box('Cannot perform analysis!',
                                     'No data was valid within the '
                                     'chosen data range.')
            return
        depth = self.LayerSpinBox.value()
        decision_tree(data, reactive, depth=depth)
        self.statusbar.showMessage('Analysis completed!')

    def show_correlation(self):  # pragma: no cover
        """Show the correlation between simulation data.

        Function which plots a correlation matrix from the simulation
        data.

        """
        self.statusbar.showMessage('Performing analysis...')
        plt.ion()
        self._get_settings()
        if not self.dataobject:
            self.display_message_box('Cannot analyse!',
                                     'No data has been loaded.')
            return
        data = self.dataobject.load_to_df(self.settings)
        correlation_matrix(data)
        plt.ioff()
        self.statusbar.showMessage('Analysis completed!')

    def refresh_data(self):  # pragma: no cover
        """Signal a walk of the data."""
        if not self.rst_file and \
                self.iofile.endswith(('.hdf5', '.zip', '.rst')):
            self.display_message_box('Cannot refresh.',
                                     'Please load PyVisA with a .rst file \n')
            return
        if self.trajfile and self.iofile.endswith('.rst'):
            self.dataobject = DataObject(iofile=self.iofile)
        elif self.rst_file:
            self.dataobject = DataObject(iofile=self.rst_file)
        else:
            self.dataobject = DataObject(iofile=self.iofile)
        self.dataobject.cycle_printed.connect(self.update_cycle)
        self.dataobject.return_coords.connect(self.update_fig)
        self.send_settings.connect(self.dataobject.get_xyz_data)
        self.dataobject.walk()


class CustomFigCanvas(FigureCanvas):  # pylint: disable=too-few-public-methods
    """
    Class definition of the custom figure canvas used in VisualApp.

    Attributes
    ----------
    surf : Placeholder for object shown in main subplot of figure.
    ax : Main subplot (Axes object) of CustomFigCanvas.
    cbar_ax : Subplot (Axes object) of colorbar in CustomFigCanvas.
    intf : Placeholder for a list of interface lines/planes in the figure.
    regl : Placeholder for 2D regression line in figure.
    legend : Placeholder for legend in figure.
    title : Placeholder for the figure title.
    {x/y/z} axis : Placeholder for the x/y/z-axis labels.

    """

    def __init__(self, ):
        """Initialize the FigureCanvas class."""
        self.fig = Figure()
        FigureCanvas.__init__(self, self.fig)
        self.surf = None
        self.cbar = None
        self.cbar_ax = None
        self.intf = None
        self.regl = None
        self.legend = None
        self.title = None
        self.ax = None
        self.xaxis, self.yaxis, self.zaxis = None, None, None
        self.set_up()

    def set_up(self, dim=2, cbar=True):
        """Define subplot projection (2D/3D).

        Function that defines the main subplot's projection (2D/3D)
        depending on chosen method from VisualApp GUI.

        Parameters
        ----------
        dim : integer, optional
              Dimension/projection of main subplot.
        cbar : boolean, optional
              If True, the space for the colorbar is considered
              in the setting.

        """
        # The Figure
        self.fig.clf()
        # Subplots, depending on dimension/projection
        if dim == 3:  # pragma: no cover
            self.ax = self.fig.add_subplot(111, projection='3d')
        elif dim == 2:
            self.ax = self.fig.add_subplot(111)
        # Axes for colorbar
        self.cbar_ax = self.fig.add_axes([0.86, 0.1, 0.03, 0.8])
        if cbar:
            self.fig.subplots_adjust(left=0.1, right=0.85, bottom=0.1, top=0.9)
        else:
            self.fig.subplots_adjust(left=0.1, right=0.90, bottom=0.1, top=0.9)
        self.cbar_ax.set_visible(False)


class DataSlave(QtCore.QObject, PathVisualize):
    """QObject class definition that holds the PathDensity data.

    By pyqtSignals this object can be called upon to either do a directory
    "walk" to collect data, or get specific data into lists by specific
    settings.

    Attributes
    ----------
    Contains the same attributes as the PathDensity() class of
    orderparam_density.py in source.

    Functions
    ---------
    load_all: function
        Load x,y and z data from all trajectories selected by the settings
        of PyVisA. Also stored the coordinates, ensemble name and cycle
        of each point for use in interactivity in the plot.
    find_cycle_ensemble : function
        Accept a MouseClick event on the
        plot and return the x, y coordinate, cycle and ensemble name of the
        closest point to wherever the user clicks on the plot.
    load_to_df : function
        Loads trajectory data into a pandas dataframe for analysis.

    Signals
    -------
    cycle_printed : PyQt Signal
        Send a list of min and max cycle number to
        the VisualApp to update the dropwidget with correct values.
    returns_coords : PyQt Signal
        Send three lists of coords (or empty)
        to the VisualApp to update the figure.

    Slots
    -----
    walk : PyQt Slot connected to VisualApp
        When activated with a string,
        initializes the walk_Dirs() function of DataObject/PathDensity
        class. At end of function, emits list of cycles to signal
        cycle_printed.
    get : PyQt Slot connected to VisualApp
        When activated with a dictionary
        with relevant plotting settings, preforms a data fetch of either
        get_Edata or get_Odata on DataObject/PathDensity class.
        At end of function, emits lists of coords to signal return_coords.

    """

    return_coords = QtCore.pyqtSignal(list, list, list)

    def __init__(self, iofile=None):
        """Initialize the QObject class and classes inherited.

        Parameters
        ----------
        iofile : string, optional
            The input file name.

        """
        self.iofile = iofile
        QtCore.QObject.__init__(self)
        PathVisualize.__init__(self)
        self.settings = {}
        self.x, self.y, self.z = [], [], []
        self.data_origin = []

    def load_all(self):
        """Load chosen data from all trajectories.

        Function that uses the load_traj function from orderparam_density()
        to load all the x, y, and z data corresponding to the chosen
        settings for plotting. The criteria for plotting is collected
        from the settings and sent to load_traj. It also stores the
        coordinates, and data_origin of each points which is used for
        the interactivity in other functions like on_pick.

        Returns
        -------
        x : list
            List of data from chosen parameter.
        y : list
            List of data from chosen parameter.
        z : list
            List of data from chosen parameter.

        """
        criteria = {'x': self.settings['op1'],
                    'y': self.settings['op2'],
                    'z': self.settings['E'],
                    'ensemble_name': self.settings['fol'],
                    'cycles': (self.settings['min_cycle'],
                               self.settings['max_cycle']),
                    'status': self.settings['ACC'],
                    'MC-move': self.settings['MC-move'],
                    'weight': self.settings['weight'],
                    'stored': self.settings['stored']}

        if self.settings['start_end'] == 'All':
            pass
        else:
            criteria['reactive'] = \
                self.settings['start_end'].lower() == 'reactive'

        x, y, z = [], [], []
        data_origin = []
        if criteria['ensemble_name'] == 'All':
            for fol in self.infos['ensemble_names']:
                criteria['ensemble_name'] = fol
                x_i, y_i, z_i, doi = self.load_traj(criteria)
                x.extend(x_i)
                y.extend(y_i)
                z.extend(z_i)
                data_origin.extend(doi)
        else:
            x, y, z, data_origin = self.load_traj(criteria)

        self.x, self.y, self.z = x, y, z
        self.data_origin = data_origin
        return x, y, z

    @QtCore.pyqtSlot(dict)
    def get_xyz_data(self, new_settings):  # pragma: no cover
        """Signal emission to get xyz data from settings in VisualApp.

        Parameters
        ----------
        new_settings : dict
            The settings provided by an user to select the data to extract.

        """
        self.settings = new_settings
        x, y, z = self.load_all()

        # If allowing data shift
        if self.settings['try_shift']:
            op1 = self.settings['op1']
            fixedx = op1 == self.infos['op_labels'][0]
            x, y = try_data_shift(x, y, fixedx)
        self.return_coords.emit(x, y, z)

    def find_cycle_ensemble(self, event):
        """Return x, y, cycle and ensemble name at the closest point.

        Parameters
        ----------
        event : Matplotlib mouse event
            The users mouse click on the plot.

        Returns
        -------
        closest : list
            List of the x and y coordinate of the closest point
            to the users click on the plot.
        cycle : int
            Corresponding cycle of the closest point.
        ensemble_name : string
            Corresponding ensemble name of the closest point.

        """
        # If no data is plotted
        if not (hasattr(self, 'x') and hasattr(self, 'y')):
            return [], '', ''
        if len(self.x) == 0 and len(self.y) == 0:
            return [], '', ''

        # Find the closest point
        data_list = list(zip(self.x, self.y))
        data_list = np.array(data_list)
        closest_index = distance.cdist([[event.xdata, event.ydata]],
                                       data_list,
                                       'seuclidean').argmin()

        closest = data_list[closest_index]
        # Find corresponding cycle and ensemble
        ensemble_name = self.data_origin[closest_index][0]
        cycle = self.data_origin[closest_index][1]
        return closest, cycle, ensemble_name

    def load_to_df(self, sim_settings):
        """Create a pandas dataframe for analysis.

        Parameters
        ----------
        sim_settings : dict
            Dict containing settings from GUI.

        Returns
        -------
        df : Pandas Dataframe
            Dataframe containing all OP, CV and energy data from
            selected ensembled.

        """
        frames = []
        number_variables = len(self.infos['op_labels']) + 3

        selection_crit = {'stored': sim_settings['stored'],
                          'reactive': sim_settings['start_end'] == 'Reactive'}
        if sim_settings['ACC'] != 'BOTH':
            selection_crit['status'] = sim_settings['ACC']
        if sim_settings['MC-move'] != 'All':
            selection_crit['MC-move'] = sim_settings['MC-move']

        ensemble_names = self.traj_data.keys() if \
            sim_settings['fol'] == 'All' else [sim_settings['fol']]
        for ensemble_name in ensemble_names:
            for cycle in self.traj_data[ensemble_name].keys():
                traj = self.traj_data[ensemble_name][int(cycle)]
                # If columns are missing from the trajectory
                if len(traj.frames.columns) == number_variables:
                    for key in selection_crit:
                        if selection_crit.get(key, 1) != traj.info.get(key, 2):
                            break
                    else:
                        frames.append(traj.frames)
        if frames:
            dataframe = pd.concat(frames, ignore_index=True, sort=False)
        else:
            dataframe = pd.DataFrame()
        dataframe = dataframe.dropna()
        return dataframe

    def load_reactive_unreactive(self, sim_settings):
        """Assign a boolean for reactive (True) and unreactive (False) path.

        Parameters
        ----------
        sim_settings : dict
            Dict containing settings from GUI.

        Returns
        -------
        reactive_data : Pandas Dataframe
            Dataframe containing True/False for each frame in the simulation.
            True if the frame is reactive, False if unreactive.

        """
        reactive_values = []
        number_variables = len(self.infos['op_labels']) + 3

        selection_crit = {'stored': sim_settings['stored'],
                          'reactive': sim_settings['start_end'] == 'Reactive'}
        if sim_settings['ACC'] != 'BOTH':
            selection_crit['status'] = sim_settings['ACC']
        if sim_settings['MC-move'] != 'All':
            selection_crit['MC-move'] = sim_settings['MC-move']

        if sim_settings['fol'] == 'All':
            ensemble_names = self.traj_data.keys()
        else:
            ensemble_names = [sim_settings['fol']]
        for ensemble_name in ensemble_names:
            for cycle in self.traj_data[ensemble_name].keys():
                traj = self.traj_data[ensemble_name][int(cycle)]
                if len(traj.frames.columns) == number_variables:
                    for key in selection_crit:
                        if selection_crit.get(key, 1) != traj.info.get(key, 2):
                            break
                    else:
                        reactive_values.extend(
                            [traj.info['reactive']] * traj.info['length'])
        return pd.DataFrame(reactive_values)

    def remove_values(self, sim_settings, min_op, max_op):
        """Remove values below and above chosen values in a dataframe.

        Parameters
        ----------
        sim_settings : dict
            Dict containing the simulation settings.
        min_op : float
            Lowest chosen value of the order parameter.
        max_op : float
            highest chosen value of the order parameter.

        Returns
        ------
        sim_reduced : Dataframe
            Dataframe containing data from the simulation with removed frames
            where the op is above or below a predetermined value.
        reactive_reduced : Dataframe
             Dataframe containing the reactive/unreactive flag for each
            frame in the simulation with removed frames
            where the op is above or below a predetermined value.

        """
        data = self.load_to_df(sim_settings)
        reactive = pd.DataFrame(np.array(self.load_reactive_unreactive(
            sim_settings), dtype=bool), columns=['reactive'])
        op_label = self.infos['op_labels'][0]
        tot_data = data.join(reactive)
        tot_data = tot_data[~(tot_data[op_label] <= min_op)]
        tot_data = tot_data[~(tot_data[op_label] >= max_op)]
        tot_data = tot_data.reset_index(drop=True)
        reactive_reduced = tot_data['reactive']
        data_reduced = tot_data.drop('reactive', axis=1)
        return data_reduced, reactive_reduced


class VisualObject(DataSlave):
    """DataSlave class used by VisualApp to load data from compressed file."""

    cycle_printed = QtCore.pyqtSignal(list)

    def __init__(self, pfile=None, trajfile=None):
        """Initialize the class and classes inherited.

        Parameters
        ----------
        pfile : string, optional
            The input file name.
        trajfile : string, optional
            The input trajectory file name.

        """
        self.pfile = pfile
        self.trajfile = trajfile
        DataSlave.__init__(self)
        self.settings = {}

    def load_data_from_single(self, iofile):  # pragma: no cover
        """Load a single trajectory file.

        Parameters
        ----------
        iofile : string
            The input file name.

        """
        # If the -data command has been used
        if self.trajfile:
            trajectory_file = self.trajfile
            search_dir = os.path.dirname(
                os.path.abspath(self.trajfile))
        else:
            # Trajectory file selected through the load data option
            trajectory_file = iofile
            search_dir = os.path.dirname(
                os.path.abspath(self.pfile))

        # Find ensemble name and cycle if possible
        ensemble_name, cycle = 'trj', 0
        for part in search_dir.split('/'):
            if part.startswith('0') and len(part) == 3:
                ensemble_name = part
            elif part.isnumeric():
                cycle = int(part)

        # Find the settings from rst files
        if iofile.endswith('.rst'):
            input_settings = settings.parse_settings_file(
                os.path.abspath(iofile))
        else:
            input_settings = settings.parse_settings_file(
                os.path.abspath(find_rst_file(search_dir)))

        # Create a composite order parameter
        functions = create_orderparameter(input_settings)
        frames = {}

        # Calculate collective variables and order parameter
        for idx, order_p in enumerate(
                recalculate_order(functions, trajectory_file, {})):
            frames[idx] = order_p

        # Create data frames
        num_cols = len(frames[0])
        labels = get_cv_names(input_settings, num_columns=num_cols)
        frames = pd.DataFrame.from_dict(frames,
                                        orient='index',
                                        columns=labels)

        # Fill information to info and infos
        info = {'ordermax': np.max(frames[labels[0]]),
                'ordermin': np.min(frames[labels[0]]),
                'length': len(frames.index),
                'ensemble_name': ensemble_name,
                'cycle': cycle}
        # Create Trajectory with frames and info
        traj = Trajectory(frames, info)
        self.traj_data = {ensemble_name: {cycle: traj}}
        infos = {'ensemble_names': [ensemble_name],
                 'cycles': [cycle, cycle],
                 'op_labels': traj.frames.columns,
                 'interfaces': input_settings['simulation']['interfaces']}
        self.infos = infos

    @QtCore.pyqtSlot(str)
    def get_data(self, pfile):  # pragma: no cover
        """Signal emission for loading data from file.

        Parameters
        ----------
        pfile : string
            The input file name.

        """
        self.pfile = pfile
        if self.pfile.endswith(tuple(TRJ_FORMATS)) or \
                (self.trajfile and self.trajfile.endswith(tuple(TRJ_FORMATS))):
            self.load_data_from_single(self.pfile)
            if not self.infos:
                return
        elif pfile.endswith(('.hdf5', '.zip')):
            self.load_whatever()
        else:
            raise ValueError('Format not recognized')
        cycles = self.infos['cycles']
        self.cycle_printed.emit(cycles)


class DataObject(DataSlave, PathDensity):
    """DataSlave class used by VisualApp to load data from PyRETIS input."""

    cycle_printed = QtCore.pyqtSignal(list)

    def __init__(self, iofile=None):  # pragma: no cover
        """Initialize the class and classes inherited.

        Parameters
        ----------
        iofile : string, optional
            The input file name.

        """
        DataSlave.__init__(self, iofile=iofile)
        PathDensity.__init__(self, iofile=iofile)
        self.traj_data = self.traj_dict

    @QtCore.pyqtSlot(str)
    def walk(self):  # pragma: no cover
        """Signal emission to begin 'walk' in the simulation directory."""
        if 'cycles' not in self.infos:
            if 'only_ops' in self.infos.keys():
                self.walk_dirs(self.infos['only_ops'])
            else:
                self.walk_dirs()
        cycles = self.infos['cycles']
        self.cycle_printed.emit(cycles)


def visualize_main(basepath, infile=None, rst_file=None,
                   trajfile=None):  # pragma: no cover
    """Run the VisualApp application.

    Parameters
    ----------
    basepath : string
        The execution folder where the input files are.
    infile : string, optional
        Path to the data that PyVisA is initialized with. Can be a
        PathDensity object, or a hdf5 file.
    rst_file : string, optional
        Path to the rst file if PyVisA has been loaded through the
        input/output rst files of the simulation.
    trajfile: string, optional
        Path to the trajectory file is PyVisA has been loaded with
        the -data command.

    """
    app = QtWidgets.QApplication(sys.argv)
    sub_files = {'rst': rst_file, 'traj': trajfile}
    window = VisualApp(folder=basepath, iofile=infile, sub_files=sub_files)
    window.show()
    sys.exit(app.exec_())

# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""A class for writing text data for the analysis.

This module defines a text plotter which supports the same method as
the generic plotter, however, the output is human-readable text.

Important classes defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TxtPlotter (:py:class:`.TxtPlotter`)
    A class for writing text output.
"""
import gzip
import logging
import os
import shutil
import numpy as np
from pyretis.inout.plotting.plotting import Plotter
from pyretis.inout.common import name_file
from pyretis.inout.common import (ENERFILES, ENERTITLE, FLUXFILES,
                                  ORDERFILES, PATHFILES, PATH_MATCH)
from pyretis.inout.formats.txt_table import txt_save_columns

np.set_printoptions(legacy='1.21')

logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


__all__ = ['TxtPlotter']


class TxtPlotter(Plotter):
    """A plotter writing text-based output.

    This class will just write text-based output. It is similar to
    (:py:class:`.MplPlotter`) in the sense that the same functions are
    supported. Here, however, we do not plot any figures, we just
    write a column based text file with the results.

    Attributes
    ----------
    out_fmt : string
        Selects format for output plots.

    """

    def __init__(self, out_fmt, backup=False, out_dir=None):
        """Initialise the text writer.

        Parameters
        ----------
        out_fmt : string
            This is the format to use for the images.
        backup : boolean, optional
            Determines if we should overwrite or backup old files.
        out_dir : string, optional
            Determines if we should write the files to a particular
            directory.

        """
        super().__init__(backup=backup, plotter_type='text',
                         out_dir=out_dir)
        self.out_fmt = out_fmt

    def output_flux(self, results):
        """Store the output from the flux analysis in text files.

        Parameters
        ----------
        results : dict
            This is the dict with the results from the flux analysis.

        Returns
        -------
        outputfile : list of dicts
            A list containing the files created for the flux and for
            the error in the flux.

        """
        outfiles = []
        # make running average txt and error txt:
        for i in range(len(results['flux'])):
            flux = results['flux'][i]
            runflux = results['runflux'][i]
            errflux = results['errflux'][i]
            outfile = name_file(FLUXFILES['runflux'].format(i + 1),
                                self.out_fmt,
                                path=self.out_dir)
            outfiles.append(outfile)
            # output running average:
            txt_save_columns(outfile, 'Time, running average',
                             (flux[:, 0], runflux), backup=self.backup)
            # output block-error results:
            outfile = name_file(FLUXFILES['block'].format(i + 1),
                                self.out_fmt,
                                path=self.out_dir)
            outfiles.append(outfile)
            _txt_block_error(outfile, 'Block error for flux analysis',
                             errflux, backup=self.backup)
        return outfiles

    def output_energy(self, results, energies):
        """Save the output from the energy analysis to text files.

        Parameters
        ----------
        results : dict
            Each item in `results` contains the results for the
            corresponding energy. It is assumed to contains the keys
            'vpot', 'ekin', 'etot', 'ham', 'temp', 'elec'.
        energies : numpy.array
            This is the raw data for the energy analysis.

        Returns
        -------
        outfiles : list
            The output files created by this function.

        """
        outfiles = []
        time = energies['time']
        # 1) Store the running average:
        header = ['Running average of energy data: time']
        data = [time]
        for key in ['vpot', 'ekin', 'etot', 'ham', 'temp', 'ext']:
            if key in results:
                data.append(results[key]['running'])
                header.append(key)
        headertxt = ' '.join(header)
        outfile = name_file(ENERFILES['run_energies'], self.out_fmt,
                            path=self.out_dir)
        outfiles.append(outfile)
        txt_save_columns(outfile, headertxt, data, backup=self.backup)
        # 2) Save block error data:
        for key in ['vpot', 'ekin', 'etot', 'temp']:
            if key in results:
                outfile = name_file(ENERFILES['block'].format(key),
                                    self.out_fmt, path=self.out_dir)
                outfiles.append(outfile)
                _txt_block_error(outfile, ENERTITLE[key],
                                 results[key]['blockerror'],
                                 backup=self.backup)
        # 3) Save histograms:
        for key in ['vpot', 'ekin', 'etot', 'temp']:
            if key in results:
                outfile = name_file(ENERFILES['dist'].format(key),
                                    self.out_fmt, path=self.out_dir)
                outfiles.append(outfile)
                _txt_histogram(outfile,
                               f'Histogram for {ENERTITLE[key]}',
                               [results[key]['distribution']],
                               backup=self.backup)
        return outfiles

    def output_orderp(self, results, orderdata):
        """Save the output from the order parameter analysis to text files.

        Parameters
        ----------
        results : dict
            Each item in `results` contains the results for the
            corresponding order parameter.
        orderdata : numpy.array
            This is the raw data for the order parameter analysis

        Returns
        -------
        outfiles : list
            The output files created by this function.

        Note
        ----
        We are here only outputting results for the first order parameter.
        I.e. other order parameters or velocities are not written here. This
        will be changed when the structure of the output order parameter
        file has been fixed. Also note that, if present, the first order
        parameter will be plotted against the second one - i.e. the second
        one will be assumed to represent the velocity here.

        """
        outfiles = []
        time = orderdata[:, 0]
        # output running average:
        outfile = name_file(ORDERFILES['run_order'], self.out_fmt,
                            path=self.out_dir)
        txt_save_columns(outfile, 'Time, running average',
                         (time, results[0]['running']),
                         backup=self.backup)
        outfiles.append(outfile)

        # output block-error results:
        outfile = name_file(ORDERFILES['block'], self.out_fmt,
                            path=self.out_dir)
        _txt_block_error(outfile, 'Block error for order param',
                         results[0]['blockerror'], backup=self.backup)
        outfiles.append(outfile)
        # output distributions:
        outfile = name_file(ORDERFILES['dist'], self.out_fmt,
                            path=self.out_dir)
        _txt_histogram(outfile, 'Order parameter',
                       [results[0]['distribution']], backup=self.backup)
        outfiles.append(outfile)
        # output msd if it was calculated:
        if 'msd' in results[0]:
            msd = results[0]['msd']
            outfile = name_file(ORDERFILES['msd'], self.out_fmt,
                                path=self.out_dir)
            txt_save_columns(outfile, 'Time MSD Std',
                             (time[:len(msd)], msd[:, 0], msd[:, 1]),
                             backup=self.backup)
            outfiles.append(outfile)
            # TODO: time c/should here be multiplied with the correct dt
        return outfiles

    def output_path(self, results, path_ensemble):
        """Output all the results obtained by the path analysis.

        Parameters
        ----------
        results : dict
            This dict contains the result from the analysis.
        path_ensemble : object like :py:class:`.PathEnsemble`
            This is the path ensemble we have analysed.

        Returns
        -------
        outfiles : list
            The output files created by this function.

        """
        ens = path_ensemble.ensemble_name  # identify the ensemble
        ens_simplified = path_ensemble.ensemble_name_simple
        detect = results['detect']
        outfiles = []
        if 'pcross' in results:
            # 1) Output pcross vs lambda:
            outfile = name_file(PATHFILES['pcross'].format(ens_simplified),
                                self.out_fmt, path=self.out_dir)
            outfiles.append(outfile)
            txt_save_columns(outfile,
                             f'Ensemble: {ens}, detect: {detect}',
                             [results['pcross'][0], results['pcross'][1]],
                             backup=self.backup)
        if 'prun' in results:
            # 2) Output the running average of p:
            outfile = name_file(PATHFILES['prun'].format(ens_simplified),
                                self.out_fmt, path=self.out_dir)
            outfiles.append(outfile)
            txt_save_columns(outfile, f'Ensemble: {ens}',
                             [results['prun']], backup=self.backup)
        if 'blockerror' in results:
            # 3) Block error results:
            outfile = name_file(PATHFILES['perror'].format(ens_simplified),
                                self.out_fmt, path=self.out_dir)
            outfiles.append(outfile)
            _txt_block_error(outfile, f'Ensemble: {ens}',
                             results['blockerror'], backup=self.backup)
        # 3) Length histograms
        outfile = name_file(PATHFILES['pathlength'].format(ens_simplified),
                            self.out_fmt, path=self.out_dir)
        outfiles.append(outfile)
        _txt_histogram(outfile, 'Histograms for acc and all',
                       [results['pathlength'][0], results['pathlength'][1]],
                       backup=self.backup)
        # 4) Shoot histograms
        outfile = name_file(PATHFILES['shoots'].format(ens_simplified),
                            self.out_fmt, path=self.out_dir)
        outfiles.append(outfile)
        _txt_shoots_histogram(outfile, results['shoots'][0],
                              results['shoots'][1], ens, backup=self.backup)
        return outfiles

    def output_pppath(self, results, path_ensemble):
        """Output all the results obtained by the repptis path analysis.

        Parameters
        ----------
        results : dict
            This dict contains the result from the analysis.
        path_ensemble : object like :py:class:`.PathEnsemble`
            This is the path ensemble we have analysed.

        Returns
        -------
        outfiles : list
            The output files created by this function.

        """
        ens = path_ensemble.ensemble_name  # identify the ensemble
        ens_simplified = path_ensemble.ensemble_name_simple
        detect = results['detect']
        outfiles = []
        if 'pcross' in results:
            # 1) Output pcross vs lambda:
            outfile = name_file(PATHFILES['pcross'].format(ens_simplified),
                                self.out_fmt, path=self.out_dir)
            outfiles.append(outfile)
            txt_save_columns(outfile,
                             f'Ensemble: {ens}, detect: {detect}',
                             [results['pcross'][0], results['pcross'][1]],
                             backup=self.backup)
        if 'prun_sl' in results:
            # 2) Output the running average of p:
            outfile = name_file(PATHFILES['prun_sl'].format(ens_simplified),
                                self.out_fmt, path=self.out_dir)
            outfiles.append(outfile)
            txt_save_columns(outfile, f'Ensemble: {ens}',
                             [results['prun_sl']], backup=self.backup)
        if 'prun_sr' in results:
            # 2) Output the running average of p:
            outfile = name_file(PATHFILES['prun_sr'].format(ens_simplified),
                                self.out_fmt, path=self.out_dir)
            outfiles.append(outfile)
            txt_save_columns(outfile, f'Ensemble: {ens}',
                             [results['prun_sr']], backup=self.backup)
        if 'blockerror_sl' in results:
            # 3) Block error results:
            outfile = name_file(PATHFILES['perror_sl'].format(ens_simplified),
                                self.out_fmt, path=self.out_dir)
            outfiles.append(outfile)
            _txt_block_error(outfile, f'Ensemble: {ens}',
                             results['blockerror_sl'], backup=self.backup)
        if 'blockerror_sr' in results:
            # 3) Block error results:
            outfile = name_file(PATHFILES['perror_sr'].format(ens_simplified),
                                self.out_fmt, path=self.out_dir)
            outfiles.append(outfile)
            _txt_block_error(outfile, f'Ensemble: {ens}',
                             results['blockerror_sr'], backup=self.backup)
        # 3) Length histograms
        outfile = name_file(PATHFILES['pathlength'].format(ens_simplified),
                            self.out_fmt, path=self.out_dir)
        outfiles.append(outfile)
        _txt_histogram(outfile, 'Histograms for acc and all',
                       [results['pathlength'][0], results['pathlength'][1]],
                       backup=self.backup)
        # 4) Shoot histograms
        outfile = name_file(PATHFILES['shoots'].format(ens_simplified),
                            self.out_fmt, path=self.out_dir)
        outfiles.append(outfile)
        _txt_shoots_histogram(outfile, results['shoots'][0],
                              results['shoots'][1], ens, backup=self.backup)
        return outfiles

    def output_matched_probability(self, path_ensembles, detect, matched,
                                   reptis=False):
        """Output the matched probabilities to a text file.

        This function will output the matched probabilities for the
        different ensembles and also output the over-all matched
        probability.

        Parameters
        ----------
        path_ensembles : list of strings
            This is the names for the path ensembles we have calculated
            the probability for.
        detect : list of floats
            These are the detect interfaces used in the analysis.
        matched : dict
            This dict contains the results from the matching of the
            probabilities. We make use of `matched['overall-prob']` and
            `matched['matched-prob']` here.

        Returns
        -------
        outfiles : list
            The files created by this function.

        """
        outfiles = []
        # start by creating the matched file, here we use a custom
        # file writer:
        if not reptis:
            outfile = name_file(PATH_MATCH['match'], self.out_fmt,
                                path=self.out_dir)
            use_gzip = outfile.endswith('.gz')
            if use_gzip:
                outfile = outfile[:-3]
            with open(outfile, 'wb') as fhandle:
                for prob, ens, idet in zip(matched['matched-prob'],
                                           path_ensembles, detect):
                    header = f'Ensemble: {ens}, idetect: {idet}'
                    np.savetxt(fhandle, prob, header=header)
            if use_gzip:
                with open(outfile, 'rb') as fhandle:
                    with gzip.open(f'{outfile}.gz', 'wb') as fhandle_out:
                        shutil.copyfileobj(fhandle, fhandle_out)
                os.remove(outfile)
                outfile = f'{outfile}.gz'
            outfiles.append(outfile)
        # output the over-all matched probability:
        if reptis:
            matched['overall-prob'] = np.array([detect,
                                                matched['pcross']]).T
        outfile = name_file(PATH_MATCH['total'], self.out_fmt,
                            path=self.out_dir)
        interf = ' , '.join([str(idet) for idet in detect])
        header = 'Total matched probability. Interfaces: {}'
        txt_save_columns(outfile, header.format(interf),
                         (matched['overall-prob'][:, 0],
                          matched['overall-prob'][:, 1]),
                         backup=self.backup)
        outfiles.append(outfile)

        # output the prun matched probability:
        if 'overall-rrun' in matched:
            outfile = name_file('overall-rrun', self.out_fmt,
                                path=self.out_dir)
            interf = ' , '.join([str(idet) for idet in detect])
            header = 'Running average overall rate. Interfaces: {}'
            txt_save_columns(outfile, header.format(interf),
                             (matched['overall-cycle'],
                              matched['overall-rrun']),
                             backup=self.backup)
            outfiles.append(outfile)
        if 'overall-prun' in matched:
            outfile = name_file('overall-prun', self.out_fmt,
                                path=self.out_dir)
            interf = ' , '.join([str(idet) for idet in detect])
            header = 'Running average crossing probability. Interfaces: {}'
            txt_save_columns(outfile, header.format(interf),
                             (matched['overall-cycle'],
                              matched['overall-prun']),
                             backup=self.backup)
            outfiles.append(outfile)

        return outfiles


def _txt_block_error(outputfile, title, error, backup=False):
    """Write the output from the error analysis to a text file.

    Parameters
    ----------
    outputfile : string
        This is the name of the output file to create.
    title : string
        This is an identifier/title to add to the header, e.g.
        'Ensemble: 001', 'Kinetic energy', etc.
    error : list
        This is the result of the error analysis.
    backup : boolean, optional
        Determines if we will do a backup of old files or not.

    """
    header = '{0}, Rel.err: {1:9.6e}, Ncor: {2:9.6f}'
    header = header.format(title, error[4], error[6])
    txt_save_columns(outputfile, header, (error[0], error[3]), backup=backup)


def _txt_histogram(outputfile, title, histograms, backup=False):
    """Write histograms to a text file.

    Parameters
    ----------
    outputfile : string
        This is the name of the output file to create.
    title : string
        A descriptive title to add to the header.
    histograms : tuple or list
        The histograms to store.
    backup : boolean, optional
        Determines if we will do a backup of old files or not.

    """
    data = []
    header = [f'{title}']
    for hist in histograms:
        header.append(f'avg: {hist[2][0]:6.2f}, std: {hist[2][1]:6.2f}')
        data.append(hist[1])
        data.append(hist[0])
    headertxt = ', '.join(header)
    txt_save_columns(outputfile, headertxt, data, backup=backup)


def _txt_shoots_histogram(outputfile, histograms, scale, ensemble,
                          backup=False):
    """Write the histograms from the shoots analysis to a text file.

    Parameters
    ----------
    histograms : dict
        These are the histograms obtained in the shoots analysis.
    scale : dict
        These are the scale factors for normalising the histograms
        obtained in the shoots analysis.
    ensemble : string
        This is the ensemble identifier, e.g. 001, 002, etc.
    outputfile : string
        This is the name of the output file to create.
    backup : boolean, optional
        Determines if we will do a backup of old files or not.

    """
    data = []
    header = [f'Ensemble: {ensemble}']
    for key in ['ACC', 'REJ', 'BWI', 'ALL']:
        try:
            mid = histograms[key][2]
            hist = histograms[key][0]
            hist_scale = hist * scale[key]
            data.append(mid)
            data.append(hist)
            data.append(hist_scale)
            header.append(f'{key} (mid, hist, hist*scale)')
        except KeyError:
            continue
    headertxt = ', '.join(header)
    txt_save_columns(outputfile, headertxt, data, backup=backup)

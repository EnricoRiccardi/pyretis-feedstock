# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Module defining functions useful in the analysis of simulation data.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

running_average (:py:func:`.running_average`)
    Method to calculate a running average.

block_error (:py:func:`.block_error`)
    Perform block error analysis.

block_error_corr (:py:func:`.block_error_corr`)
    Method to run a block error analysis and calculate relative
    errors and correlation length.
"""
import numpy as np
from pyretis.analysis.histogram import histogram_and_avg

__all__ = ['running_average', 'block_error', 'block_error_corr',
           'analyse_data', 'mean_square_displacement']

np.seterr(divide='ignore', invalid='ignore')
np.set_printoptions(legacy='1.21')


def running_average(data):
    """Create a running average of the given data.

    The running average will be calculated over the rows.

    Parameters
    ----------
    data : numpy.array
        This is the data we will average.

    Returns
    -------
    out : numpy.array
        The running average.

    """
    if len(data.shape) == 2:
        w_data = data[:, 0]*data[:, 1]
        return w_data.cumsum() / data[:, 1].cumsum()
    return data.cumsum() / np.ones(data.shape[0]).cumsum()


def block_error(data, maxblock=None, blockskip=1, weights=None):
    """Perform block error analysis.

    This function will estimate the standard deviation in the input
    data by performing a block analysis. The number of blocks
    to consider can be specified or it will be taken as the
    half of the length of the input data. Averages and variance are
    calculated using an on-the-fly algorithm [1]_.

    Parameters
    ----------
    data : numpy.array (or iterable with data points)
        The data to analyse.
    maxblock : int, optional
        Can be used to set the maximum length of the blocks to
        consider. Note that the `maxblock` will never be set longer
        than half the length in data.
    blockskip : int, optional
        This can be used to skip certain block lengths, i.e.
        `blockskip = 1` will consider all blocks up to `maxblock`, while
        `blockskip = n` will consider every n'th block up to `maxblock`,
        i.e. it will use block lengths equal to `1`, `1 + n`, `1 + 2*n`,
        and so on.

    Returns
    -------
    blocklen : numpy.array
        These contain the block lengths considered.
    block_avg : numpy.array
        The averages as a function of the block length.
    block_err : numpy.array
        Estimate of errors as a function of the block length.
    block_err_avg : float
        Average of the error estimate using blocks where
        ``length > maxblock//2``.

    References
    ----------
    .. [1] Wikipedia, "Algorithms for calculating variance",
       http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance

    """
    if maxblock is None or maxblock < 1:
        maxblock = len(data) // 2
    else:
        maxblock = min(maxblock, len(data) // 2)
    indata = data
    if weights is None:
        if len(data.shape) == 1:
            weights = np.ones(len(data))*2
        else:
            data = (i[0] for i in indata)
            weights = (i[1] for i in indata)
    # define helper variables:
    blocklen = np.arange(0, maxblock, blockskip, dtype='i4')
    # blocklen contains the lengths of the blocks
    blocklen += 1
    # +1 to make blocklen[i] = length of block no i where numbering
    # starts at 0 -> blocklen[0] = 1 and so on. Note that arange does
    # create [0, ..., maxblock).
    block = np.zeros(len(blocklen))  # to accumulate values for a block
    blockw = np.zeros(len(blocklen))  # accumulate weights for a block
    tot_w = np.zeros(block.shape)  # accumulate total weights
    tot_w_s = np.zeros(block.shape)  # accumulate total weights**2
    block_avg = np.zeros(block.shape)  # to store averages in block
    block_var = np.zeros(block.shape)  # estimator of variance

    # Algorithm taken from
    # https://markusthill.github.io/math/stats/ml/
    # online-estimation-of-weighted-sample-mean-and-coviarance-matrix/
    # equation 39 where we delay the normalization to the end
    for i, (datai, weighti) in enumerate(zip(data, weights)):
        block += weighti*datai  # accumulate the value to all blocks
        blockw += weighti
        # next pick out blocks which are "full":
        k = np.where((i + 1) % blocklen == 0)[0]
        # update estimate of average and variance
        block[k] = np.nan_to_num(block[k] / blockw[k])  # catch 0/0 = NaN
        tot_w[k] += blockw[k]
        tot_w_s[k] += blockw[k]**2
        delta1 = blockw[k] * (block[k] - block_avg[k])
        block_avg[k] = block_avg[k] + np.nan_to_num(delta1 / tot_w[k])
        block_var[k] = block_var[k] + delta1 * (block[k] - block_avg[k])
        # reset these blocks
        block[k] = 0.0
        blockw[k] = 0

    nsamp = tot_w * tot_w / tot_w_s
    block_var *= nsamp / (nsamp - 1) / tot_w
    block_err = np.sqrt(block_var / nsamp)  # estimate of error
    block_err_avg = np.average(
        block_err[np.where(blocklen > maxblock // 2)[0]])
    return blocklen, block_avg, block_err, block_err_avg


def block_error_corr(data, maxblock=None, blockskip=1):
    """Run block error analysis on the given data.

    This will run the block error analysis and return the relative
    errors and correlation length.

    Parameters
    ----------
    data : numpy.array
        Data to analyse.
    maxblock : int, optional
        The maximum block length to consider.
    blockskip : int, optional
        This can be used to skip certain block lengths, i.e.
        `blockskip = 1` will consider all blocks up to `maxblock`, while
        `blockskip = n` will consider every n'th block up to `maxblock`,
        i.e. it will use block lengths equal to `1`, `1 + n`, `1 + 2*n`,
        and so on.

    Returns
    -------
    out[0] : numpy.array
        These contains the block lengths considered (`blen`).
    out[1] : numpy.array
        Estimate of errors as a function of the block length (`berr`).
    out[2] : float
        Average of the error estimate for blocks (`berr_avg`)
        with ``length > maxblock // 2``.
    out[3] : numpy.array
        Estimate of relative errors normalised by the overall average
        as a function of block length (`rel_err`).
    out[4] : float
        The average relative error (`avg_rel_err`), for blocks
        with ``length > maxblock // 2``.
    out[5] : numpy.array
        The estimated correlation length as a function of the block
        length (`ncor`).
    out[6] : float
        The average (for blocks with length > maxblock // 2) estimated
        correlation length (`avg_ncor`).

    """
    blen, bavg, berr, berr_avg = block_error(data, maxblock=maxblock,
                                             blockskip=blockskip)
    if len(bavg) == 0 or abs(bavg[0]) == 0:
        nan = np.full_like(berr, np.nan)
        return blen, berr, berr_avg, nan, np.nan, nan, np.nan
    # also calculate some relative errors:
    rel_err = np.divide(berr, abs(bavg[0]))
    avg_rel_err = np.divide(berr_avg, abs(bavg[0]))
    ncor = np.divide(berr**2, berr[0]**2)
    avg_ncor = np.divide(berr_avg**2, berr[0]**2)
    return blen, berr, berr_avg, rel_err, avg_rel_err, ncor, avg_ncor


def mean_square_displacement(data, ndt=None):
    """Calculate the mean square displacement for the given data.

    Parameters
    ----------
    data : numpy.array, 1D
        This numpy.array contains the data as a function of time.
    ndt : int, optional
        This parameter is the number of time origins. I.e. points up to
        ndt will be used as time origins. If not specified the value of
        the input ``data.size // 5`` will be used.

    Returns
    -------
    msd : numpy.array, 2D
        The first column is the mean squared displacement and the
        second column is the corresponding standard deviation.

    """
    length = data.shape[0]
    if ndt is None or ndt < 1:
        ndt = length // 5
    msd = []
    for i in range(1, ndt):
        delta = (data[i:] - data[:-i])**2
        msd.append((delta.mean(), delta.std()))
    return np.array(msd)


def analyse_data(data, settings):
    """Analyse the given data and run some common analysis procedures.

    Specifically, it will:

    1) Calculate a running average.

    2) Obtain a histogram.

    3) Run a block error analysis.


    Parameters
    ----------
    data : numpy.array, 1D
        This numpy.array contains the data as a function of time.
    settings : dict
        This dictionary contains settings for the analysis.

    Returns
    -------
    result : dict
        This dict contains the results.

    """
    result = {}
    asett = settings['analysis']
    # 1) Do the running average
    result['running'] = running_average(data)
    # 2) Obtain distributions:
    result['distribution'] = histogram_and_avg(data, asett['bins'],
                                               density=True)
    # 3) Do the block error analysis:
    result['blockerror'] = block_error_corr(data,
                                            maxblock=asett['maxblock'],
                                            blockskip=asett['blockskip'])
    return result

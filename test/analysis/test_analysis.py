# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test block error analysis."""
import pytest
import numpy as np
from pyretis.analysis import block_error
from pyretis.inout.analysisio.analysisio import analyse_file


class TestAnalysis:
    """Test that we can do block error analysis."""

    def test_block_error_weights(self):
        data = np.random.random(100)
        weights = np.ones(data.shape)*(np.random.random() + 0.001)  # catch 0.0
        r1 = block_error(data)
        r2 = block_error(data, weights=weights)
        r3 = block_error(np.stack([data, weights], axis=-1))
        for i, j, k in zip(r1, r2, r3):
            assert np.allclose(i, j)
            assert np.allclose(i, k)

    def test_block_error_zero_weights(self):
        # This has to be even to not alter the number of blocks
        data = np.random.random(100)
        weights = np.ones(data.shape)

        # add a giant data point that should not matter
        data2 = np.array([i for i in data]+[1e10])
        weights2 = np.array([i for i in weights]+[0])
        r1 = block_error(data, weights=weights)
        r2 = block_error(data2, weights=weights2)
        for i, j in zip(r1, r2):
            assert np.allclose(i, j)

    def test_block_error_zero_weights_start(self):
        data = np.random.random(100)
        weights = np.ones(data.shape)
        weights[0] = 0  # this would lead to nan's
        r1 = block_error(data, weights=weights)
        assert not np.isnan(r1[1][0])

    def test_analyse_file_fail(self):
        """Test the function to fail."""
        with pytest.raises(ValueError) as err:
            analyse_file('no_ext', 'non_existent', {})
        assert 'Unknown file type' in str(err.value)

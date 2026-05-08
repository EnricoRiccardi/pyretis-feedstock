# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Tests for pyretis.testing.simulation_comparison (coverage completion).

These tests cover functions and branches not exercised by test_compare.py:
- compare_numerical_mse
- compare_restarted_text_files (all branches)
- compare_reports_normalized (all branches)
- compare_restarted_cross_files (all branches)
- skip_keys branch in compare_text_line_by_line
- token-length mismatch branch in _compare_block_comments
"""
import os
import numpy as np
import pytest

from pyretis.testing.simulation_comparison import (
    compare_numerical_mse,
    compare_restarted_text_files,
    compare_reports_normalized,
    compare_restarted_cross_files,
    compare_text_line_by_line,
    _compare_block_comments,
    read_files,
)

HERE = os.path.abspath(os.path.dirname(__file__))
TEST_PATH = os.path.join(HERE, 'files_for_comparison')
ENERGY1 = os.path.join(TEST_PATH, 'energy1.txt')
ENERGY2 = os.path.join(TEST_PATH, 'energy2.txt')
ENERGY5 = os.path.join(TEST_PATH, 'energy5.txt')


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write(path, content):
    """Write text content to a file."""
    path.write_text(content, encoding='utf-8')
    return str(path)


# Cross-file header + data lines used in multiple tests
_CROSS_HEADER = '#     Step  Int Dir\n'
_CROSS_ROW1 = '       592    1   -\n'
_CROSS_ROW2 = '      1171    1   +\n'
_CROSS_ROW3 = '      1642    1   -\n'


# ---------------------------------------------------------------------------
# compare_numerical_mse
# ---------------------------------------------------------------------------

class TestCompareNumericalMse:
    """Tests for compare_numerical_mse."""

    def test_equal_files(self, tmp_path):
        """MSE of a file compared with itself is 0."""
        np.savetxt(str(tmp_path / 'a.txt'), np.array([[1.0, 2.0], [3.0, 4.0]]))
        f = str(tmp_path / 'a.txt')
        ok, msg = compare_numerical_mse(f, f)
        assert ok
        assert 'within tolerance' in msg

    def test_different_values(self, tmp_path):
        """MSE above tolerance returns False."""
        np.savetxt(str(tmp_path / 'a.txt'), np.array([[1.0, 2.0]]))
        np.savetxt(str(tmp_path / 'b.txt'), np.array([[1.0, 5.0]]))
        ok, msg = compare_numerical_mse(str(tmp_path / 'a.txt'),
                                        str(tmp_path / 'b.txt'), tol=1e-12)
        assert not ok
        assert 'MSE' in msg

    def test_shape_mismatch(self, tmp_path):
        """Differing shapes return False with shape message."""
        np.savetxt(str(tmp_path / 'a.txt'), np.array([[1.0, 2.0]]))
        np.savetxt(str(tmp_path / 'b.txt'), np.array([[1.0, 2.0], [3.0, 4.0]]))
        ok, msg = compare_numerical_mse(str(tmp_path / 'a.txt'),
                                        str(tmp_path / 'b.txt'))
        assert not ok
        assert 'Shapes differ' in msg


# ---------------------------------------------------------------------------
# compare_restarted_text_files
# ---------------------------------------------------------------------------

class TestCompareRestartedTextFiles:
    """Tests for compare_restarted_text_files (all branches)."""

    def test_matching_restart(self, tmp_path):
        """file11 + file12 (minus overlap) equals file2."""
        part1 = _write(tmp_path / 'p1.txt',
                       '# header\nline A\nline B\n')
        # part2 starts with a header, then the overlap (last line of p1),
        # then new data
        part2 = _write(tmp_path / 'p2.txt',
                       '# header2\nline B\nline C\n')
        full = _write(tmp_path / 'full.txt',
                      '# header\nline A\nline B\nline C\n')
        ok, msg = compare_restarted_text_files(part1, part2, full)
        assert ok
        assert 'match' in msg

    def test_part2_no_data(self, tmp_path):
        """Part2 with only header lines → False."""
        part1 = _write(tmp_path / 'p1.txt', '# header\nline A\n')
        part2 = _write(tmp_path / 'p2.txt', '# header only\n')
        full = _write(tmp_path / 'full.txt', '# header\nline A\n')
        ok, msg = compare_restarted_text_files(part1, part2, full)
        assert not ok
        assert 'no data' in msg.lower()

    def test_overlap_mismatch(self, tmp_path):
        """Last line of part1 ≠ first data line of part2 → False."""
        part1 = _write(tmp_path / 'p1.txt', '# header\nline A\nline B\n')
        part2 = _write(tmp_path / 'p2.txt', '# header\nline X\nline C\n')
        full = _write(tmp_path / 'full.txt',
                      '# header\nline A\nline B\nline C\n')
        ok, msg = compare_restarted_text_files(part1, part2, full)
        assert not ok
        assert 'Overlapping' in msg

    def test_line_count_mismatch(self, tmp_path):
        """Combined line count ≠ full file line count → False."""
        part1 = _write(tmp_path / 'p1.txt', '# header\nline A\nline B\n')
        part2 = _write(tmp_path / 'p2.txt', '# header\nline B\nline C\n')
        # Full has an extra line
        full = _write(tmp_path / 'full.txt',
                      '# header\nline A\nline B\nline C\nextra\n')
        ok, msg = compare_restarted_text_files(part1, part2, full)
        assert not ok
        assert 'mismatch' in msg.lower()

    def test_content_mismatch(self, tmp_path):
        """Matching line count but differing content → False."""
        part1 = _write(tmp_path / 'p1.txt', '# header\nline A\nline B\n')
        part2 = _write(tmp_path / 'p2.txt', '# header\nline B\nline C\n')
        full = _write(tmp_path / 'full.txt',
                      '# header\nline A\nline B\nline DIFFERENT\n')
        ok, msg = compare_restarted_text_files(part1, part2, full)
        assert not ok
        assert 'Mismatch' in msg


# ---------------------------------------------------------------------------
# compare_reports_normalized
# ---------------------------------------------------------------------------

class TestCompareReportsNormalized:
    """Tests for compare_reports_normalized (all branches)."""

    def test_equal_reports(self, tmp_path):
        """Identical reports are equal."""
        content = '<html><body><p>Hello</p></body></html>\n'
        f1 = _write(tmp_path / 'r1.html', content)
        f2 = _write(tmp_path / 'r2.html', content)
        ok, msg = compare_reports_normalized(f1, f2)
        assert ok
        assert 'equal' in msg.lower()

    def test_generator_meta_skipped(self, tmp_path):
        """Lines with 'meta name="generator"' are ignored."""
        body = '<html>\n<p>text</p>\n</html>\n'
        f1 = _write(tmp_path / 'r1.html',
                    '<meta name="generator" content="Docutils 0.17">\n' + body)
        f2 = _write(tmp_path / 'r2.html',
                    '<meta name="generator" content="Docutils 0.18">\n' + body)
        ok, msg = compare_reports_normalized(f1, f2)
        assert ok

    def test_timestamp_skipped(self, tmp_path):
        """Lines with 'generated by PyRETIS' are ignored."""
        body = '<html>\n<p>text</p>\n</html>\n'
        f1 = _write(tmp_path / 'r1.html',
                    body + 'generated by PyRETIS on 2024-01-01\n')
        f2 = _write(tmp_path / 'r2.html',
                    body + 'generated by PyRETIS on 2025-06-15\n')
        ok, msg = compare_reports_normalized(f1, f2)
        assert ok

    def test_gray_grey_normalized(self, tmp_path):
        """'color: gray' and 'color: grey' are treated as equal."""
        f1 = _write(tmp_path / 'r1.css', 'a { color: gray; }\n')
        f2 = _write(tmp_path / 'r2.css', 'a { color: grey; }\n')
        ok, msg = compare_reports_normalized(f1, f2)
        assert ok

    def test_different_content(self, tmp_path):
        """Reports with genuinely different content are not equal."""
        f1 = _write(tmp_path / 'r1.html', '<p>Hello</p>\n')
        f2 = _write(tmp_path / 'r2.html', '<p>World</p>\n')
        ok, msg = compare_reports_normalized(f1, f2)
        assert not ok
        assert 'Mismatch' in msg

    def test_different_line_count(self, tmp_path):
        """Reports with different numbers of content lines are not equal."""
        f1 = _write(tmp_path / 'r1.html', '<p>A</p>\n<p>B</p>\n')
        f2 = _write(tmp_path / 'r2.html', '<p>A</p>\n')
        ok, msg = compare_reports_normalized(f1, f2)
        assert not ok
        assert 'number' in msg.lower()

    def test_id_line_skipped(self, tmp_path):
        """':Id: $Id: html4css1.css' lines are ignored."""
        body = '<p>content</p>\n'
        f1 = _write(tmp_path / 'r1.html',
                    body + ':Id: $Id: html4css1.css v1.0\n')
        f2 = _write(tmp_path / 'r2.html',
                    body + ':Id: $Id: html4css1.css v2.0\n')
        ok, _ = compare_reports_normalized(f1, f2)
        assert ok


# ---------------------------------------------------------------------------
# compare_restarted_cross_files
# ---------------------------------------------------------------------------

class TestCompareRestartedCrossFiles:
    """Tests for compare_restarted_cross_files."""

    def _make_cross(self, tmp_path, name, rows):
        """Write a cross file with the given data rows."""
        content = _CROSS_HEADER + ''.join(rows)
        return _write(tmp_path / name, content)

    def test_matching_with_overlap(self, tmp_path):
        """part1 + part2 (deduplicating shared last/first row) == full."""
        f11 = self._make_cross(tmp_path, 'c11.txt',
                               [_CROSS_ROW1, _CROSS_ROW2])
        f12 = self._make_cross(tmp_path, 'c12.txt',
                               [_CROSS_ROW2, _CROSS_ROW3])
        f2 = self._make_cross(tmp_path, 'c2.txt',
                              [_CROSS_ROW1, _CROSS_ROW2, _CROSS_ROW3])
        ok, msg = compare_restarted_cross_files(f11, f12, f2)
        assert ok
        assert 'matches' in msg

    def test_matching_without_overlap(self, tmp_path):
        """When last row of part1 differs from first row of part2."""
        f11 = self._make_cross(tmp_path, 'c11.txt', [_CROSS_ROW1])
        f12 = self._make_cross(tmp_path, 'c12.txt', [_CROSS_ROW2])
        f2 = self._make_cross(tmp_path, 'c2.txt',
                              [_CROSS_ROW1, _CROSS_ROW2])
        ok, msg = compare_restarted_cross_files(f11, f12, f2)
        assert ok

    def test_shape_mismatch(self, tmp_path):
        """Combined data shape ≠ full data shape returns False."""
        f11 = self._make_cross(tmp_path, 'c11.txt', [_CROSS_ROW1])
        f12 = self._make_cross(tmp_path, 'c12.txt', [_CROSS_ROW2, _CROSS_ROW3])
        # Full has only one row
        f2 = self._make_cross(tmp_path, 'c2.txt', [_CROSS_ROW1])
        ok, msg = compare_restarted_cross_files(f11, f12, f2)
        assert not ok
        assert 'shape' in msg.lower() or 'mismatch' in msg.lower()

    def test_data_mismatch(self, tmp_path):
        """Same shape but different data returns False."""
        f11 = self._make_cross(tmp_path, 'c11.txt', [_CROSS_ROW1])
        f12 = self._make_cross(tmp_path, 'c12.txt', [_CROSS_ROW2])
        # Full has different content
        f2 = self._make_cross(tmp_path, 'c2.txt',
                              [_CROSS_ROW1, _CROSS_ROW3])
        ok, msg = compare_restarted_cross_files(f11, f12, f2)
        assert not ok
        assert 'mismatch' in msg.lower()

    def test_empty_part1(self, tmp_path):
        """Empty part1: combined = part2."""
        f11 = self._make_cross(tmp_path, 'c11.txt', [])
        f12 = self._make_cross(tmp_path, 'c12.txt', [_CROSS_ROW1, _CROSS_ROW2])
        f2 = self._make_cross(tmp_path, 'c2.txt', [_CROSS_ROW1, _CROSS_ROW2])
        ok, msg = compare_restarted_cross_files(f11, f12, f2)
        assert ok


# ---------------------------------------------------------------------------
# Branches in existing functions not yet covered
# ---------------------------------------------------------------------------

class TestUncoveredBranches:
    """Cover the remaining branches in already-tested functions."""

    def test_compare_text_skip_keys(self, tmp_path):
        """skip_keys filters out matching lines before comparison."""
        f1 = _write(tmp_path / 'a.txt',
                    'exe_path = /old/path\nvalue = 1\n')
        f2 = _write(tmp_path / 'b.txt',
                    'exe_path = /new/path\nvalue = 1\n')
        ok, msg = compare_text_line_by_line(f1, f2, skip_keys=['exe_path'])
        assert ok

    def test_compare_text_skip_keys_line_count_differs(self, tmp_path):
        """skip_keys can cause line count to differ → not equal."""
        f1 = _write(tmp_path / 'a.txt',
                    'exe_path = /old\nvalue = 1\nextra = x\n')
        f2 = _write(tmp_path / 'b.txt',
                    'exe_path = /new\nvalue = 1\n')
        ok, msg = compare_text_line_by_line(f1, f2, skip_keys=['exe_path'])
        assert not ok

    def test_compare_block_comments_token_length_mismatch(self):
        """_compare_block_comments: different token count per line → False."""
        c1 = ['# Cycle: 0, status: ACC']
        c2 = ['# Cycle: 0']   # fewer tokens
        assert not _compare_block_comments(c1, c2)

    def test_compare_block_comments_nonfloat_token_differs(self):
        """_compare_block_comments: non-float differing token → False."""
        c1 = ['# status: ACC extra']
        c2 = ['# status: REJ extra']
        assert not _compare_block_comments(c1, c2)

    def test_read_files_skip_comments(self, tmp_path):
        """read_files with read_comments=False skips '#' lines."""
        f = _write(tmp_path / 'f.txt', '# comment\ndata line\n')
        result = read_files(str(f), read_comments=False)
        assert result == [['data line\n']]

    def test_compare_block_comments_different_lengths(self):
        """_compare_block_comments: different number of lines → False."""
        c1 = ['# line1', '# line2']
        c2 = ['# line1']
        assert not _compare_block_comments(c1, c2)

    def test_compare_block_comments_first_line_equal_second_differs(self):
        """_compare_block_comments: first line equal, second differs."""
        c1 = ['# same header', '# status: ACC']
        c2 = ['# same header', '# status: REJ']
        assert not _compare_block_comments(c1, c2)

    def test_compare_block_comments_equal_inner_token(self):
        """_compare_block_comments: some tokens equal, floats too far."""
        # '#' and 'val:' tokens match (inner continue), '1.0' vs '2.0' differ
        c1 = ['# val: 1.0']
        c2 = ['# val: 2.0']
        assert not _compare_block_comments(c1, c2)

    def test_compare_block_comments_close_floats_returns_true(self):
        """_compare_block_comments: '1.0' vs '1' are same float, True."""
        # Lines differ as strings but '1.0' and '1' are equal as floats,
        # so math.isclose returns True and the function returns True.
        c1 = ['# val: 1.0']
        c2 = ['# val: 1']
        assert _compare_block_comments(c1, c2)

    def test_compare_restarted_cross_files_empty_part2(self, tmp_path):
        """compare_restarted_cross_files: non-empty part1, empty part2."""
        f11 = _write(tmp_path / 'c11.txt',
                     _CROSS_HEADER + _CROSS_ROW1 + _CROSS_ROW2)
        # part2 has only header (no data rows)
        f12 = _write(tmp_path / 'c12.txt', _CROSS_HEADER)
        f2 = _write(tmp_path / 'c2.txt',
                    _CROSS_HEADER + _CROSS_ROW1 + _CROSS_ROW2)
        ok, msg = compare_restarted_cross_files(f11, f12, f2)
        assert ok

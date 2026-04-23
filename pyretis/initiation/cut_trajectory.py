# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Utility functions for cutting GROMACS and XYZ trajectories.

This module was created as part of the PyRETIS trajectory loading
parallelization.
"""
import os
import re
import shutil
import subprocess


def get_trajectory_timestep(filename, gmx_exe='gmx'):
    """Get the time step between the first two frames of a .trr/.xtc file.

    Uses ``gmx dump`` and kills the process as soon as two timestamps
    have been parsed, so only a few KB of I/O are needed.

    Parameters
    ----------
    filename : str
        Path to the trajectory file.
    gmx_exe : str, optional
        Path to the GROMACS executable.

    Returns
    -------
    dt : float or None
        Time step in ps, or None if it cannot be determined.
    """
    try:
        proc = subprocess.Popen(
            [gmx_exe, 'dump', '-f', filename],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        times = _read_two_timestamps(proc)
        proc.kill()
        proc.wait()
        if len(times) >= 2:
            dt = times[1] - times[0]
            if dt > 0:
                return dt
    except (OSError, subprocess.SubprocessError):
        pass
    return None


def _read_two_timestamps(proc):
    """Read up to two timestamps from a running gmx dump process.

    Parameters
    ----------
    proc : subprocess.Popen
        Running gmx dump process with stdout pipe.

    Returns
    -------
    times : list of float
        Up to two timestamp values parsed from the output.
    """
    times = []
    for raw_line in proc.stdout:
        line = raw_line.decode('utf-8', errors='replace')
        if '   t=' in line:
            parts = line.split('t=')
            if len(parts) >= 2:
                try:
                    t_val = float(parts[-1].strip().split()[0])
                    times.append(t_val)
                except (ValueError, IndexError):
                    pass
            if len(times) >= 2:
                break
    return times


def do_cut_trajectory_split(filename, ext, frame_indices, chunk_frames,
                            tmpdir, engine_info):
    """Cut a trajectory into chunks with a single ``gmx trjconv -split``.

    A single trjconv call reads the source file once and writes all
    chunks in a single pass, avoiding repeated full-file reads.

    Parameters
    ----------
    filename : str
        Path to the source trajectory file.
    ext : str
        File extension ('.trr' or '.xtc').
    frame_indices : list of int
        All 0-based frame indices to extract.
    chunk_frames : int
        Desired number of frames per output chunk.
    tmpdir : str
        Temporary directory for the output files.
    engine_info : dict
        Must contain 'gmx' key with path to the GROMACS executable.

    Returns
    -------
    list of str or None
        Sorted output file paths, or None on failure (caller should
        fall back to the serial cutting approach).
    """
    if ext not in ['.trr', '.xtc']:
        return None

    gmx_exe = engine_info.get('gmx', 'gmx')

    dt = get_trajectory_timestep(filename, gmx_exe)
    if dt is None or dt <= 0:
        return None

    ndx_file = os.path.join(tmpdir, 'all_frames.ndx')
    with open(ndx_file, 'w', encoding='utf-8') as ndxf:
        ndxf.write('[ frames ]\n')
        for idx in sorted(frame_indices):
            ndxf.write(f'{idx + 1}\n')  # 1-based for GROMACS

    split_time = chunk_frames * dt
    output_base = os.path.join(tmpdir, f'seg{ext}')

    subprocess.run(
        [gmx_exe, 'trjconv',
         '-f', filename,
         '-o', output_base,
         '-fr', ndx_file,
         '-split', str(split_time)],
        input=b'0\n',
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False
    )

    output_files = [
        os.path.join(tmpdir, f)
        for f in os.listdir(tmpdir)
        if f.startswith('seg') and f.endswith(ext)
    ]
    if not output_files:
        return None

    def _sort_key(path):
        nums = re.findall(r'(\d+)', os.path.basename(path).replace(ext, ''))
        return int(nums[-1]) if nums else 0

    output_files.sort(key=_sort_key)
    return output_files


def _cut_gromacs(filename, ext, minidx, maxidx, cut_filename, gmx_exe):
    """Cut a GROMACS trajectory segment.

    Parameters
    ----------
    filename : str
        Path to the source trajectory file.
    ext : str
        File extension ('.trr' or '.xtc').
    minidx : int
        First frame index (0-based) to include.
    maxidx : int
        Last frame index (0-based) to include.
    cut_filename : str
        Output trajectory filename.
    gmx_exe : str
        Path to the GROMACS executable.
    """
    ndx_file = cut_filename + '.ndx'
    with open(ndx_file, 'w', encoding='utf-8') as ndxf:
        ndxf.write('[ frames ]\n')
        # trjconv uses 1-based frame indices; offset by 1 so GROMACS
        # does not drop frame 0.
        for i in range(minidx, maxidx + 1):
            ndxf.write(f'{i + 1}\n')

    subprocess.run(
        [gmx_exe, 'trjconv',
         '-f', filename,
         '-o', cut_filename,
         '-fr', ndx_file],
        input=b'0\n',
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False)


def _cut_xyz(filename, minidx, maxidx, cut_filename):
    """Cut an XYZ trajectory segment.

    Parameters
    ----------
    filename : str
        Path to the source XYZ trajectory file.
    minidx : int
        First frame index (0-based) to include.
    maxidx : int
        Last frame index (0-based) to include.
    cut_filename : str
        Output trajectory filename.
    """
    with open(filename, 'r', encoding='utf-8') as f_in:
        first_line = f_in.readline()
        if not first_line:
            return
        try:
            natoms = int(first_line.strip())
        except ValueError:
            return
        lines_per_frame = natoms + 2

        f_in.seek(0)
        for _ in range(minidx * lines_per_frame):
            if not f_in.readline():
                break

        with open(cut_filename, 'w', encoding='utf-8') as f_out:
            for _ in range((maxidx - minidx + 1) * lines_per_frame):
                line = f_in.readline()
                if not line:
                    break
                f_out.write(line)


def do_cut_trajectory(filename, ext, minidx, maxidx,
                      cut_filename, engine_info):
    """Physically cut a trajectory segment into a temporary file.

    Parameters
    ----------
    filename : string
        The name of the original trajectory file.
    ext : string
        The file extension (e.g., '.trr', '.xtc', '.xyz').
    minidx : int
        The first frame index to include.
    maxidx : int
        The last frame index to include.
    cut_filename : string
        The name of the new temporary trajectory file.
    engine_info : dict
        A dictionary with engine-specific tools (e.g., GMX executable).
    """
    if ext in ['.trr', '.xtc']:
        gmx_exe = engine_info.get('gmx', 'gmx')
        _cut_gromacs(filename, ext, minidx, maxidx, cut_filename, gmx_exe)
    elif ext == '.xyz':
        _cut_xyz(filename, minidx, maxidx, cut_filename)
    else:
        shutil.copy(filename, cut_filename)

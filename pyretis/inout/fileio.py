# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Module defining the base classes for the PyRETIS output.

Important classes defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

FileIO (:py:class:`.FileIO`)
    A generic class for handling input & output with files.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

read_some_lines (:py:func:`.read_some_lines`)
    Method to read lines from PyRETIS data files.

"""
import os
import logging
from pyretis.inout.common import OutputBase, create_backup


logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


__all__ = ['FileIO', 'read_some_lines']


class FileIO(OutputBase):
    """A generic class for handling IO with files.

    This class defines how PyRETIS stores and reads data.
    Formatting is handled by an object like :py:class:`.OutputFormatter`

    Attributes
    ----------
    filename : string
        Name (e.g. path) to the file to read or write.
    file_mode : string
        Specifies the mode in which the file is opened.
    backup : boolean
        Determines the behavior if we want to write to a file
        that is already existing.
    fileh : object like :py:class:`io.IOBase`
        The file handle we are interacting with.
    last_flush : object like :py:class:`datetime.datetime`
        The previous time for flushing to the file.

    """

    target = 'file'

    def __init__(self, filename, file_mode, formatter, backup=True):
        """Set up the file object.

        Parameters
        ----------
        filename : string
            The path to the file to open or read.
        file_mode : string
            Specifies the mode for opening the file.
        formatter : object like py:class:`.OutputFormatter`
            The object responsible for formatting output.
        backup : boolean, optional
            Defines how we handle cases where we write to a
            file which is already existing.

        """
        super().__init__(formatter)
        self.filename = filename
        self.file_mode = file_mode
        if backup not in (True, False):
            logger.info('Setting backup to default: True')
            self.backup = True
        else:
            self.backup = backup
        self.fileh = None
        if self.file_mode.startswith('a') and self.formatter is not None:
            self.formatter.print_header = False
        self.last_flush = None

    def open_file_read(self):
        """Open a file for reading."""
        if not self.file_mode.startswith('r'):
            raise ValueError(
                f'Inconsistent file mode "{self.file_mode}" for reading'
            )
        try:
            encoding = 'utf-8' if 'b' in self.file_mode else None
            self.fileh = open(self.filename, self.file_mode, encoding=encoding)
        except (OSError, IOError) as error:
            logger.warning(
                'Could not open file "%s" for reading', self.filename
            )
            logger.warning(
                'I/O error (%d): %s', error.errno, error.strerror
            )
        return self.fileh

    def open_file_write(self):
        """Open a file for writing.

        In this method, we also handle the possible backup settings.
        """
        if not self.file_mode[0] in ('a', 'w'):
            raise ValueError(
                f'Inconsistent file mode "{self.file_mode}" for writing'
            )
        msg = []
        try:
            if os.path.isfile(self.filename):
                msg = ''
                if self.file_mode.startswith('a'):
                    logger.info(
                        'Appending to existing file "%s"', self.filename
                    )
                else:
                    if self.backup:
                        msg = create_backup(self.filename)
                        logger.debug(msg)
                    else:
                        logger.debug(
                            'Overwriting existing file "%s"', self.filename
                        )
            encoding = 'utf-8' if 'b' in self.file_mode else None
            self.fileh = open(self.filename, self.file_mode, encoding=encoding)
        except (OSError, IOError) as error:  # pragma: no cover
            logger.warning(
                'Could not open file "%s" for writing', self.filename
            )
            logger.warning(
                'I/O error (%d): %s', error.errno, error.strerror
            )
        return self.fileh

    def open(self):
        """Open a file for reading or writing."""
        if self.fileh is not None:
            logger.debug(
                '%s asked to open file, but it has already opened a file.',
                self.__class__.__name__
            )
            return self.fileh
        if self.file_mode[0] in ('r',):
            return self.open_file_read()
        if self.file_mode[0] in ('a', 'w'):
            return self.open_file_write()
        raise ValueError(f'Unknown file mode "{self.file_mode}"')

    def load(self):
        """Read blocks or lines from the file."""
        return self.formatter.load(self.filename)

    def write(self, towrite, end='\n'):
        """Write a string to the file.

        Parameters
        ----------
        towrite : string
            The string to output to the file.
        end : string, optional
            Appended to `towrite` when writing, can be used to print a
            new line after the input `towrite`.

        Returns
        -------
        status : boolean
            True if we managed to write, False otherwise.

        """
        status = False
        if towrite is None:
            return status
        if self.fileh is not None and not self.fileh.closed:
            try:
                if end is not None:
                    self.fileh.write(f'{towrite}{end}')
                    status = True
                else:
                    self.fileh.write(towrite)
                    status = True
            except (OSError, IOError) as error:  # pragma: no cover
                msg = f'Write I/O error ({error.errno}): {error.strerror}'
                logger.critical(msg)
            self.flush()
            return status
        if self.fileh is not None and self.fileh.closed:
            logger.warning('Ignored writing to closed file %s', self.filename)
        if self.fileh is None:
            logger.critical(
                'Attempting to write to empty file handle for file %s',
                self.filename
            )
        return status

    def close(self):
        """Close the file."""
        if self.fileh is not None and not self.fileh.closed:
            try:
                self.flush()
            finally:
                self.fileh.close()

    def flush(self):
        """Flush file buffers to file."""
        if self.fileh is not None and not self.fileh.closed:
            self.fileh.flush()
            os.fsync(self.fileh.fileno())

    def output(self, step, data):
        """Open file before first write."""
        if self.first_write:
            self.open()
        return super().output(step, data)

    def __del__(self):
        """Close the file in case the object is deleted."""
        self.close()

    def __enter__(self):
        """Context manager for opening the file."""
        self.open()
        return self

    def __exit__(self, *args):
        """Context manager for closing the file."""
        self.close()

    def __iter__(self):
        """Make it possible to iterate over lines in the file."""
        return self

    def __next__(self):
        """Let the file object handle the iteration."""
        if self.fileh is None:
            raise StopIteration
        if self.fileh.closed:
            raise StopIteration
        return next(self.fileh)

    def __str__(self):
        """Return basic info."""
        msg = [f'FileIO (file: "{self.filename}")']
        if self.fileh is not None and not self.fileh.closed:
            msg += ['\t* File is open']
            msg += [f'\t* Mode: {self.fileh.mode}']
        msg += [f'\t* Formatter: {self.formatter}']
        return '\n'.join(msg)


def _read_line_data(ncol, stripline, line_parser):
    """Read data for :py:func:`.read_some_lines.`.

    Parameters
    ----------
    ncol : integer
        The expected number of columns to read. If this is less than 1
        it is not yet set. Note that we skip data which appear
        inconsistent. A warning will be issued about this.
    stripline : string
        The line to read. Note that we assume that leading and
        trailing spaces have been removed.
    line_parser : callable
        A method we use to parse a single line.

    """
    if line_parser is None:
        # Just return data without any parsing:
        return stripline, True, ncol
    try:
        linedata = line_parser(stripline)
    except (ValueError, IndexError):
        return None, False, -1
    newcol = len(linedata)
    if ncol == -1:  # first item
        ncol = newcol
    if newcol == ncol:
        return linedata, True, ncol
    # We assume that this is line is malformed --- skip it!
    return None, False, -1


def read_some_lines(filename, line_parser, block_label='#'):
    """Open a file and try to read as many lines as possible.

    This method will read a file using the given `line_parser`.
    If the given `line_parser` fails at a line in the file,
    `read_some_lines` will stop here. Further, this method
    will read data in blocks and yield a block when a new
    block is found. A special string (`block_label`) is assumed to
    identify the start of blocks.

    Parameters
    ----------
    filename : string
        This is the name/path of the file to open and read.
    line_parser : function, optional
        This is a function which knows how to translate a given line
        to a desired internal format. If not given, a simple float
        will be used.
    block_label : string, optional
        This string is used to identify blocks.

    Yields
    ------
    data : list
        The data read from the file, arranged in dicts.

    """
    ncol = -1  # The number of columns
    new_block = {'comment': [], 'data': []}
    yield_block = False
    read_comment = False
    with open(filename, 'r', encoding='utf-8') as fileh:
        for i, line in enumerate(fileh):
            stripline = line.strip()
            if stripline.startswith(block_label):
                # this is a comment, then a new block will follow,
                # unless this is a multi-line comment.
                if read_comment:  # part of multi-line comment...
                    new_block['comment'].append(stripline)
                else:
                    if yield_block:
                        # Yield the current block
                        yield_block = False
                        yield new_block
                    new_block = {'comment': [stripline], 'data': []}
                    yield_block = True  # Data has been added
                    ncol = -1
                    read_comment = True
            else:
                read_comment = False
                data, _yieldb, _ncol = _read_line_data(ncol, stripline,
                                                       line_parser)
                if data:
                    new_block['data'].append(data)
                    ncol = _ncol
                    yield_block = _yieldb
                else:
                    logger.warning('Skipped malformed data in "%s", line: %i',
                                   filename, i)
    # if the block has not been yielded, yield it
    if yield_block:
        yield_block = False
        yield new_block

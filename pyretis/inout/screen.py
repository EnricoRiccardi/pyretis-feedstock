# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Module defining the base classes for the PyRETIS output.

Important classes defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ScreenOutput (:py:class:`.FileIO`)
    A generic class for handling output to the screen.

Important constants defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROGRESS : int
    Custom log level (25) between INFO (20) and WARNING (30).
    Used for user-facing progress messages (green on console).

BANNER : int
    Custom log level (26) between PROGRESS (25) and WARNING (30).
    Used for decorative/banner text (cyan on console): logo,
    version info, references.

"""
import logging
from pyretis.inout.common import OutputBase


logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())

# Custom log levels for user-facing console output.
# Console handler is set to PROGRESS so only PROGRESS, BANNER,
# WARNING, ERROR, CRITICAL appear on screen.
# logger.info() and logger.debug() go to the log file only.
PROGRESS = 25
BANNER = 26
REFERENCE = 27
logging.addLevelName(PROGRESS, 'PROGRESS')
logging.addLevelName(BANNER, 'BANNER')
logging.addLevelName(REFERENCE, 'REFERENCE')


def log_progress(self, message, *args, **kwargs):
    """Log a message at the PROGRESS level."""
    self.log(PROGRESS, message, *args, **kwargs)


def log_banner(self, message, *args, **kwargs):
    """Log a message at the BANNER level."""
    self.log(BANNER, message, *args, **kwargs)


def log_reference(self, message, *args, **kwargs):
    """Log a message at the REFERENCE level (white on console)."""
    self.log(REFERENCE, message, *args, **kwargs)


logging.Logger.progress = log_progress
logging.Logger.banner = log_banner
logging.Logger.reference = log_reference


__all__ = ['ScreenOutput', 'PROGRESS', 'BANNER', 'REFERENCE']  # noqa: F401


class ScreenOutput(OutputBase):
    """A class for handling output to the screen."""

    target = 'screen'

    def write(self, towrite, end=None):
        """Write a string to the file.

        Parameters
        ----------
        towrite : string
            The string to output to the file.
        end : string, optional
            Override how the print statements ends.

        Returns
        -------
        status : boolean
            True if we managed to write, False otherwise.

        """
        if towrite is None:
            return False
        if end is not None:
            print(towrite, end=end)
            return True
        print(towrite)
        return True

# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""PyRETIS style for use with pygments.

This style is based on `Xcode` from the Pygments.
For the original style:

    :pyretis_copyright: Copyright 2006-2029 by the Pygments team, see AUTHORS.
    :license: BSD, see LICENSE for details.

See: http://pygments.org/
"""

from pygments.style import Style
from pygments.token import Keyword, Name, Comment, String, Error, \
     Number, Operator, Literal, Generic, Punctuation


def setup(app):
    """Register the script."""
    pass


class PyretisStyle(Style):
    """Style similar to the Xcode default colouring theme."""

    default_style = '#262626'

    styles = {
        Comment:                '#177500',
        Comment.Preproc:        '#633820',
        String:                 '#C41A16',
        String.Char:            '#2300CE',
        Operator:               '#262626',
        Keyword:                '#A90D91',
        Name:                   '#262626',
        Name.Attribute:         '#836C28',
        Name.Class:             '#3F6E75',
        Name.Function:          '#262626',
        Name.Builtin:           '#A90D91',
        Name.Builtin.Pseudo:    '#5B269A',
        Name.Variable:          '#262626',
        Name.Tag:               '#262626',
        Name.Decorator:         '#262626',
        Name.Label:             '#262626',
        Punctuation:            '#262626',
        Literal:                '#1C01CE',
        Number:                 '#1C01CE',
        Error:                  '#262626',
        Generic:                "#262626",
        Generic.Deleted:        "#A00000",
        Generic.Emph:           "italic",
        Generic.Error:          "#FF0000",
        Generic.Heading:        "bold #000080",
        Generic.Inserted:       "#00A000",
        Generic.Output:         "#262626",
        Generic.Prompt:         "bold #000080",
        Generic.Strong:         "bold",
        Generic.Subheading:     "bold #800080",
        Generic.Traceback:      "#0044DD",
    }

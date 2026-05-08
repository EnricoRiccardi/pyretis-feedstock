# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Extension to color the name of the app + commands in code snippets."""
from pygments.lexer import inherit
from pygments.lexers.shell import BashLexer
from pygments.token import Name
from sphinx import highlighting


KEYWORDS = (
    'pip',
    'pip3',
    'git',
    'f2py',
    'pyretisrun',
    'pyretisanalyse',
    'virtualenv',
    'conda',
    'sudo',
    'mkdir',
    'pylint',
    'pycodestyle',
    'pydocstyle',
)


class PyretisLexer(BashLexer):
    """Add a few keyword tokes. These will be colored."""

    name = 'Pyretis in bash'
    aliases = ['pyretis']
    mimetypes = ['application/x-sh', 'application/x-shellscript']
    merged = '|'.join(KEYWORDS)
    tokens = {
        'basic': [
            (r'\b({})(?=[\s)`])'.format(merged),
             Name.Builtin), inherit,
        ],
    }


def setup(app):
    pass


highlighting.lexers['pyretis'] = PyretisLexer()

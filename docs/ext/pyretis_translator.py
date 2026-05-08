# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Remove the container class from <div> items.

This is just to remove some unwanted containers in the HTML
output.

"""
from sphinx.writers import html


class PyretisTranslator(html.HTMLTranslator):
    """Create a translator that removes the container class."""

    def visit_container(self, node):
        self.body.append(self.starttag(node, 'div', CLASS=''))


def setup(app):
    app.set_translator('html', PyretisTranslator)

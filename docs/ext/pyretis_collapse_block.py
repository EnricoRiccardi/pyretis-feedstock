# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Defining a collapse block for the PyRETIS web-page."""
from docutils import nodes
from pathlib import Path
from docutils.parsers.rst import Directive
from docutils.parsers.rst.directives import unchanged
from jinja2 import Template


HTML_PRE = """
<p>
<a class="btn btn-primary btn-sm" data-toggle="collapse" href="#{{ id }}">
{{ heading }} &raquo;</a>
  <div id="{{ id }}" class="collapse">"""
TEMPLATE_PRE = Template(HTML_PRE)
HTML_POST = """  </div></p>"""
TEMPLATE_POST = Template(HTML_POST)


class PyretisCollapseNode(nodes.General, nodes.Element):
    """Create a new node for the collapse block contents."""
    pass


class PyretisCollapseBlock(Directive):
    """Define a new directive."""

    required_arguments = 0
    optional_arguments = 0
    option_spec = {'heading': unchanged}
    has_content = True

    def run(self):
        """Set up a new node for a collapse block."""
        env = self.state.document.settings.env
        rst_source = Path(self.state_machine.node.document['source'])
        rst_filename = rst_source.name
        target_id = "%s.pyretiscb-%d" % (
            rst_filename,
            env.new_serialno('pyretis-ccb')
        )
        target_id = target_id.replace(".", "-")
        # target_node = nodes.target('', '', ids=[target_id])

        node = PyretisCollapseNode()
        node['heading'] = self.options.get('heading', 'title')
        node['target_id'] = target_id
        self.state.nested_parse(self.content, self.content_offset, node)
        return [node]


def html_visit_block(self, node):
    """Render the template when the block is visited."""
    self.body.append(
        TEMPLATE_PRE.render(id=node['target_id'],
                            heading=node['heading'])
    )


def html_depart_block(self, _):
    """End the HTML code when the block is done."""
    self.body.append(TEMPLATE_POST.render())


def setup(app):
    """Register the collable block."""
    app.add_node(
        PyretisCollapseNode,
        html=(
            html_visit_block,
            html_depart_block
        )
    )
    app.add_directive('pyretis-collapse-block', PyretisCollapseBlock)

# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Create a directive for showing input examples."""
from docutils import nodes
from docutils.parsers.rst import Directive
from docutils.parsers.rst.directives import unchanged
from jinja2 import Template


HTML_PRE = """
<div class="container">
    <div class="panel panel-keyword">
      <div class="panel-heading">{{ heading }}</div>
      <div class="panel-body">"""
TEMPLATE_PRE = Template(HTML_PRE)
HTML_POST = """</div>
  </div>
</div>"""
TEMPLATE_POST = Template(HTML_POST)


class PyretisExampleNode(nodes.General, nodes.Element):
    pass


class PyretisInputExample(Directive):
    """Create a new directive for PyRETIS examples."""

    required_arguments = 0
    optional_arguments = 1
    option_spec = {'class-name': unchanged}
    has_content = True

    def run(self):
        node = PyretisExampleNode()
        class_name = self.options.get('class-name', None)
        if not self.arguments:
            section = 'Example section'
        else:
            section = 'Example {} section'.format(self.arguments[0])
        if class_name:
            node['heading'] = '{} for {}:'.format(section, class_name)
        else:
            node['heading'] = '{}:'.format(section)
        self.state.nested_parse(self.content, self.content_offset, node)
        return [node]


def html_visit_input_example_block(self, node):
    self.body.append(
        TEMPLATE_PRE.render(
            heading=node['heading']
        )
    )


def html_depart_input_example_block(self, node):
    self.body.append(TEMPLATE_POST.render())


def setup(app):
    app.add_node(
        PyretisExampleNode,
        html=(
            html_visit_input_example_block,
            html_depart_input_example_block
        )
    )
    app.add_directive('pyretis-input-example', PyretisInputExample)

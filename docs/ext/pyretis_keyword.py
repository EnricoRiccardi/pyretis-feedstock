# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Define a new directive for PyRETIS keywords."""
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


class PyretisKeywordNode(nodes.General, nodes.Element):
    pass


class PyretisKeyword(Directive):

    required_arguments = 2
    optional_arguments = 5
    option_spec = {'specific': unchanged}
    has_content = True

    def run(self):
        node = PyretisKeywordNode()
        if len(self.arguments) == 2:
            arguments = [i for i in self.arguments]
        else:
            if self.arguments[0].lower() == 'parameter':
                arguments = [
                    '{} {}'.format(
                        self.arguments[0],
                        self.arguments[1]
                    )
                ]
                arguments += [' '.join(self.arguments[2:])]
            else:
                arguments = ([self.arguments[0]] +
                             [' '.join(self.arguments[1:])])
        if self.options.get('specific', 'no') == 'yes':
            node['heading'] = '{} = {}'.format(*arguments)
        else:
            node['heading'] = '{} = <i>{}</i>'.format(*arguments)
        self.state.nested_parse(self.content, self.content_offset, node)
        return [node]


def html_visit_keyword_block(self, node):
    self.body.append(
        TEMPLATE_PRE.render(
            heading=node['heading']
        )
    )


def html_depart_keyword_block(self, node):
    self.body.append(TEMPLATE_POST.render())


def setup(app):
    app.add_node(
        PyretisKeywordNode,
        html=(
            html_visit_keyword_block,
            html_depart_keyword_block
        )
    )
    app.add_directive('pyretis-keyword', PyretisKeyword)

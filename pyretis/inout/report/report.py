# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""General functions for generating reports.

This module contains some general functions for report generation. These
functions are used by the specific report generators to format the
reports.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

get_template (:py:func:`.get_template`)
    Returns the template for a specific output format and report type.

render_report (:py:func:`.render_report`)
    Render a report using a template and jinja2.

generate_report (:py:func:`.generate_report`)
    Generate a specific report from analysis output.
"""
import datetime
import logging
import os
# for converting rst to html and/or latex:
import docutils.core
from docutils.writers.html4css1 import Writer as HTMLWriter
from docutils.writers.html4css1 import HTMLTranslator
# for using templates
import jinja2
from pyretis import __version__ as VERSION
from pyretis.info import PROGRAM_NAME
from pyretis.inout.report.markup import latexify_number
from pyretis.inout.report.report_md import generate_report_mdflux
from pyretis.inout.report.report_path import (generate_report_retis,
                                              generate_report_retis0,
                                              generate_report_tis,
                                              generate_report_tis_path,
                                              generate_report_repptis,
                                              generate_report_pptis_path)
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


__all__ = ['get_template', 'render_report', 'generate_report']


_DATE_FMT = '%d.%m.%Y %H:%M:%S'
# File names for pre-defined templates.
# - html is done via rst (i.e. there is no html template)
# - htm is assumed to be equal to html
# - tex is assumed to be equal to latex
_TEMPLATES = {'rst': 'report_{}.rst',
              'html': 'report_{}.rst',
              'htm': 'report_{}.rst',
              'latex': 'report_{}.tex',
              'tex': 'report_{}.tex',
              'txt': 'report_{}.txt'}

_TEMPLATE_NAMES = {'md-flux': 'mdflux',
                   'retis': 'retis',
                   'retis0': 'retis0',
                   'tis': 'tis',
                   'repptis': 'repptis',
                   'pptis': 'pptis',
                   'make-tis-files': 'tis-multiple'}

_REPORT_MAP = {'md-flux': generate_report_mdflux,
               'retis': generate_report_retis,
               'retis0': generate_report_retis0,
               'make-tis-files': generate_report_tis,
               'repptis': generate_report_repptis,
               'pptis': generate_report_pptis_path,
               'tis': generate_report_tis_path}

# Table for file extensions:
_EXT = {'rst': 'rst',
        'html': 'html',
        'htm': 'htm',  # in case some people prefer it
        'latex': 'tex',
        'tex': 'tex',
        'txt': 'txt'}


def _rst_to_html(rst):
    """Convert a reStrcuturedText string to simple HTML.

    Parameters
    ----------
    rst : string
        The string to convert.

    Returns
    -------
    out : string
        A HTML document corresponding to the input reStructuredText.

    """
    htmlwriter = HTMLWriter()
    htmlwriter.translator_class = HTMLTranslator
    override = {'output_encoding': 'unicode'}
    # custom css can be added by: 'stylesheet_path': '/path/to/style.css'
    html = docutils.core.publish_string(rst, writer=htmlwriter,
                                        settings_overrides=override)
    return html


def get_template(output, report_type, template=None):
    """Return the template to use for a specified output format.

    The output is one of the defined output types, for instance 'rst'
    for restructured text or 'latex' for latex. Different report types
    will have different templates and the report types must also be
    specified here.

    Parameters
    ----------
    output : string
        This string selects the output format for the template, i.e.,
        rst, html, latex, tex.
    template : string, optional
        The full path to the template to use. If not given/found, the
        defaults in :py:const:`._TEMPLATES` will be used.
    report_type : string
        This is the type of report we are doing, e.g. TIS or MD.

    Returns
    -------
    out[0] : string
        Filename of template to use.
    out[1] : string
        The path to the template to use.

    """
    if template is None or not os.path.isfile(template):
        # Use default template, this is located in the templates dir:
        path = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(path, 'templates')
        ltype = report_type.lower()
        template = _TEMPLATES[output].format(_TEMPLATE_NAMES[ltype])
        path_to_template = os.path.join(path, template)
        if not os.path.isfile(path_to_template):
            msg = f'Could not locate template "{path_to_template}"!'
            raise ValueError(msg)
    else:
        # user specified full path to template:
        path = os.path.dirname(template)
        template = os.path.basename(template)
    return template, path


def render_report(report, output, template, path):
    """Render a report using a template and jinja2.

    The report is given as a dictionary which is used to fill in a
    template with jinja2. The template is given as a string
    (a file name) with a path to the template. The output can also be
    specified here and this is only used to convert to HTML if that is
    the desired output.

    Parameters
    ----------
    report : dict
        This dict contains the data to be reported. It is assumed that
        this dict matches the specified template.
    output : string
        This is the desired output format. Here it's used only for
        generating HTML as this is done via rst.
    template : string
        This is the template to use (the file name).
    path :  string
        This is the template file to use (its full path).

    Returns
    -------
    out[0] : string
        The generated report in the desired format.
    out[1] : string
        The file extension (i.e. file type) for the generated report.

    Note
    ----
    The parameters `template` and `path` are typically obtained by a
    call to :py:func:`.get_template`.

    """
    env = jinja2.Environment(
        block_start_string='@{%',
        block_end_string='%}@',
        variable_start_string='@{{',
        variable_end_string='}}@',
        autoescape=jinja2.select_autoescape(['html', 'htm', 'xml']),
        loader=jinja2.FileSystemLoader(path)
    )
    # pylint: disable=maybe-no-member
    render = env.get_template(template).render(report)
    # pylint: enable=maybe-no-member
    if output == 'html':
        return _rst_to_html(render), _EXT[output]
    return render, _EXT[output]


def generate_report(report_type, analysis_results, output, template=None):
    """Generate a report of a given type with the given analysis results.

    Parameters
    ----------
    report_type : string
        Selects the kind of report we want.
    analysis_results : dict
        The results from running the analysis.
    output : string
        Output format for the report.
    template : string, optional
        The full path to the template to use. If not given/found, the
        defaults in :py:const`._TEMPLATES` will be used. This is handled
        by :py:func:`.get_template`.

    Returns
    -------
    out[0] : string
        This is the generated report in the desired format
    out[1] : string
        This is an extension which can be used when writing the report
        to a file.

    """
    report = {'version': VERSION,
              'program': PROGRAM_NAME,
              'date': datetime.datetime.now().strftime(_DATE_FMT),
              'figures': {}, 'tables': {}, 'numbers': {}}
    if report_type not in _REPORT_MAP:
        raise ValueError(f'Unkown report type {report_type}')
    generator = _REPORT_MAP[report_type]
    # Check if the output is a valid format
    if output not in _TEMPLATES:
        msg = 'Format {} not defined for {} report. Defaulting to rst'
        msg = msg.format(output, report_type)
        logger.warning(msg)
        output = 'rst'
    template, path = get_template(output, report_type, template=template)
    generated = generator(analysis_results, output=output)
    report.update(generated)
    # Remove white-space from numbers:
    for key, value in report['numbers'].items():
        if isinstance(value, str):
            report['numbers'][key] = value.strip()

    if output in ('latex', 'tex', 'html', 'htm'):
        # Latexify numbers:
        for key, value in report['numbers'].items():
            if isinstance(value, str):
                report['numbers'][key] = latexify_number(value)
    return render_report(report, output, template, path)

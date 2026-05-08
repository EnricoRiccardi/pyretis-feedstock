# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Just plot some results using bokeh."""
import numpy as np
from bokeh.plotting import figure, output_file, show
from bokeh.models import Range1d, LinearAxis
from bokeh.palettes import all_palettes
from bokeh.layouts import gridplot


def make_plot():
    """Create plot using bokeh."""
    colors = all_palettes['Set1'][4]

    raw_data = np.loadtxt('thermo.txt')
    data = {}
    for i, key in enumerate(('step', 'temp', 'pot', 'kin', 'tot', 'press')):
        data[key] = raw_data[:, i]

    output_file('thermo_plot.html')
    tot0 = data['tot'][0]
    plt1 = figure(x_range=Range1d(data['step'][0], data['step'][-1]*1.05),
                  width=400, height=400, toolbar_location='above')
    plt1.line(data['step'], data['pot']-tot0, line_width=3,
              line_color=colors[0], alpha=0.8, legend='Potential')
    plt1.line(data['step'], data['kin']-tot0, line_width=3,
              line_color=colors[1], alpha=0.8, legend='Kinetic')
    plt1.line(data['step'], data['tot']-tot0, line_width=3,
              line_color=colors[2], alpha=0.8, legend='Total')

    plt1.xaxis.axis_label = 'Step'
    plt1.yaxis.axis_label = 'Energy per particle'
    plt1.legend.location = 'center_left'
    plt1.legend.label_text_font_size = '12pt'
    plt1.background_fill_color = '#E8DDCB'
    plt1.background_fill_alpha = 0.5
    plt1.xaxis.axis_label_text_font_size = '12pt'
    plt1.yaxis.axis_label_text_font_size = '12pt'
    plt1.xaxis.major_label_text_font_size = '12pt'
    plt1.yaxis.major_label_text_font_size = '12pt'

    plt2 = figure(x_range=Range1d(data['step'][0], data['step'][-1]),
                  width=400, height=400, toolbar_location='above')
    plt2.line(data['step'], data['press'], line_width=3,
              line_color=colors[1], alpha=0.8, legend='Pressure')
    plt2.xaxis.axis_label = 'Step'
    plt2.yaxis.axis_label = 'Pressure'
    plt2.legend.location = 'center_left'
    plt2.background_fill_color = '#E8DDCB'
    plt2.background_fill_alpha = 0.5

    plt2.extra_y_ranges = {'temp': Range1d(start=data['temp'].min(),
                                           end=data['temp'].max()*1.05)}
    plt2.line(data['step'], data['temp'], line_width=3,
              line_color=colors[2], alpha=0.8, legend='Temperature',
              y_range_name='temp')
    plt2.add_layout(LinearAxis(y_range_name='temp', axis_label='Temperature'),
                    'right')
    plt2.xaxis.axis_label_text_font_size = '12pt'
    plt2.yaxis.axis_label_text_font_size = '12pt'
    plt2.legend.label_text_font_size = '12pt'

    plt2.xaxis.major_label_text_font_size = '12pt'
    plt2.yaxis.major_label_text_font_size = '12pt'

    subplot = gridplot(
        (plt1, plt2),
        ncols=2,
        plot_width=400,
        plot_height=400,
    )
    show(subplot)


if __name__ == '__main__':
    make_plot()

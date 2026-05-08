# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Plot the object relation graphs."""
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle


mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'DejaVu Serif'
mpl.rcParams['savefig.directory'] = ''
mpl.rcParams['savefig.format'] = 'svg'


OBJECTS = [
    {
        'title': 'System',
        'methods': ['add_particle()'],
        'attributes-ref': ['particles', 'box', 'forcefield'],
    },
    {
        'title': 'Box',
        'methods': ['pbc_wrap()'],
    },
    {
        'title': 'Particles',
        'methods': ['get_pos()', 'set_pos()'],
        'attributes': ['pos', 'vel'],
    },
    {
        'title': 'ForceField',
        'attributes-ref': ['potential'],
    },
    {
        'title': 'PathEnsemble',
        'methods': ['store_path()'],
        'attributes': ['paths'],
        'attributes-ref': ['last_path'],
    },
    {
        'title': 'Path',
        'methods': ['append()'],
        'attributes-ref': ['phasepoints'],
    },
    {
        'title': 'PotentialFunction',
        'methods': ['potential()', 'force()', 'potential_and_force()'],
    },
    {
        'title': 'System',
        'attributes-ref': ['particles', 'forcefield'],
    },
    {
        'title': 'Simulation',
        'methods': ['step()', 'run()'],
        'attributes-ref': ['system', 'engine'],
    },
    {
        'title': 'ExternalMDEngine',
        'methods': ['step()', 'propagate()']
    },
    {
        'title': 'MDEngine',
        'methods': ['integration_step()', 'propagate()'],
    },
    {
        'title': 'ParticlesExt',
        'methods': ['get_pos()', 'set_pos()'],
        'attributes': ['config', 'vel_rev'],
    },
    {
        'title': 'System',
        'attributes-ref': ['particles', 'box', 'forcefield'],
    },
    {
        'title': 'EngineBase',
        'methods': ['integration_step()', 'propagate()'],
    },
]


def draw_object(title, methods, attributes, attributes_ref):
    """Draw a simple box representing an object."""
    origin = (0, 0)

    longest = max(
        len(title)*1.15,
        max([len(i) for i in methods + attributes + attributes_ref])
    )
    print('Longest text:', longest)

    width = longest * 6

    linew = 2
    txtheight = 18

    meth_height = 10
    meth_sep = 2
    attr_height = 10
    attr_sep = 2
    attr_ref_height = 15
    attr_ref_sep = 3

    initial_sep = 0

    if methods:
        initial_sep = 1
    else:
        if attributes:
            initial_sep = 1
        else:
            initial_sep = 6

    last_sep = 0
    if attributes_ref:
        last_sep = 3
    else:
        if attributes:
            last_sep = 3
        else:
            last_sep = 3
    if attributes_ref:
        extra = 2
    else:
        extra = 0

    height = (linew + txtheight + initial_sep + last_sep + extra +
              len(methods) * meth_height +
              len(attributes) * attr_height +
              len(attributes_ref) * attr_ref_height +
              max(len(attributes) - 1, 0) * attr_sep +
              max(len(methods) - 1, 0) * meth_sep +
              max(len(attributes_ref) - 1, 0) * attr_ref_sep)

    if attributes_ref:
        last_sep = 20
    else:
        if attributes:
            last_sep = 15
        else:
            last_sep = 15

    rtitle = Rectangle(origin, width, height, color='#43ac6a')
    origin2 = (origin[0] + linew, origin[1] + linew)
    width2 = width - 2*linew
    height2 = height - txtheight - linew
    rcontents = Rectangle(origin2, width2, height2, color='white')

    ytxt = origin[1] + height - txtheight*0.5
    xtxt = origin[0] + 2

    txt_title = mpl.text.Text(xtxt, ytxt, title, fontsize=22, color='white',
                              verticalalignment='center',
                              horizontalalignment='left')

    txt_box = [txt_title]

    xtxt += linew
    ytxt = height2 - initial_sep

    for met in methods:
        txt = mpl.text.Text(xtxt, ytxt, met, fontsize=18,
                            style='normal',
                            color='#262626',
                            verticalalignment='top',
                            horizontalalignment='left')
        txt_box.append(txt)
        ytxt -= meth_height
        ytxt -= meth_sep

    for met in attributes:
        txt = mpl.text.Text(xtxt, ytxt, met, fontsize=18,
                            style='italic',
                            color='#262626',
                            verticalalignment='top',
                            horizontalalignment='left')
        txt_box.append(txt)
        ytxt -= attr_height
        ytxt -= attr_sep

    xtxt += 2.0
    ytxt -= extra
    for attr in attributes_ref:
        txt = mpl.text.Text(
            xtxt, ytxt, attr, fontsize=18,
            style='italic', color='#262626',
            verticalalignment='top', horizontalalignment='left',
            bbox={'boxstyle':
                  'roundtooth,pad=0.25',
                  'fc': '#e6e6e6',
                  'ec': '#333333'}
        )
        # Different options for the bbox:
        # {'boxstyle': 'round,pad=0.25',
        # {'boxstyle': 'round4,pad=0.25',
        # {'boxstyle': 'sawtooth,pad=0.25',
        ytxt -= attr_ref_height
        ytxt -= attr_ref_sep
        txt_box.append(txt)
    return (rtitle, rcontents), txt_box


def make_figure(obj):
    """Just plot a figure."""
    patches, txt = draw_object(
        obj['title'],
        obj.get('methods', []),
        obj.get('attributes', []),
        obj.get('attributes-ref', [])
    )

    fig = plt.figure()
    ax1 = fig.add_subplot(111, aspect='equal')
    for patch in patches:
        ax1.add_patch(patch)
    for box in txt:
        ax1.add_artist(box)
    ax1.set_xlim(-10, 150)
    ax1.set_ylim(-10, 110)
    plt.axis('off')
    plt.show()


def make_txt_figure(txts):
    """Just plot a figure."""

    fig = plt.figure()
    ax1 = fig.add_subplot(111, aspect='equal')
    xtxt = 10
    ytxt = 10
    for txt in txts:
        txtbox = mpl.text.Text(
            xtxt, ytxt, txt, fontsize=14,
            color='#ffffff',
            style='italic',
            verticalalignment='top',
            horizontalalignment='left',
            bbox={'boxstyle': 'round,pad=0.25',
                  'fc': '#34495e',
                  'ec': '#262626'}
        )
        ytxt += 20
        ax1.add_artist(txtbox)
    ax1.set_xlim(-10, 110)
    ax1.set_ylim(-10, 110)
    plt.axis('off')
    plt.show()


def main():
    """Do all plotting."""
    for obj in OBJECTS:
        make_figure(obj)

    var = ["tuple: ('file.xyz', 123)",
           "boolean: True/False",
           "numpy.array",
           r"return copy of numpy.array",
           r"return ('file.xyz', 123)"]
    make_txt_figure(var)


if __name__ == '__main__':
    main()

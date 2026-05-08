# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Example of running a MD NVE simulation."""
# pylint: disable=invalid-name
import numpy as np
# imports for the plotting:
from matplotlib import pyplot as plt
from matplotlib import animation
import matplotlib as mpl
from matplotlib import gridspec
from pyretis.core import System, create_box, Particles
from pyretis.core.units import CONVERT, create_conversion_factors
from pyretis.inout.plotting import COLORS, COLOR_SCHEME
from pyretis.setup import (create_system, create_engine,
                           create_force_field, create_simulation)


PCOLOR = {'A': 'blue', 'B': 'magenta'}  # Colors for drawing
# Define potential parameters:
WCA_PARAMETERS = {
    0: {'sigma': 1.0, 'epsilon': 1.0, 'factor': 2.**(1./6.)},
    1: {'sigma': 1.0, 'epsilon': 1.0, 'factor': 2.**(1./6.)},
}
DWCA_PARAMETERS = {
    'types': [(1, 1)],
    'rzero': 1.0 * (2.0**(1.0/6.0)),
    'height': 6.0,
    'width': 0.25,
}
# Give simulation settings:
SETTINGS = {
    'system': {
        'temperature': 2.0,
        'units': 'lj',
    },
    'box': {
        'low': [0.0, 0.0],
        'high': [3.6, 3.6],
    },
    'simulation': {
        'task': 'md-nve',
        'steps': 1100,
    },
    'engine': {
        'class': 'velocityverlet',
        'timestep': 0.0025,
    },
    'output': {
        'backup': 'overwrite',
        'energy-file': 1,
        'screen': 10,
        'trajectory-file': 10,
    },
    'potential': [
        {
            'class': 'PairLennardJonesCutnp',
            'dim': 2,
            'shift': True,
            'mixing': 'geometric',
            'parameter': WCA_PARAMETERS,
        },
        {
            'class': 'DoubleWellWCA',
            'dim': 2,
            'parameter': DWCA_PARAMETERS,
        },
    ],
    'particles': {
        'position': {'generate': 'sq', 'repeat': [3, 3], 'lcon': 1.0},
        'velocity': {'generate': 'maxwell', 'momentum': True, 'seed': 0},
        'ptype': [0, 1, 1, 0],
        'name': ['A', 'B', 'B', 'A'],
        'mass': {'A': 1.0, 'B': 1.0},
    },
}
UNIT = SETTINGS['system']['units']
create_conversion_factors(UNIT)
print('# Creating system from settings.')
settings = SETTINGS
system = create_system(settings)
system.forcefield = create_force_field(settings)
system.particles.pos -= (np.average(system.particles.pos, axis=0) -
                         0.5 * system.box.length)  # center in box
print('# Creating simulation from settings.')

settings['system']['obj'] = system
simulation = create_simulation(settings)
print('# Creating output tasks from settings.')
simulation.set_up_output(settings, progress=False)
size = system.box.bounds()
BIDX = [i for i, ptype in enumerate(system.particles.ptype) if ptype == 1]
dwca = system.forcefield.potential[1]
# Some additional set-up for the animation
timeunit = (settings['engine']['timestep'] *
            CONVERT['time'][UNIT, 'fs'])
timeendfs = settings['simulation']['steps'] * timeunit

time, step, v_pot, e_kin, e_tot, temperature = [], [], [], [], [], []
SIGMA = CONVERT['length'][UNIT, 'A']
ECONV = CONVERT['energy'][UNIT, 'kcal/mol']
# We will in this example animate on the fly. Here we do some additional
# set-up to be able to do just that :-)
mpl.rc('axes', labelsize='large')
mpl.rc('font', family='serif')
fig = plt.figure(figsize=(12, 6))
grid = gridspec.GridSpec(2, 2)
# This adds the first axis. Here we will plot the
# particles with velocity and force vectors.
ax1 = fig.add_subplot(grid[:, 0])
ax1.set_xlim((size[0] + np.array([-0.6, 0.6])) * SIGMA)
ax1.set_ylim((size[1] + np.array([-0.6, 0.6])) * SIGMA)
ax1.set_aspect('equal', 'datalim')
ax1.set_xlabel(u'x / Å')
ax1.set_ylabel(u'y / Å')
ax1.set_xticks([size[0][0] * SIGMA, size[0][1] * SIGMA])
ax1.set_yticks([size[1][0] * SIGMA, size[1][1] * SIGMA])
ax1.xaxis.labelpad = -10
ax1.yaxis.labelpad = -10

pos0 = system.box.pbc_wrap(system.particles.pos)
# Set up circles to represent the particles:
circles = []
for _ in range(system.particles.npart):
    circles.append(plt.Circle((0, 0), radius=SIGMA * 0.5, alpha=0.5))
    circles[-1].set_visible(False)
    ax1.add_patch(circles[-1])
# Add arrows for the forces and velocities:
force_arrow = plt.quiver(pos0[:, 0], pos0[:, 1],
                         color=COLORS['almost_black'], zorder=4)
vel_arrow = plt.quiver(pos0[:, 0], pos0[:, 1],
                       color=COLOR_SCHEME['colorblind_10'][1], zorder=4)
# Also add arrows for a "legend":
plt.quiverkey(force_arrow, 3, -4.7, 9, 'Forces', coordinates='data',
              color=COLORS['almost_black'], fontproperties={'size': 'large'})
plt.quiverkey(vel_arrow, 9, -4.7, 9, 'Velocities', coordinates='data',
              color=COLOR_SCHEME['colorblind_10'][1],
              fontproperties={'size': 'large'})
# Also add a line representing the bond
linebond, = ax1.plot([], [], lw=3, ls='-', color=PCOLOR['B'], alpha=0.8)
# Draw lines representing the box boundaries:
ax1.axhline(y=size[1][0] * SIGMA, lw=4, ls=':', alpha=0.5,
            color=COLORS['almost_black'])
ax1.axhline(y=size[1][1] * SIGMA, lw=4, ls=':', alpha=0.5,
            color=COLORS['almost_black'])
ax1.axvline(x=size[0][0] * SIGMA, lw=4, ls=':', alpha=0.5,
            color=COLORS['almost_black'])
ax1.axvline(x=size[0][1] * SIGMA, lw=4, ls=':', alpha=0.5,
            color=COLORS['almost_black'])
# Add second axis for displaying energies:
ax2 = fig.add_subplot(grid[0, 1])
ax2.set_xlim(0, timeendfs)
ax2.set_ylim(-0.6, 1.1)
ax2.set_xlabel('Time / fs')
ax2.set_ylabel('Energy / (kcal/mol)')
time_text = ax2.text(0.02, 0.90, '', transform=ax2.transAxes)
linepot, = ax2.plot([], [], lw=4, ls='-', color=COLOR_SCHEME['deep'][0],
                    alpha=0.8, label='Potential')
linekin, = ax2.plot([], [], lw=4, ls='-', color=COLOR_SCHEME['deep'][1],
                    alpha=0.8, label='Kinetic')
linetot, = ax2.plot([], [], lw=4, ls='-', color=COLORS['almost_black'],
                    alpha=0.8, label='Total')
ax2.legend(loc='lower left', ncol=3, frameon=False,
           columnspacing=1, labelspacing=1)


def plot_dwca_potential():
    """Plot the double well WCA potential.

    This is a helper function to plot the potential.
    It is creating a fake system with a fake box and moves a
    particle relative to another one in order to obtain the potential.

    Returns
    -------
    out[0] : numpy.array
        Positions, can be used as an x-coordinate in a plot.
    out[1] : numpy.array
        The potential energy as a function of `out[0]`, can be used as the
        y-coordinate in a plot.

    """
    rpos = np.linspace(0.1, 5, 500)
    potdwca = []
    fakebox = create_box(low=[0.0, 0.0], high=[10.0, 10.0])
    fakesys = System(units='lj', box=fakebox)
    fakesys.particles = Particles(dim=system.get_dim())
    fakesys.add_particle(name='B', pos=np.zeros(2), ptype=1)
    fakesys.add_particle(name='B', pos=np.zeros(2), ptype=1)
    for ri in rpos:
        fakesys.particles.pos[-1] = np.array([ri, 0.0])
        potdwca.append(dwca.potential(fakesys))
    return rpos, np.array(potdwca)


# add third axis for plotting the potential and order parameter:
ax3 = fig.add_subplot(grid[1, 1])
rbond, pot_dwca = plot_dwca_potential()
linedwpot, = ax3.plot(rbond, pot_dwca, lw=3, ls='-',
                      color=COLORS['almost_black'])
ax3.set_ylim(0, dwca.params['height'] + 1)
min_max = dwca.min_max()
ax3.set_xlim(min_max[0] - 0.2, min_max[1] + 0.2)
ax3.set_ylabel('Double well potential')
ax3.set_xlabel('Bond length (circle)')
orderscatter = ax3.scatter([], [], alpha=0.5, color=PCOLOR['B'],
                           marker='o', s=200)
plt.subplots_adjust(left=0.08, top=0.95, bottom=0.15, right=0.95, hspace=0.3)


def init():
    """Declare what to re-draw when clearing the animation frame.

    Returns
    -------
    out : list
        list of the patches to be drawn.

    """
    patches = []
    for ci in circles:
        ci.set_visible(False)
        patches.append(ci)
    linebond.set_data([], [])
    patches.append(linebond)
    force_arrow.set_visible(False)
    patches.append(force_arrow)
    vel_arrow.set_visible(False)
    patches.append(vel_arrow)
    time_text.set_text('')
    patches.append(time_text)
    orderscatter.set_offsets(([], []))
    patches.append(orderscatter)
    return patches


def get_max_vector(vectors):
    """Determine the longest vector in a list of vectors.

    Parameters
    ----------
    vectors : numpy.array
        Numpy array of vectors to analyse, using `numpy.dot`.

    Returns
    -------
    vmax : float
        The length of the largest vector.

    """
    vmax = None
    for vi in vectors:
        vs = np.sqrt(np.dot(vi, vi))
        if vmax is None or vs > vmax:
            vmax = vs
    return vmax


def get_velocity_force_arrows(forces, vels):
    """Obtain the force and velocity vectors.

    Parameters
    ----------
    forces : numpy.array
        The forces on the particles.
    vels : numpy.array
        The velocity of the particles.

    Returns
    -------
    out[0] : numpy.array
        The x-component of the forces, normalised.
    out[1] : numpy.array
        The y-component of the forces, normalised.
    out[2] : numpy.array
        The x-component of the velocities, normalised.
    out[3] : numpy.array
        The y-component of the velocities, normalised.

    """
    fmax, vmax = get_max_vector(forces), get_max_vector(vels)
    FU, FV, VU, VV = [], [], [], []
    for fi, vi in zip(forces, vels):
        fj = 10.0 * fi / fmax
        vj = 10.0 * vi / vmax
        FU.append(fj[0])
        FV.append(fj[1])
        VU.append(vj[0])
        VV.append(vj[1])
    return FU, FV, VU, VV


def spring_bond(delta, dr, part1, part2):
    """Create positions for a zig-zag line.

    This is a function that will create positions which can be used to
    create a zig-zag bond. It is used here to illustrate a spring bond
    between two atoms

    Parameters
    ----------
    delta : numpy.array
        Distance vector between `part1` and `part2`, subjected to periodic
        boundary conditions.
    dr : float
        Length of `delta` vector.
    part1 : numpy.array
        Particle 1 position. Bond is drawn from `part1`.
    part2 : numpy.array
        Particle 2 position. Bond is drawn to `part2`.

    Returns
    -------
    out[0] : numpy.array
        X-coordinates for the line.
    out[1] : numpy.array
        Y-coordinates for the line.

    """
    delta_u = delta / dr
    xpos = [part1[0]]
    ypos = [part1[1]]
    for pidx, add in enumerate(np.linspace(0.0, dr-1, 11)):
        point = part1 + (add + 0.5) * delta_u
        if pidx in [2, 4, 6, 8]:
            if delta_u[0] == 0:
                dperp = np.array([0.0, 0.0])
            else:
                dperp_v = np.array([-delta_u[1] / delta_u[0], 1.0])
                dperp = dperp_v / np.sqrt(np.dot(dperp_v, dperp_v))
            sig = 1 if delta_u[0] > 0.0 else -1.0
            if pidx in [2, 6]:
                dvec = sig*0.2*dperp
            else:
                dvec = -sig*0.2*dperp
            point = point + dvec
        xpos.append(point[0])
        ypos.append(point[1])
    xpos.append(part2[0])
    ypos.append(part2[1])
    xpos = np.array(xpos) * SIGMA
    ypos = np.array(ypos) * SIGMA
    return xpos, ypos


def update(frame, sys, sim):
    """Update plots for the animation.

    This function will be running the simulation and updating the plots.
    It is called one time per step, and we choose to update the simulation
    inside this function

    Parameters
    ----------
    frame : int
        The current frame number, supplied by `animation.FuncAnimation`.
    sys : object
        The system object we are simulating
    sim : object
        The simulation we are running

    Returns
    -------
    out : list
        List of the patches to be drawn.

    """
    pos = sys.particles.pos
    patches = []
    for ci, pi, itype in zip(circles, pos, sys.particles.name):
        ci.center = np.array([pi[0], pi[1]]) * SIGMA
        ci.set_color(PCOLOR[itype])
        ci.set_visible(True)
        patches.append(ci)
    # update the force and velocity vectors:
    FU, FV, VU, VV = get_velocity_force_arrows(sys.particles.force,
                                               sys.particles.vel)
    force_arrow.set_offsets(pos * SIGMA)
    force_arrow.set_UVC(FU, FV)
    force_arrow.set_visible(True)
    patches.append(force_arrow)
    vel_arrow.set_offsets(pos * SIGMA)
    vel_arrow.set_UVC(VU, VV)
    vel_arrow.set_visible(True)
    patches.append(vel_arrow)

    if not sim.is_finished():
        result = sim.step()
        # Since we are not using run - manually output:
        for task in sim.output_tasks:
            task.output(result)
        # reaction coordinate:
        delta = sys.box.pbc_dist_coordinate(sys.particles.pos[BIDX[1]] -
                                            sys.particles.pos[BIDX[0]])
        dr = np.sqrt(np.dot(delta, delta))
        points = [dr, dwca.potential(sys)]
        orderscatter.set_offsets(points)
        patches.append(orderscatter)
        # draw the bond:
        xpos, ypos = spring_bond(delta, dr,
                                 sys.particles.pos[BIDX[0]],
                                 sys.particles.pos[BIDX[1]])
        linebond.set_data(xpos, ypos)
        patches.append(linebond)
        # here we calculate some energies and updates the energy plots:
        step.append(result['cycle']['step'])
        time.append(result['cycle']['step'] * timeunit)
        temperature.append(result['thermo']['temp'])
        v_pot.append(ECONV * result['thermo']['vpot'])
        e_kin.append(ECONV * result['thermo']['ekin'])
        e_tot.append(e_kin[-1] + v_pot[-1])
        # update plots with energies:
        linepot.set_data(time, (np.array(v_pot) - v_pot[0]))
        patches.append(linepot)
        linekin.set_data(time, (np.array(e_kin)))
        patches.append(linekin)
        linetot.set_data(time, (np.array(e_tot) - v_pot[0]))
        patches.append(linetot)
        # also display current simulation time;
        time_text.set_text(f'Time: {time[-1]:6.2f} fs (frame: {frame})')
        patches.append(time_text)
        return patches
    print('Simulation is done')
    return patches


# This will run the animation/simulation:
anim = animation.FuncAnimation(fig, update,
                               frames=settings['simulation']['steps']+1,
                               fargs=[system, simulation],
                               repeat=False, interval=2, blit=True,
                               init_func=init)
# for making a movie:
# anim.save('particles.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
plt.show()

# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Example of running a MD NVE simulation.

In this example, we animate the output.
"""
# pylint: disable=invalid-name
import numpy as np
# imports for the plotting:
from matplotlib import pyplot as plt
from matplotlib import animation
import matplotlib as mpl
from pyretis.core.units import CONVERT, create_conversion_factors
from pyretis.inout.plotting import COLORS, COLOR_SCHEME
from pyretis.setup import create_simulation

ljparams = {0: {'sigma': 1.0, 'epsilon': 1.0, 'rcut': 2.5}}
# Define simulation settings:
settings = {}
settings['system'] = {
    'temperature': 1.0,
    'dimensions': 2,
    'units': 'lj'
}
settings['box'] = {
    'low': [0.0, 0.0],
    'high': [1.1*3.405, 1.1*3.405],
    'periodic': [True, True]
}
settings['simulation'] = {
    'task': 'md-nve',
    'steps': 950
}
settings['engine'] = {
    'class': 'velocityverlet',
    'timestep': 0.0025
}
settings['output'] = {
    'backup': 'overwrite',
    'energy-file': 1,
    'screen': 10,
    'trajectory-file': 10
}
settings['potential'] = [
    {'class': 'PairLennardJonesCutnp',
     'dim': 2, 'shift': True,
     'mixing': 'geometric',
     'parameter': ljparams}
]
settings['particles'] = {
    'position': {'input_file': 'initial.xyz'},
    'velocity': {'generate': 'maxwell', 'momentum': True, 'seed': 0}
}
UNIT = settings['system']['units']
create_conversion_factors(UNIT)
SIGMA = CONVERT['length'][UNIT, 'A']
ECONV = CONVERT['energy'][UNIT, 'kcal/mol']
print('# Creating system from settings.')
simulation_nve = create_simulation(settings)
print('# Creating output tasks from settings.')
simulation_nve.set_up_output(settings, progress=False)
ljsystem = simulation_nve.system
size = ljsystem.box.bounds()
npart = ljsystem.particles.npart
msg = 'Added {:d} particles to a simple square lattice'
print(msg.format(npart))
npart = float(npart)

# We will in this example animate on the fly.
# We will then have to do some extra set up. The actual simulation is
# carried out by calling `simulation_nve.step()` in the `update`
# function which is executed by the animation.FuncAnimation() call.
# In effect animation.FuncAnimation will run the simulation one step,
# update the plot and display it and continue this loop until the
# simulation is done.
timeunit = (settings['engine']['timestep'] *
            CONVERT['time'][UNIT, 'fs'])
timeendfs = settings['simulation']['steps'] * timeunit
time, step, v_pot, e_kin, e_tot, temperature = [], [], [], [], [], []
mpl.rc('axes', labelsize='large')
mpl.rc('font', family='serif')
fig = plt.figure(figsize=(12, 6))
# This adds the first axis. Here we will plot the
# particles with velocity and force vectors.
ax1 = fig.add_subplot(121)
ax1.set_xlim((size[0] + np.array([-0.5, 0.5])) * SIGMA)
ax1.set_ylim((size[1] + np.array([-0.5, 0.5])) * SIGMA)
ax1.set_aspect('equal', 'datalim')
ax1.set_xlabel(u'x / Å')
ax1.set_ylabel(u'y / Å')
ax1.set_xticks([size[0][0] * SIGMA, size[0][1] * SIGMA])
ax1.set_yticks([size[1][0] * SIGMA, size[1][1] * SIGMA])
ax1.xaxis.labelpad = -5
ax1.yaxis.labelpad = -5

pos0 = ljsystem.box.pbc_wrap(ljsystem.particles.pos)
# Set up circles to represent the particles:
circles = []
for _ in range(int(npart)):
    circles.append(plt.Circle((0, 0), radius=SIGMA * 0.5, alpha=0.5,
                              color='blue'))
    circles[-1].set_visible(False)
    ax1.add_patch(circles[-1])
# Add arrows for the forces and velocities:
force_arrow = plt.quiver(pos0[:, 0], pos0[:, 1],
                         color=COLORS['almost_black'], zorder=4)
vel_arrow = plt.quiver(pos0[:, 0], pos0[:, 1],
                       color=COLOR_SCHEME['colorblind_10'][1], zorder=4)
# Also add arrows for a "legend":
plt.quiverkey(force_arrow, 3, -3.5, 9, 'Forces', coordinates='data',
              color=COLORS['almost_black'], fontproperties={'size': 'large'})
plt.quiverkey(vel_arrow, 9, -3.5, 9, 'Velocities', coordinates='data',
              color=COLOR_SCHEME['colorblind_10'][1],
              fontproperties={'size': 'large'})
# Draw the lines representing the box boundaries:
ax1.axhline(y=size[1][0] * SIGMA, lw=4, ls=':', alpha=0.5,
            color=COLORS['almost_black'])
ax1.axhline(y=size[1][1] * SIGMA, lw=4, ls=':', alpha=0.5,
            color=COLORS['almost_black'])
ax1.axvline(x=size[0][0] * SIGMA, lw=4, ls=':', alpha=0.5,
            color=COLORS['almost_black'])
ax1.axvline(x=size[0][1] * SIGMA, lw=4, ls=':', alpha=0.5,
            color=COLORS['almost_black'])
# Add second axis for plotting the energies
ax2 = fig.add_subplot(122)
ax2.set_xlim(0, timeendfs)
ax2.set_ylim(-0.05, 0.25)
ax2.set_xlabel('Time / fs')
ax2.set_ylabel('Energy / (kcal/mol)')
time_text = ax2.text(0.02, 0.95, '', transform=ax2.transAxes)
linepot, = ax2.plot([], [], lw=4, ls='-', color=COLOR_SCHEME['deep'][0],
                    alpha=0.8, label='Potential')
linekin, = ax2.plot([], [], lw=4, ls='-', color=COLOR_SCHEME['deep'][1],
                    alpha=0.8, label='Kinetic')
linetot, = ax2.plot([], [], lw=4, ls='-', color=COLORS['almost_black'],
                    alpha=0.8, label='Total')
ax2.legend(loc='lower left', ncol=2, frameon=False)
plt.subplots_adjust(left=0.08, right=0.97, top=0.95, bottom=0.15)


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
    forceu, forcev, velu, velv = [], [], [], []
    for fi, vi in zip(forces, vels):
        fj = 10.0 * fi / fmax
        vj = 10.0 * vi / vmax
        forceu.append(fj[0])
        forcev.append(fj[1])
        velu.append(vj[0])
        velv.append(vj[1])
    return forceu, forcev, velu, velv


def update(frame, sim):
    """Update plots for the animation.

    This function will be running the simulation and updating the plots.
    It is called one time per step, and we choose to update the simulation
    inside this function

    Parameters
    ----------
    frame : int
        The current frame number, supplied by `animation.FuncAnimation`.
    sim : object like `Simulation`
        The simulation we are running.

    Returns
    -------
    out : list
        list of the patches to be drawn.

    """
    particles = sim.system.particles
    pos = sim.system.box.pbc_wrap(particles.pos)
    patches = []
    # update positions of the circles according to the particles:
    for ci, pi in zip(circles, pos):
        ci.center = np.array([pi[0], pi[1]]) * SIGMA
        ci.set_visible(True)
        patches.append(ci)
    # update the force and velocity vectors:
    forceu, forcev, velu, velv = get_velocity_force_arrows(particles.force,
                                                           particles.vel)
    force_arrow.set_offsets(pos * SIGMA)
    force_arrow.set_UVC(forceu, forcev)
    force_arrow.set_visible(True)
    patches.append(force_arrow)
    vel_arrow.set_offsets(pos * SIGMA)
    vel_arrow.set_UVC(velu, velv)
    vel_arrow.set_visible(True)
    patches.append(vel_arrow)

    if not sim.is_finished():
        result = sim.step()
        for tsk in sim.output_tasks:
            tsk.output(result)
        # here we calculate some energies and updates the energy plots:
        step.append(result['cycle']['step'])
        time.append(step[-1] * timeunit)
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
    print('Simulation is done.')
    return patches


def init():
    """Declare what to re-draw when clearing the animation frame.

    Returns
    -------
    out : list
        list of the patches to be drawn.

    """
    patches = []
    force_arrow.set_visible(False)
    patches.append(force_arrow)
    vel_arrow.set_visible(False)
    patches.append(vel_arrow)
    time_text.set_text('')
    patches.append(time_text)
    for ci in circles:
        ci.set_visible(False)
        patches.append(ci)
    return patches


# This will run the animation/simulation:
anim = animation.FuncAnimation(fig, update,
                               frames=settings['simulation']['steps']+1,
                               fargs=[simulation_nve],
                               repeat=False, interval=2, blit=True,
                               init_func=init)
# for making a movie:
# anim.save('particles.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
plt.show()

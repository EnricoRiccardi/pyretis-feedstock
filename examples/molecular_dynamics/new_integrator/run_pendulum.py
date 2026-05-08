# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Double pendulum example."""
# pylint: disable=invalid-name
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib import gridspec
from pyretis.core.units import create_conversion_factors
from pyretis.setup import (create_system, create_engine,
                           create_force_field, create_simulation)

# Define potential parameters:
L1 = 1.5
L2 = 2.0
M1 = 1.0
M2 = 1.0
L1SQ = L1**2
L2SQ = L2**2
THETA1 = 45. * np.pi / 180.
THETA2 = 90. * np.pi / 180.
DTHETA1 = 0.0
DTHETA2 = 0.0
STEPS_PER_FRAME = 10
# Give simulation settings:
settings = {}
settings['simulation'] = {
    'task': 'md-nve',
    'steps': 100000,
    'exe_path': None
}
settings['system'] = {
    'units': 'gromacs',
    'temperature': 1.,
    'dimensions': 2
}
settings['engine'] = {
    'class': 'VVIntegrator',
    'timestep': 0.005,
    'module': 'myintegrator.py'
}
settings['output'] = {
    'backup': 'overwrite',
    'energy-file': 1,
    'screen': 10,
    'trajectory-file': 10
}
settings['potential'] = [
    {'class': 'DoublePendulum',
     'module': 'doublependulum.py',
     'g': 1.0, 'l1': L1, 'l2': L2}
]
settings['particles'] = {
    'position': {'generate': 'sq', 'repeat': [2, 1], 'lcon': 1},
    'ptype': [0, 0],
    'name': ['A', 'B'],
    'mass': {'A': M1, 'B': M2}
}

UNIT = settings['system']['units']
create_conversion_factors(UNIT)
print('# Creating system from settings.')
system = create_system(settings)
system.forcefield = create_force_field(settings)
system.particles.pos[0][0] = THETA1
system.particles.pos[1][0] = THETA2
system.particles.vel[0][0] = DTHETA1
system.particles.vel[1][0] = DTHETA2
settings['system']['obj'] = system
simulation = create_simulation(settings)
mpl.rc('axes', labelsize='large')
mpl.rc('font', family='serif')
fig = plt.figure(figsize=(12, 6))
grid = gridspec.GridSpec(3, 2)
ax1 = fig.add_subplot(grid[:, 0])
ax1.set_xlim((-5, 5))
ax1.set_ylim((-5, 2))
ax1.set_aspect('equal')
time_text = ax1.text(0.02, 0.90, '', transform=ax1.transAxes)

# set up circles to represent the two objects:
circ1 = plt.Circle((0, 0), radius=0.25, alpha=0.6)
circ1.set_visible(False)
ax1.add_patch(circ1)
circ2 = plt.Circle((0, 0), radius=0.25, alpha=0.6)
circ2.set_visible(False)
ax1.add_patch(circ2)
# Add one for the orrigin
circ0 = plt.Circle((0, 0), radius=0.1, alpha=0.6, color='#262626')
ax1.add_patch(circ0)
# Add a lines for the connection between them
line1, = ax1.plot([], [], lw=3, ls='-', color='#262626', alpha=0.5)
line2, = ax1.plot([], [], lw=3, ls='-', color='#262626', alpha=0.5)
# add second axis for displaying energies:
ax2 = fig.add_subplot(grid[0, 1])
linepot, = ax2.plot([], [], lw=3, ls='-', alpha=0.8)
ax2.set_ylim(-7.5, 0)
ax3 = fig.add_subplot(grid[1, 1])
linekin, = ax3.plot([], [], lw=3, ls='-', alpha=0.8)
ax3.set_ylim(0, 3.5)
ax4 = fig.add_subplot(grid[2, 1])
linetot, = ax4.plot([], [], lw=3, ls='-', alpha=0.8)
ax4.set_ylim(-2.5, -0.5)
ax4.set_xlabel('Steps')
ax2.set_xticks([])
ax3.set_xticks([])


def init():
    """Declare what to re-draw when clearing the animation frame.

    Returns
    -------
    out : list
        List of the patches to be drawn.

    """
    patches = []
    circ1.set_visible(False)
    circ2.set_visible(False)
    patches.append(circ1)
    patches.append(circ2)
    line1.set_data([], [])
    patches.append(line1)
    line2.set_data([], [])
    patches.append(line2)
    time_text.set_text('')
    patches.append(time_text)
    return patches


def update(frame, sys, sim):
    """Update plots for the animation.

    This function will be running the simulation and updating the plots.
    It is called one time per step, and we choose to update the simulation
    inside this function

    Parameters
    ----------
    frame : int
        The current frame number, supplied by `animation.FuncAnimation`.
    sys : object like :py:class:`.System`
        The system object we are simulating.
    sim : object like :py:class:`.Simulation`
        The simulation we are running.

    Returns
    -------
    out : list
        List of the patches to be drawn.

    """
    patches = []
    if not sim.is_finished():
        for _ in range(STEPS_PER_FRAME):
            result = sim.step()
        pos = sys.particles.pos
        theta1 = pos[0][0]
        costheta1 = np.cos(theta1)
        sintheta1 = np.sin(theta1)
        x1, y1 = L1 * sintheta1, -L1 * costheta1
        circ1.center = np.array([x1, y1])
        circ1.set_visible(True)
        patches.append(circ1)
        theta2 = pos[1][0]
        costheta2 = np.cos(theta2)
        sintheta2 = np.sin(theta2)
        x2, y2 = x1 + L2*sintheta2, y1 - L2 * costheta2
        circ2.center = np.array([x2, y2])
        circ2.set_visible(True)
        patches.append(circ2)
        line1.set_data([0.0, x1], [0.0, y1])
        patches.append(line1)
        line2.set_data([x1, x2], [y1, y2])
        patches.append(line2)
        vpot = result['thermo']['vpot']
        vel = sys.particles.vel
        dtheta1 = vel[0][0]
        dtheta2 = vel[1][0]
        dt1sq = dtheta1**2
        dt2sq = dtheta2**2
        ekin = (0.5*(M1+M2)*L1SQ*dt1sq + 0.5*M2*L2SQ*dt2sq +
                M2*L1*L2*dtheta1*dtheta2*np.cos(theta1-theta2))
        vpot = 2.0 * result['thermo']['vpot']
        etot = ekin + vpot
        linekin.set_xdata(np.append(linekin.get_xdata(), frame))
        linekin.set_ydata(np.append(linekin.get_ydata(), ekin))
        patches.append(linekin)
        linepot.set_xdata(np.append(linepot.get_xdata(), frame))
        linepot.set_ydata(np.append(linepot.get_ydata(), vpot))
        patches.append(linepot)
        linetot.set_xdata(np.append(linetot.get_xdata(), frame))
        linetot.set_ydata(np.append(linetot.get_ydata(), etot))
        patches.append(linetot)
        ax2.set_xlim(0, frame + 1)
        ax3.set_xlim(0, frame + 1)
        ax4.set_xlim(0, frame + 1)
        time_text.set_text(f"Step: {result['cycle']['step']}")
        patches.append(time_text)
        return patches
    print('Simulation is done')
    return patches


# This will run the animation/simulation:
FMAX = int(float(settings['simulation']['steps']) / float(STEPS_PER_FRAME))
anim = animation.FuncAnimation(fig, update,
                               frames=FMAX+1,
                               fargs=[system, simulation],
                               repeat=False, interval=2, blit=True,
                               init_func=init)
# for making a movie:
# anim.save('particles.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
plt.show()

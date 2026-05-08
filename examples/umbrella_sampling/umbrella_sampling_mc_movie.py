# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""
This is an example of how we can create an animation with PyRETIS.

It will run the umbrella sampling defined in `umbrella_sampling_mc.py` and
draw the frames as an animation.

"""
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from umbrella_sampling_mc import (
    UMBRELLA_WINDOWS,
    set_up_simulation,
    calculate_unbiased_potential,
)


def main_amination(windows, save=False):
    """Run the animation."""
    # For plotting unbiased potential:
    xpos = np.linspace(-2, 2, 250)
    vpot = calculate_unbiased_potential(xpos)

    fig = plt.figure()
    axs = plt.axes(xlim=(-1.05, 1.05), ylim=(-0.3, 0.05))
    axs.set_xlabel(r'Position ($x$)', fontsize='large')
    axs.set_ylabel('Potential energy ($V(x)$) / eV', fontsize='large')
    # Plot the line we have to cross:
    linec = axs.axvline(x=-10, lw=2, ls=':', color='#262626')
    # Plot the umbrella region:
    axv = axs.axvspan(xmin=-10, xmax=-10, alpha=0.1)
    # Plot the potential:
    axs.plot(xpos, vpot, lw=3)
    # Plot the particle:
    scat = axs.scatter(-10, -10, s=150)
    # Add text with info:
    txt = axs.text(0.02, 0.95, '', transform=axs.transAxes,
                   fontsize='large')

    def init():
        """Return what we have to re-draw."""
        return [scat, txt]

    n_umb = len(windows)
    simulations = []
    systems = []
    for i, window in enumerate(windows):
        over = windows[min(i + 1, n_umb - 1)][0]
        simulation = set_up_simulation(
            window, over, 1, initial_positions=np.array([window[0]]),
        )
        simulations.append(simulation)
        systems.append(simulation.system)

    anim = animation.FuncAnimation(
        fig,
        update_animation,
        fargs=[simulations, systems, scat, linec, axv, txt],
        repeat=False,
        interval=10,
        blit=True,
        init_func=init
    )
    if save:
        anim.save('movie.mp4')
    else:
        plt.show(block=True)


def update_animation(frame, simulations, systems, scat, linec, axv, txt):
    """Update the animation."""
    patches = []
    finished = 0
    for simulation in simulations:
        if simulation.is_finished():
            finished += 1
    if finished == len(simulations):
        print(f'Stopping animation at frame: {frame}')
        raise StopIteration
    simulation = simulations[finished]
    system = systems[finished]

    if simulation.cycle['step'] == 0:
        over = simulation.overlap
        window = simulation.umbrella
        linec.set_data([over, over], linec.get_data()[1])
        region = np.array([[window[0], 0.0], [window[0], 1.0],
                           [window[1], 1.0], [window[1], 0.0],
                           [window[0], 0.0]])
        axv.set_xy(region)
    patches.append(linec)
    patches.append(axv)
    for _ in range(50):
        if not simulation.is_finished():
            simulation.step()
    if simulation.is_finished() and finished + 1 == len(simulations):
        txt.set_text('Done!')
        patches.append(txt)
    else:
        txt.set_text(
            f"Window {finished}, step {simulation.cycle['step']}"
        )
    patches.append(txt)
    ener = system.particles.vpot
    offsets = [(xi, ener) for xi in system.particles.pos]
    scat.set_offsets(offsets)
    patches.append(scat)
    return patches


if __name__ == '__main__':
    main_amination(UMBRELLA_WINDOWS, save=False)

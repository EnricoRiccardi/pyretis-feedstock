# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Example of performing particle swarm optimization."""
# pylint: disable=invalid-name
import numpy as np
from pyretis.core import create_box, Particles, System
from pyretis.simulation import Simulation
from pyretis.forcefield import ForceField
from psoengine import PSOEngine
from ackley import Ackley, ackley_potential
from matplotlib import pyplot as plt
from matplotlib import animation, cm


NPART = 10
STEPS = 1000
MINX, MAXX = -10, 10
TXT = 'Step: {:5d}: Best: (x, y) = ({:10.3e}, {:10.3e}), pot = {:10.3e}'


def set_up():
    """Just create system and simulation."""
    box = create_box(low=[MINX, MINX], high=[MAXX, MAXX],
                     periodic=[True, True])
    print('Created a box:')
    print(box)

    print(f'Creating system with {NPART} particles')
    system = System(units='reduced', box=box)
    system.particles = Particles(dim=2)
    for _ in range(NPART):
        pos = np.random.uniform(low=MINX, high=MAXX, size=(1, 2))
        system.add_particle(pos)

    ffield = ForceField('Single Ackley function',
                        potential=[Ackley()])
    system.forcefield = ffield
    print(f'Force field is:\n{system.forcefield}')

    print('Creating simulation:')
    engine = PSOEngine(0.7, 1.5, 1.5)
    simulation = Simulation(settings=None, controls={'steps': STEPS})
    task_integrate = {'func': engine.integration_step,
                      'args': [system],
                      'result': 'gbest', 'first': True}
    simulation.add_task(task_integrate)
    return simulation, system


def main():
    """Just run the optimization, no plotting."""
    simulation, _ = set_up()
    for result in simulation.run():
        step = result['cycle']['step']
        best = result['gbest']
        if step % 10 == 0:
            print(TXT.format(step, best[1][0], best[1][1], best[0]))


def evaluate_potential_grid():
    """Evaluate the Ackley potential on a grid."""
    X, Y = np.meshgrid(np.linspace(MINX, MAXX, 100),
                       np.linspace(MINX, MAXX, 100))
    Z = ackley_potential(X, Y)
    return X, Y, Z


def update_animation(frame, system, simulation, scatter):
    """Update animation."""
    patches = []
    if not simulation.is_finished() and frame > 0:
        results = simulation.step()
        best = results['gbest']
        if frame % 10 == 0:
            print(TXT.format(frame, best[1][0], best[1][1], best[0]))
    scatter.set_offsets(system.particles.pos)
    patches.append(scatter)
    return patches


def main_animation():
    """Run the simulation and update for animation."""
    simulation, system = set_up()
    fig = plt.figure()
    ax1 = fig.add_subplot(111, aspect='equal')
    ax1.set_xlim((MINX, MAXX))
    ax1.set_ylim((MINX, MAXX))
    ax1.axes.get_xaxis().set_visible(False)
    ax1.axes.get_yaxis().set_visible(False)
    X, Y, pot = evaluate_potential_grid()
    ax1.contourf(X, Y, pot, cmap=cm.viridis, zorder=1)
    scatter = ax1.scatter(system.particles.pos[:, 0],
                          system.particles.pos[:, 1], marker='o', s=50,
                          edgecolor='#262626', facecolor='white')

    def init():
        """Just return what to re-draw."""
        return [scatter]
    # This will run the animation/simulation:
    anim = animation.FuncAnimation(fig, update_animation,
                                   frames=STEPS+1,
                                   fargs=[system, simulation, scatter],
                                   repeat=False, interval=30, blit=True,
                                   init_func=init)
    plt.show()
    return anim


if __name__ == '__main__':
    main_animation()

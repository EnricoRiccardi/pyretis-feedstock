# -*- coding: utf-8 -*-
# Copyright (c) 2023, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Definition of numerical MD integrators.

These integrators are representations of engines for performing
molecular dynamics. Typically they will propagate Newtons
equations of motion in time numerically.

Important classes defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MDEngine (:py:class:`.MDEngine`)
    Base class for internal MDEngines.

RandomWalk (:py:class:`.RandomWalk`)
    A Random Walk integrator.

Verlet (:py:class:`.Verlet`)
    A Verlet MD integrator.

VelocityVerlet (:py:class:`.VelocityVerlet`)
    A Velocity Verlet MD integrator.

Langevin (:py:class:`.Langevin`)
    A Langevin MD integrator.
"""
import logging
import numpy as np
from pyretis.core.random_gen import create_random_generator
from pyretis.core.particlefunctions import (calculate_kinetic_energy,
                                            calculate_thermo,
                                            calculate_thermo_path,
                                            reset_momentum)
from pyretis.engines.engine import EngineBase

logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


__all__ = ['MDEngine', 'RandomWalk', 'Verlet', 'VelocityVerlet', 'Langevin']


class MDEngine(EngineBase):
    """Base class for internal MD integrators.

    This class defines an internal MD integrator. This class of
    integrators work with the positions and velocities of the system
    object directly. Further, we make use of the system object in
    order to update forces etc.

    Attributes
    ----------
    timestep : float
        Time step for the integration.
    description : string
        Description of the MD integrator.
    dynamics : str
        A short string to represent the type of dynamics produced
        by the integrator (NVE, NVT, stochastic, ...).

    """

    engine_type = 'internal'

    def __init__(self, timestep, description, dynamics=None):
        """Set up the integrator.

        Parameters
        ----------
        timestep : float
            The time step for the integrator in internal units.
        description : string
            A short description of the integrator.
        dynamics : string or None, optional
            Description of the kind of dynamics the integrator does.

        """
        super().__init__(description)
        self.timestep = timestep
        self.dynamics = dynamics

    def integration_step(self, _):
        """Perform a single time step of the integration.

        Parameters
        ----------
        - : place holder

        Returns
        -------
        out : None
            Does not return anything, in derived classes it will
            typically update the given `System`.

        """
        raise NotImplementedError

    def select_thermo_function(self, thermo='full'):
        """Select function for calculating thermodynamic properties.

        Parameters
        ----------
        thermo : string, or None, optional
            String which selects the kind of thermodynamic output.

        Returns
        -------
        thermo_func : callable or None
            The function matching the requested thermodynamic output.

        """
        thermo_func = None
        if thermo is not None:
            if thermo not in ('full', 'path'):
                logger.debug(
                    'Unknown thermo "%s" requested: Using "path"!',
                    thermo,
                )
                thermo = 'path'
            logger.debug('Thermo output for %s.integrate(): "%s"',
                         self.__class__.__name__, thermo)
            if thermo == 'full':
                thermo_func = calculate_thermo
            elif thermo == 'path':
                thermo_func = calculate_thermo_path
        return thermo_func

    def integrate(self, ensemble, steps, thermo='full'):
        """Perform several integration steps.

        This method will perform several integration steps, but it will
        also calculate order parameter(s) if requested and energy
        terms.

        Parameters
        ----------
        ensemble: dict
            It contains:

            * `system` : object like :py:class:`.System`
              The system we are integrating.
            * `order_function` : object like :py:class:`.OrderParameter`
              An order function can be specified if we want to
              calculate the order parameter along with the simulation.

        steps : integer
            The number of steps we are going to perform. Note that we
            do not integrate on the first step (e.g. step 0) but we do
            obtain the other properties. This is to output the starting
            configuration.
        thermo : string, optional
            Select the thermodynamic properties we are to calculate.

        Yields
        ------
        results : dict
            The result of a MD step. This contains the state of the
            system and also the order parameter(s) (if calculated) and
            the thermodynamic quantities (if calculated).

        """
        thermo_func = self.select_thermo_function(thermo=thermo)
        system = ensemble['system']
        order_function = ensemble.get('order_function')
        system.potential_and_force()
        for i in range(steps):
            if i == 0:
                pass
            else:
                self.integration_step(system)
            results = {'system': system}
            if order_function is not None:
                results['order'] = self.calculate_order(ensemble)
            if thermo_func:
                results['thermo'] = thermo_func(system)
            yield results

    def invert_dt(self):
        """Invert the time step for the integration.

        Returns
        -------
        out : boolean
            True if the time step is positive, False otherwise.

        """
        self.timestep *= -1.0
        return self.timestep > 0.0

    def propagate(self, path, ensemble, reverse=False):
        """Generate a path by integrating until a criterion is met.

        This function will generate a path by calling the function
        specifying the integration step repeatedly. The integration is
        carried out until the order parameter has passed the specified
        interfaces or if we have integrated for more than a specified
        maximum number of steps. The given system defines the initial
        state and the system is reset to its initial state when this
        method is done.

        Parameters
        ----------
        path : object like :py:class:`.PathBase`
            This is the path we use to fill in phase-space points.
            We are here not returning a new path, this since we want
            to delegate the creation of the path (type) to the method
            that is running `propagate`.
        ensemble: dict
            It contains:

            * `system` : object like :py:class:`.System`
              The system object gives the initial state for the
              integration. The initial state is stored and the system is
              reset to the initial state when the integration is done.
            * `order_function` : object like :py:class:`.OrderParameter`
              An order function can be specified if we want to
              calculate the order parameter along with the simulation.
            * `interfaces` : list of floats
              These interfaces define the stopping criterion.

        reverse : boolean, optional
            If True, the system will be propagated backward in time.

        Returns
        -------
        success : boolean
            This is True if we generated an acceptable path.
        status : string
            A text description of the current status of the propagation.

        """
        status = 'Propagate w/internal engine'
        logger.debug(status)
        success = False
        # Copy the system, so that we can propagate without altering it:
        system = ensemble['system']
        system.potential_and_force()  # Make sure forces are set.
        left, _, right = ensemble['interfaces']
        for i in range(path.maxlen):
            system.order = self.calculate_order(ensemble)
            ekin = calculate_kinetic_energy(system.particles)[0]
            system.particles.ekin = ekin
            status, success, stop, _ = self.add_to_path(path, system.copy(),
                                                        left, right)
            if stop:
                logger.debug('Stopping propagate at step: %i ', i)
                break
            if reverse:
                system.particles.vel *= -1.0
                self.integration_step(system)
                system.particles.vel *= -1.0
            else:
                self.integration_step(system)
        logger.debug('Propagate done: "%s" (success: %s)', status, success)
        return success, status

    @staticmethod
    def modify_velocities(ensemble, vel_settings):
        """Modify the velocities of the current state.

        This method will modify the velocities of a time slice.
        And it is part of the integrator since it, conceptually,
        fits here:  we are acting on the system and modifying it.

        Parameters
        ----------
        ensemble : dict
            It contains:

            * `system : object like :py:class:`.System`
              This is the system that contains the particles we are
              investigating
            * `rgen` : object like :py:class:`.RandomGenerator`
              This is the random generator that will be used.

        vel_settings: dict.
            It contains all the info for the velocity:

            * `sigma_v` : numpy.array, optional
              These values can be used to set a standard deviation (one
              for each particle) for the generated velocities.
            * `aimless` : boolean, optional
              Determines if we should do aimless shooting or not.
            * `momentum` : boolean, optional
              If True, we reset the linear momentum to zero after
              generating.
            * `rescale or rescale_energy` : float, optional
              In some NVE simulations, we may wish to re-scale the
              energy to a fixed value. If `rescale` is a float > 0,
              we will re-scale the energy (after modification of
              the velocities) to match the given float.

        Returns
        -------
        dek : float
            The change in the kinetic energy.
        kin_new : float
            The new kinetic energy.

        """
        rgen = ensemble['rgen']
        system = ensemble['system']
        rescale = vel_settings.get('rescale_energy',
                                   vel_settings.get('rescale'))
        particles = system.particles
        if rescale is not None and rescale is not False:
            if rescale > 0:
                kin_old = rescale - particles.vpot
                do_rescale = True
            else:
                logger.warning('Ignored re-scale 6.2%f < 0.0.', rescale)
                return 0.0, calculate_kinetic_energy(particles)[0]
        else:
            kin_old = calculate_kinetic_energy(particles)[0]
            do_rescale = False
        if vel_settings.get('aimless', False):
            vel, _ = rgen.draw_maxwellian_velocities(system)
            particles.vel = vel
        else:  # Soft velocity change, from a Gaussian distribution:
            dvel, _ = rgen.draw_maxwellian_velocities(
                system, sigma_v=vel_settings['sigma_v'])
            particles.vel = particles.vel + dvel
        if vel_settings.get('momentum', False):
            reset_momentum(particles)
        if do_rescale:
            system.rescale_velocities(rescale)
        kin_new = calculate_kinetic_energy(particles)[0]
        dek = kin_new - kin_old
        return dek, kin_new

    def __call__(self, system):
        """To allow calling `MDEngine(system)`.

        Here, we are just calling `self.integration_step(system)`.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The system we are integrating.

        Returns
        -------
        out : None
            Does not return anything, but will update the particles.

        """
        return self.integration_step(system)

    @staticmethod
    def calculate_order(ensemble, xyz=None, vel=None, box=None):
        """Return the order parameter.

        This method is just to help to calculate the order parameter
        in cases where only the engine can do it.

        Parameters
        ----------
        ensemble : dict
            It contains:

            * `system` : object like :py:class:`.System`
              This is the system that contains the particles we are
              investigating
            * `order_function` : object like :py:class:`.OrderParameter`
              The object used for calculating the order parameter.

        xyz : numpy.array, optional
            The positions to use. Typically for internal engines, this
            is not needed. It is included here as it can be used for
            testing and also to be compatible with the generic function
            defined by the parent.
        vel : numpy.array, optional
            The velocities to use.
        box : numpy.array, optional
            The current box vectors.

        Returns
        -------
        out : list of floats
            The calculated order parameter(s).

        """
        system = ensemble['system']
        order_function = ensemble['order_function']
        if any((xyz is None, vel is None, box is None)):
            return order_function.calculate(system)
        system.particles.pos = xyz
        system.particles.vel = vel
        system.box.update_size(box)
        return order_function.calculate(system)

    def kick_across_middle(self, ensemble, middle, tis_settings):
        """Force a phase point across the middle interface.

        This is accomplished by repeatedly kicking the phase point so
        that it crosses the middle interface.

        Parameters
        ----------
        ensemble : dict
            It contains:

            * `system` : object like :py:class:`.System`
              This is the system that contains the particles we are
              investigating
            * `order_function` : object like :py:class:`.OrderParameter`
              The object used for calculating the order parameter.
            * `rgen` : object like :py:class:`.RandomGenerator`
              This is the random generator that will be used.

        middle : float
            This is the value for the middle interface.
        tis_settings : dict
            This dictionary contains settings for TIS. Explicitly used here:

            * `zero_momentum`: boolean, determines if the momentum is zeroed.
            * `rescale_energy`: boolean, determines if energy is re-scaled.

        Returns
        -------
        out[0] : object like :py:class:`.System`
            The phase-point just before the interface.
        out[1] : object like :py:class:`.System`
            The phase-point just after the interface.

        Note
        ----
        This function will update the system state.

        """
        # We search for crossing with the middle interface and do this
        # by sequentially kicking the initial phase point:
        system = ensemble['system']
        previous = system.copy()
        system.potential_and_force()  # Make sure forces are set.
        curr = self.calculate_order(ensemble)[0]
        logger.info('Kicking from: %9.6f', curr)
        while True:
            # Save current state:
            previous = system.copy()
            # Modify velocities:
            self.modify_velocities(
                ensemble, {'sigma_v': None,
                           'aimless': True,
                           'momentum': tis_settings['zero_momentum'],
                           'rescale': tis_settings['rescale_energy']})
            # Update order parameter in case it is velocity dependent:
            curr = self.calculate_order(ensemble)[0]
            previous.order = curr
            # Store modified velocities:
            previous.particles.set_vel(system.particles.get_vel())
            # Integrate forward one step:
            self.integration_step(system)
            # Compare previous order parameter and the new one:
            prev = curr
            curr = self.calculate_order(ensemble)[0]
            if curr == middle:
                # By construction we want two points, one left and one
                # right of the interface, and these two points should
                # be connected by a MD step. If we hit exactly on the
                # interface we just fall back:
                system.particles = previous.particles.copy()
                curr = previous.order
            else:
                if (prev <= middle < curr) or (curr < middle <= prev):
                    # Middle interface was crossed, just stop the loop.
                    logger.info('Crossing found: %9.6f %9.6f ', prev, curr)
                    break
                if (prev <= curr < middle) or (middle < curr <= prev):
                    # We are getting closer, keep the new point.
                    pass
                else:  # We did not get closer, fall back to previous point.
                    system.particles = previous.particles.copy()
                    curr = previous.order
        system.order = curr
        return previous, system

    def dump_phasepoint(self, phasepoint, deffnm=None):
        """For compatibility with external integrators."""
        return

    def clean_up(self):
        """Clean up after using the engine.

        Currently, this is only included for compatibility with external
        integrators.

        """
        return


class Verlet(MDEngine):
    """The Verlet MD integrator.

    This class defines the Verlet MD integrator.

    Attributes
    ----------
    half_idt : float
        Half of the inverse time step: `0.5 / timestep`.
    timestepsq : float
        Squared time step: `timestep**2`.
    previous_pos : numpy.array
        Stores the previous positions of the particles.

    """

    def __init__(self, timestep):
        """Set up the Verlet MD integrator.

        Parameters
        ----------
        timestep : float
            The time step in internal units.

        """
        super().__init__(timestep, 'Verlet MD integrator', dynamics='NVE')
        self.half_idt = 0.5 / self.timestep
        self.timestepsq = self.timestep**2
        self.previous_pos = None

    def set_initial_positions(self, particles):
        """Get initial positions for the Verlet integration.

        Parameters
        ----------
        particles : object like :py:class:`.Particles`
            The initial configuration. Positions and velocities are
            required.

        """
        self.previous_pos = particles.pos - particles.vel * self.timestep

    def integration_step(self, system):
        """Perform one Verlet integration step.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The system to integrate/act on. Assumed to have a particle
            list in ``system.particles``.

        Returns
        -------
        out : None
            Does not return anything, but alters the state of the given
            system.

        """
        particles = system.particles
        acc = particles.force * particles.imass
        pos = 2.0 * particles.pos - self.previous_pos + acc * self.timestepsq
        particles.vel = (pos - self.previous_pos) * self.half_idt
        self.previous_pos, particles.pos = particles.pos, pos
        system.potential_and_force()

    def restart_info(self):
        """Return restart info.

        The restart info for the engine Verlet.

        Returns
        -------
        info : dict,
            Contains all the updated simulation settings and counters.

        """
        info = super().restart_info()
        info['class'] = 'Verlet'
        for attr, value in self.__dict__.items():
            info[attr] = value

        return info

    def load_restart_info(self, info):
        """Load restart information.

        Parameters
        ----------
        info : dict
            The dictionary with the restart information, should be
            similar to the dict produced by :py:func:`.restart_info`.

        """
        super().load_restart_info(info)
        self.timestep = info.get('timestep')
        self.half_idt = info.get('half_idt')
        self.timestepsq = info.get('timestepsq')
        self.previous_pos = info.get('previous_pos')


class RandomWalk(MDEngine):
    """A Random Walker integrator.

    This class defines a Random walker integrator.

    Attributes
    ----------
    timestep : float
        The length of the step.

    """

    def __init__(self, timestep, rgen=None, seed=0):
        """Set up the Random walker integrator.

        Parameters
        ----------
        timestep : float
            The time step in internal units.
        rgen : string, optional
            This string can be used to pick a particular random
            generator, which is useful for testing.
        seed : integer, optional
            A seed for the random generator.

        """
        super().__init__(timestep, 'Random Walker integrator',
                         dynamics='NoNe')
        rgen_settings = {'seed': seed, 'rgen': rgen}
        self.rgen = create_random_generator(rgen_settings)

    def integration_step(self, system):
        """Random Walker integration, one time step.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The system to integrate/act on. Assumed to have a particle
            list in ``system.particles``.

        Returns
        -------
        out : None
            Does not return anything, but alters the state of the given
            `system`.

        """
        particles = system.particles
        particles.vel, _ = self.rgen.draw_maxwellian_velocities(system)
        particles.pos += self.timestep * particles.vel


class VelocityVerlet(MDEngine):
    """The Velocity Verlet MD integrator.

    This class defines the Velocity Verlet integrator.

    Attributes
    ----------
    half_timestep : float
        Half of the timestep.

    """

    def __init__(self, timestep):
        """Set up the Velocity Verlet integrator.

        Parameters
        ----------
        timestep : float
            The time step in internal units.

        """
        super().__init__(timestep, 'Velocity Verlet MD integrator',
                         dynamics='NVE')
        self.half_timestep = self.timestep * 0.5

    def integration_step(self, system):
        """Velocity Verlet integration, one time step.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The system to integrate/act on. Assumed to have a particle
            list in ``system.particles``.

        Returns
        -------
        out : None
            Does not return anything, but alters the state of the given
            `system`.

        """
        particles = system.particles
        imass = particles.imass
        particles.vel += self.half_timestep * particles.force * imass
        particles.pos += self.timestep * particles.vel
        system.potential_and_force()
        particles.vel += self.half_timestep * particles.force * imass

    def restart_info(self):
        """Return some restart info.

        The restart info for the engine VelocityVerlet.

        Returns
        -------
        info : dict,
            Contains all the updated simulation settings and counters.

        """
        info = super().restart_info()
        info['class'] = 'VelocityVerlet'
        for attr, value in self.__dict__.items():
            info[attr] = value

        return info

    def load_restart_info(self, info):
        """Load restart information.

        Parameters
        ----------
        info : dict
            The dictionary with the restart information, should be
            similar to the dict produced by :py:func:`.restart_info`.

        """
        super().load_restart_info(info)
        self.timestep = info.get('timestep')
        self.half_timestep = info.get('half_timestep')
        self.dynamics = info.get('dynamics')


class Langevin(MDEngine):
    """The Langevin MD integrator.

    This class defines a Langevin integrator.

    Attributes
    ----------
    rgen : object like :py:class:`.RandomGenerator`
        This is the class that handles the generation of random numbers.
    gamma : float
        The friction parameter.
    high_friction : boolean
        Determines if we are in the high friction limit and should
        do the over-damped version.
    init_params : boolean
        If true, we will initiate parameters for the Langevin
        integrator when `integrate_step` is invoked.
    param_high : dict
        This contains the parameters for the high friction limit. Here
        we integrate the equations of motion according to:
        ``r(t + dt) = r(t) + dt * f(t)/m*gamma + dr``.
        The items in the dict are:

        * `sigma` : float
          standard deviation for the positions, used when drawing dr.
        * `bddt` : numpy.array
          Equal to ``dt*gamma/masses``, since the masses is a
          numpy.array this will have the same shape.

    param_iner : dict
        This dict contains the parameters for the non-high friction
        limit where we integrate the equations of motion according to:
        ``r(t + dt) = r(t) + c1 * dt * v(t) + c2*dt*dt*a(t) + dr``
        and
        ``v(r + dt) = c0 * v(t) + (c1-c2)*dt*a(t) + c2*dt*a(t+dt) + dv``.
        The dict contains:

        * `c0` : float
          Corresponds to ``c0`` in the equation above.
        * `a1` : float
          Corresponds to ``c1*dt`` in the equation above.
        * `a2` : numpy.array
          Corresponds to ``c2*dt*dt/mass`` in the equation above.
          Here we divide by the masses in order to use the forces rather
          than the acceleration. Since the masses might be different for
          different particles, this will result in a numpy.array with
          shape equal to the shape of the masses.
        * `b1` : numpy.array
          Corresponds to ``(c1-c2)*dt/mass`` in the equation above.
          Here we also divide by the masses, resulting in a numpy.array.
        * `b2` : numpy.array
          Corresponds to ``c2*dt/mass`` in the equation above.
          Here we also divide by the masses, resulting in a numpy.array.
        * `mean` : numpy.array (2,)
          The means for the bivariate Gaussian distribution.
        * `cov` : numpy.array (2,2)
          This array contains the covariance for the bivariate Gaussian
          distribution. `param_iner['mean']` and `param_iner['cov']` are
          used as parameters when drawing ``dr`` and ``dv`` from the
          bivariate distribution.

    Note
    ----
    Currently, we are using a multi-normal distribution from numpy.
    Consider replacing this one as it seems somewhat slow.

    """

    def __init__(self, timestep, gamma, rgen=None, seed=0,
                 high_friction=False):
        """Set up the Langevin integrator.

        Actually, it is very convenient to set some variables for the
        different particles. However, to have a uniform initialisation
        for the different integrators, we postpone this.
        This initialisation can be done later by calling explicitly the
        function `self._init_parameters(system)` or it will be called
        the first time `self.integration_step` is invoked.

        Parameters
        ----------
        timestep : float
            The time step in internal units.
        gamma : float
            The gamma parameter for the Langevin integrator.
        rgen : string, optional
            This string can be used to pick a particular random
            generator, which is useful for testing.
        seed : integer, optional
            A seed for the random generator.
        high_friction : boolean, optional
            Determines if we are in the high_friction limit and should
            do the over-damped version.

        """
        super().__init__(timestep, 'Langevin MD integrator',
                         dynamics='stochastic')
        self.gamma = gamma
        self.high_friction = high_friction
        rgen_settings = {'seed': seed, 'rgen': rgen}
        self.rgen = create_random_generator(rgen_settings)
        self.param_high = {'sigma': None, 'bddt': None}
        self.param_iner = {'c0': None, 'a1': None, 'a2': None,
                           'b1': None, 'b2': None, 'mean': None, 'cov': None}
        self.init_params = True

    def _init_parameters(self, system):
        """Extra initialisation of the Langevin integrator.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The system to integrate/act on. Assumed to have a particle
            list in ``system.particles``.

        Returns
        -------
        out : None
            Does not return anything, but updates ``self.param``.

        """
        beta = system.temperature['beta']
        imasses = system.particles.imass
        if self.high_friction:
            self.param_high['sigma'] = np.sqrt(2.0 * self.timestep *
                                               imasses/(beta * self.gamma))
            self.param_high['bddt'] = self.timestep * imasses / self.gamma
        else:
            gammadt = self.gamma * self.timestep
            exp_gdt = np.exp(-gammadt)
            if self.gamma > 0.0:
                c_0 = exp_gdt
                c_1 = (1.0 - c_0) / gammadt
                c_2 = (1.0 - c_1) / gammadt
            else:
                raise ValueError(
                    'Langevin: Found gamma = {:6.2f} < 0.'.format(self.gamma)
                )

            self.param_iner['c0'] = c_0
            self.param_iner['a1'] = c_1 * self.timestep
            self.param_iner['a2'] = c_2 * self.timestep**2 * imasses
            self.param_iner['b1'] = (c_1 - c_2) * self.timestep * imasses
            self.param_iner['b2'] = c_2 * self.timestep * imasses

            self.param_iner['mean'] = []
            self.param_iner['cov'] = []
            self.param_iner['cho'] = []

            for imass in imasses:
                sig_ri2 = ((self.timestep * imass[0] / (beta * self.gamma)) *
                           (2. - (3. - 4.*exp_gdt + exp_gdt**2) / gammadt))
                sig_vi2 = (1.0 - exp_gdt**2) * imass[0] / beta
                cov_rvi = (imass[0]/(beta * self.gamma)) * (1.0 - exp_gdt)**2
                cov_matrix = np.array([[sig_ri2, cov_rvi],
                                       [cov_rvi, sig_vi2]])
                self.param_iner['cov'].append(cov_matrix)
                self.param_iner['cho'].append(np.linalg.cholesky(cov_matrix))
                self.param_iner['mean'].append(np.zeros(2))

    def integration_step(self, system):
        """Langevin integration, one time step.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The system to integrate/act on. Assumed to have a particle
            list in ``system.particles``.

        Returns
        -------
        out : None
            Does not return anything, but alters the state of the given
            system.

        """
        if self.init_params:
            self._init_parameters(system)
            self.init_params = False
        if self.high_friction:
            return self.integration_step_overdamped(system)
        return self.integration_step_inertia(system)

    def integration_step_overdamped(self, system):
        """Over damped Langevin integration, one time step.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The system to integrate/act on. Assumed to have a particle
            list in ``system.particles``.

        Returns
        -------
        out : None
            Does not return anything, but alters the state of the given
            `system`.

        """
        system.potential_and_force()  # update forces
        particles = system.particles
        rands = self.rgen.normal(loc=0.0, scale=self.param_high['sigma'],
                                 size=particles.vel.shape)
        particles.pos += self.param_high['bddt'] * particles.force + rands
        particles.vel = rands

    def integration_step_inertia(self, system):
        """Langevin integration, one time step.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The system to integrate/act on. Assumed to have a particle
            list in ``system.particles``.

        Returns
        -------
        out : None
            Does not return anything, but alters the state of the given
            `system`.

        """
        particles = system.particles
        pos_rand = np.zeros(particles.pos.shape)
        vel_rand = np.zeros(particles.vel.shape)
        if self.gamma > 0.0:
            mean, cov = self.param_iner['mean'], self.param_iner['cov']
            cho = self.param_iner['cho']
            for i, (meani, covi, choi) in enumerate(zip(mean, cov, cho)):
                randxv = self.rgen.multivariate_normal(meani, covi, cho=choi,
                                                       size=system.get_dim())
                pos_rand[i] = randxv[:, 0]
                vel_rand[i] = randxv[:, 1]
        particles.pos += (self.param_iner['a1'] * particles.vel +
                          self.param_iner['a2'] * particles.force + pos_rand)

        vel2 = (self.param_iner['c0'] * particles.vel +
                self.param_iner['b1'] * particles.force + vel_rand)

        system.potential_and_force()  # Update forces.

        particles.vel = vel2 + self.param_iner['b2'] * particles.force

    def restart_info(self):
        """Return restart info.

        The restart info for the engine Langevin.

        Returns
        -------
        info : dict
            Contains all the updated simulation settings and counters.

        """
        info = super().restart_info()
        info['class'] = 'Langevin'
        for attr, value in self.__dict__.items():
            if attr == 'rgen':
                info[attr] = value.get_state()
            else:
                info[attr] = value
        return info

    def load_restart_info(self, info):
        """Load restart information.

        Parameters
        ----------
        info : dict
            The dictionary with the restart information, should be
            similar to the dict produced by :py:func:`.restart_info`.

        """
        super().load_restart_info(info)
        self.timestep = info.get('timestep')
        self.dynamics = info.get('dynamics')
        self.gamma = info.get('gamma')
        self.high_friction = info.get('high_friction')
        self.param_high = info.get('param_high')
        self.param_iner = info.get('param_iner')
        self.init_params = info.get('init_params')

        if 'rgen' in info:
            self.rgen = create_random_generator(info['rgen'])

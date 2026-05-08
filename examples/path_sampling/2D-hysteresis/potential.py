# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This is a 2D example potential."""
import logging
import numpy as np
from pyretis.forcefield.potential import PotentialFunction
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


class Hyst2D(PotentialFunction):
    r"""Hyst2D(PotentialFunction).

    This class defines a 2D dimensional potential with two
    stable states.
    The potential energy (:math:`V_\text{pot}`) is given by

    .. math::

       V_\text{pot}(x, y) = \gamma_1 (x^2 + y^2)^2 +
       \gamma_2 \exp(\alpha_1 (x - x_0)^2 + \alpha_2 (y - y_0)^2) +
       \gamma_3 \exp(\beta_1 (x + x_0)^2 + \beta_2(y + y_0)^2)

    where :math:`x` and :math:`y` gives the positions and :math:`\gamma_1`,
    :math:`\gamma_2`, :math:`\gamma_3`, :math:`\alpha_1`, :math:`\alpha_2`,
    :math:`\beta_1`, :math:`\beta_2`, :math:`x_0` and :math:`y_0` are
    potential parameters.

    Attributes
    ----------
    params : dict
        Contains the parameters. The keys are:

        * `gamma1`: The :math:`\gamma_1` parameter for the potential.
        * `gamma2`: The :math:`\gamma_2` parameter for the potential.
        * `gamma3`: The :math:`\gamma_3` parameter for the potential.
        * `alpha1`: The :math:`\alpha_1` parameter for the potential.
        * `alpha2`: The :math:`\alpha_2` parameter for the potential.
        * `beta1`: The :math:`\beta_1` parameter for the potential.
        * `beta2`: The :math:`\beta_2` parameter for the potential.
        * `x0`: The :math:`x_0` parameter for the potential.
        * `y0`: The :math:`y_0` parameter for the potential.

    """

    def __init__(self, desc='2D hysteresis'):
        """Set up the potential.

        Parameters
        ----------
        a : float, optional
            Parameter for the potential.
        b : float, optional
            Parameter for the potential.
        c : float, optional
            Parameter for the potential.
        desc : string, optional
            Description of the force field.

        """
        super().__init__(dim=2, desc=desc)
        self.params = {'gamma1': 0.0, 'gamma2': 0.0, 'gamma3': 0.0,
                       'alpha1': 0.0, 'alpha2': 0.0, 'beta1': 0.0,
                       'beta2': 0.0, 'x0': 0.0, 'y0': 0.0}

    def potential(self, system):
        """Evaluate the potential.

        Parameters
        ----------
        system : object like `System`
            The system we evaluate the potential for. Here, we
            make use of the positions only.

        Returns
        -------
        out : float
            The potential energy.

        """
        x = system.particles.pos[:, 0]  # pylint: disable=invalid-name
        y = system.particles.pos[:, 1]  # pylint: disable=invalid-name
        gam1 = self.params['gamma1']
        gam2 = self.params['gamma2']
        gam3 = self.params['gamma3']
        alf1 = self.params['alpha1']
        alf2 = self.params['alpha2']
        bet1 = self.params['beta1']
        bet2 = self.params['beta2']
        x_0 = self.params['x0']
        y_0 = self.params['y0']
        v_pot = (gam1 * (x**2 + y**2)**2 +
                 gam2 * np.exp(alf1 * (x - x_0)**2 + alf2 * (y - y_0)**2) +
                 gam3 * np.exp(bet1 * (x + x_0)**2 + bet2 * (y + y_0)**2))
        return v_pot.sum()

    def force(self, system):
        """Evaluate forces.

        Parameters
        ----------
        system : object like `System`
            The system we evaluate the potential for. Here, we
            make use of the positions only.

        Returns
        -------
        out[0] : numpy.array
            The calculated force.
        out[1] : numpy.array
            The virial, currently not implemented for this potential.

        """
        x = system.particles.pos[:, 0]  # pylint: disable=invalid-name
        y = system.particles.pos[:, 1]  # pylint: disable=invalid-name
        gam1 = self.params['gamma1']
        gam2 = self.params['gamma2']
        gam3 = self.params['gamma3']
        alf1 = self.params['alpha1']
        alf2 = self.params['alpha2']
        bet1 = self.params['beta1']
        bet2 = self.params['beta2']
        x_0 = self.params['x0']
        y_0 = self.params['y0']
        term = 4.0 * gam1 * (x**2 + y**2)
        exp1 = gam2 * np.exp(alf1 * (x - x_0)**2 + alf2 * (y - y_0)**2)
        exp2 = gam3 * np.exp(bet1 * (x + x_0)**2 + bet2 * (y + y_0)**2)
        forces = np.zeros_like(system.particles.pos)
        forces[:, 0] = -(x * term +
                         2.0 * alf1 * (x - x_0) * exp1 +
                         2.0 * bet1 * (x + x_0) * exp2)
        forces[:, 1] = -(y * term +
                         2.0 * alf2 * (y - y_0) * exp1 +
                         2.0 * bet2 * (y + y_0) * exp2)
        virial = np.zeros((self.dim, self.dim))  # just return zeros here
        return forces, virial

    def potential_and_force(self, system):
        """Evaluate the potential and the force.

        Parameters
        ----------
        system : object like `System`
            The system we evaluate the potential for. Here, we
            make use of the positions only.

        Returns
        -------
        out[0] : float
            The potential energy as a float.
        out[1] : numpy.array
            The force as a numpy.array of the same shape as the
            positions in `particles.pos`.
        out[2] : numpy.array
            The virial, currently not implemented for this potential.

        """
        virial = np.zeros((self.dim, self.dim))  # just return zeros here
        x = system.particles.pos[:, 0]  # pylint: disable=invalid-name
        y = system.particles.pos[:, 1]  # pylint: disable=invalid-name
        gam1 = self.params['gamma1']
        gam2 = self.params['gamma2']
        gam3 = self.params['gamma3']
        alf1 = self.params['alpha1']
        alf2 = self.params['alpha2']
        bet1 = self.params['beta1']
        bet2 = self.params['beta2']
        x_0 = self.params['x0']
        y_0 = self.params['y0']
        term0 = (x**2 + y**2)
        exp1 = gam2 * np.exp(alf1 * (x - x_0)**2 + alf2 * (y - y_0)**2)
        exp2 = gam3 * np.exp(bet1 * (x + x_0)**2 + bet2 * (y + y_0)**2)
        v_pot = gam1 * term0**2 + exp1 + exp2
        term = 4.0 * gam1 * term0
        forces = np.zeros_like(system.particles.pos)
        forces[:, 0] = -(x * term +
                         2.0 * alf1 * (x - x_0) * exp1 +
                         2.0 * bet1 * (x + x_0) * exp2)
        forces[:, 1] = -(y * term +
                         2.0 * alf2 * (y - y_0) * exp1 +
                         2.0 * bet2 * (y + y_0) * exp2)
        return v_pot.sum(), forces, virial

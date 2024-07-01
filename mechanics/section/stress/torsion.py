from typing import Type
from math import pi
import sympy as sp
from mechanics import Quantity
from mechanics.section.geometry import Shape, Circle, Annulus

Q_ = Quantity


class Torsion:
    """
    Determination of the maximum shear stress in a circular or annular
    cross-section due to torsional loading.
    """
    def __init__(
        self,
        T: Quantity,
        shape: Shape,
    ) -> None:
        """Creates a `Torsion` object.

        Parameters
        ----------
        T:
            Resultant internal torque in a circular or annular cross-section.
        shape:
            Shape of the cross-section (either a `Circle` or `Annulus` object).
        """
        self.T = T
        self.shape = shape
        self.tau_max = self.__calc_max_shear_stress()

    @classmethod
    def design(
        cls,
        shape: Type[Shape],
        T: Quantity,
        t: Quantity | None = None,
        tau_allow: Quantity | None = None,
        tau_fail: Quantity | None = None,
        safety_factor: float | None = None
    ) -> Quantity:
        """Returns the required outer radius of a round shaft for a given
        resultant internal torque.

        Parameters
        ----------
        shape:
            Type of shape of the shaft's cross-section. Only class `Circle` (for
            a solid shaft) or class `Annulus` (for a hollow shaft) are
            recognized.
        T:
            The resultant internal torque in the cross-section.
        t:
            Wall thickness in case of a hollow shaft.
        tau_allow:
            Allowable shear stress.
        tau_fail:
            Failure shear stress.
        safety_factor:
            Factor of safety to be applied to the failure stress.

        Notes
        -----
        Either parameter `tau_allow` alone, or parameter `tau_fail` together
        with parameter `safety_factor` must be specified.
        """
        if (tau_fail is not None) and (safety_factor is not None):
            tau_allow = cls.__calc_tau_allow(tau_fail, safety_factor)
        if issubclass(shape, Circle):
            r = (2 * T / (tau_allow * pi)) ** (1 / 3)
            return r
        if issubclass(shape, Annulus):
            r_o = cls.__calc_outer_radius(T, tau_allow, t)
            return r_o

    @classmethod
    def allowable_torque(
        cls,
        shape: Shape,
        tau_allow: Quantity | None = None,
        tau_fail: Quantity | None = None,
        safety_factor: float | None = None
    ) -> Quantity:
        """Returns the allowable torque that can be applied to a round shaft with
        known dimensions.

        Parameters
        ----------
        shape:
            Shape of the shaft's cross-section. Only instances of class `Circle`
            (in case of a solid shaft) or class `Annulus` (in case of a hollow
            shaft or tube) are recognized.
        tau_allow:
            Allowable shear stress.
        tau_fail:
            Failure shear stress.
        safety_factor:
            Factor of safety.

        Notes
        -----
        Either parameter `tau_allow` alone, or parameter `tau_fail` together
        with parameter `safety_factor` must be specified.
        """
        if (tau_fail is not None) and (safety_factor is not None):
            tau_allow = cls.__calc_tau_allow(tau_fail, safety_factor)
        if isinstance(shape, (Circle, Annulus)):
            if isinstance(shape, Circle):
                r = shape.dim.radius
            else:
                r = shape.dim.outer_radius
            J = shape.polar_moment_of_inertia
            T_allow = tau_allow * J / r
            return T_allow

    def __calc_max_shear_stress(self) -> Quantity:
        """Returns the shear stress at the outer radius of a round shaft or
        tube.
        """
        if isinstance(self.shape, (Circle, Annulus)):
            if isinstance(self.shape, Circle):
                r = self.shape.dim.radius
            else:
                r = self.shape.dim.outer_radius
            J = self.shape.polar_moment_of_inertia
            tau_max = self.T * r / J
            return tau_max

    @staticmethod
    def __calc_tau_allow(
        tau_fail: Quantity,
        safety_factor: float
    ) -> Quantity:
        """Calculates the allowable shear stress from the given failure stress
        and safety factor.
        """
        tau_allow = tau_fail / safety_factor
        return tau_allow

    @staticmethod
    def __calc_outer_radius(
        T: Quantity,
        tau_allow: Quantity,
        t: Quantity
    ) -> Quantity:
        """Calculates the required outer radius of a round hollow shaft or
        tube.
        """
        T = T.to('N * m').m
        tau_allow = tau_allow.to('Pa').m
        t = t.to('m').m
        r_o = sp.Symbol('r_o')
        r_i = r_o - t
        k = tau_allow * pi / (2 * T)
        lhs = k * (r_o ** 4 - r_i ** 4)
        eq = sp.Eq(lhs, r_o)
        solutions = sp.solve(eq)
        solutions = [s for s in solutions if float(s) > t]
        r_o = Q_(solutions[0], 'm')
        return r_o

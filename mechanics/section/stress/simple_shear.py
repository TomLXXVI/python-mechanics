from mechanics import Quantity
from mechanics.section.geometry import Shape


class SimpleShear:
    """
    Determination of the average shear stress in a cross-section due to shear
    loading.

    Notes
    -----
    This class can be used for small elements like bolts, pins, welds, etc.
    """
    def __init__(
        self,
        V: Quantity,
        shape: Shape | None = None,
        A: Quantity | None = None
    ) -> None:
        """Creates a `SimpleShear` instance.

        Parameters
        ----------
        V:
            Resultant internal shear force in a cross-section.
        shape:
            Shape of the cross-section.
        A:
            Area of the cross-section.

        Notes
        -----
        Either `shape` or `A` must be specified. `A` can be used if the shape is
        not known, but only the area of the cross-section is given.
        """
        self.V = V
        if shape is not None:
            self.shape = shape
            self.A = shape.area
        else:
            self.shape = None
            self.A = A
        self.tau = self.__calc_avg_shear()

    @classmethod
    def design(
        cls,
        V: Quantity,
        tau_allow: Quantity | None = None,
        tau_fail: Quantity | None = None,
        safety_factor: float | None = None
    ) -> Quantity:
        """Returns the required area of the cross-section.

        Parameters
        ----------
        V:
            Shear force acting at the cross-section.
        tau_allow:
            The allowable shear stress.
        tau_fail:
            The failure shear stress of the material.
        safety_factor:
            Factor of safety to be applied to the failure stress.

        Notes
        -----
        Either parameter `tau_allow` alone, or parameter `tau_fail` together
        with parameter `safety_factor` must be specified.
        """
        if (tau_fail is not None) and (safety_factor is not None):
            tau_allow = cls.__calc_tau_allow(tau_fail, safety_factor)
        A = V / tau_allow
        return A

    def __calc_avg_shear(self) -> Quantity:
        """Returns the average shear stress in the cross-section."""
        tau = self.V / self.A
        return tau

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

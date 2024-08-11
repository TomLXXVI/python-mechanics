from mechanics import Quantity
from mechanics.geometry.shapes import Shape


class AxialLoading:
    """
    Determine the (average) normal stress at a cross-section due to axial
    loading.
    """
    def __init__(
        self,
        N: Quantity,
        shape: Shape | None = None,
        A: Quantity | None = None
    ) -> None:
        """Creates an `AxialLoading` instance.

        Parameters
        ----------
        N:
            Resultant internal normal force in a cross-section.
        shape:
            Shape of the cross-section.
        A:
            Area of the cross-section.

        Notes
        -----
        Either `shape` or `A` must be specified. `A` can be used if the shape is
        not known, but only the area of the cross-section is given.
        """
        self.N = N
        if shape is not None:
            self.shape = shape
            self.A = self.shape.area
        else:
            self.shape = None
            self.A = A
        self.sigma = self.__calc_avg_stress()

    @classmethod
    def design(
        cls,
        N: Quantity,
        sigma_allow: Quantity | None = None,
        sigma_fail: Quantity | None = None,
        safety_factor: float | None = None
    ) -> Quantity:
        """Returns the required area of the cross-section.

        Parameters
        ----------
        N:
            Normal force acting at the cross-section.
        sigma_allow:
            The allowable normal stress.
        sigma_fail:
            The failure normal stress of the material.
        safety_factor:
            Factor of safety to be applied to the failure stress.

        Notes
        -----
        Either parameter `sigma_allow` alone, or parameter `sigma_fail` together
        with parameter `safety_factor` must be specified.
        """
        if (sigma_fail is not None) and (safety_factor is not None):
            sigma_allow = cls.__calc_sigma_allow(sigma_fail, safety_factor)
        A = N / sigma_allow
        return A

    def __calc_avg_stress(self) -> Quantity:
        """Returns the average normal stress in the cross-section."""
        sigma = self.N / self.A
        return sigma

    @staticmethod
    def __calc_sigma_allow(
        sigma_fail: Quantity,
        safety_factor: float
    ) -> Quantity:
        """Calculates the allowable normal stress from the given failure stress
        and safety factor.
        """
        sigma_allow = sigma_fail / safety_factor
        return sigma_allow

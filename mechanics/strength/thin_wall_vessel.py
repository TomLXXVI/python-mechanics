# Thin-walled pressure vessels
# The ratio of inner radius to wall thickness must be 10 or more.

from mechanics import Quantity


class CylindricalVessel:

    def __init__(
        self,
        inner_radius: Quantity,
        thickness: Quantity,
        pressure: Quantity | None = None,
        sigma_allow: Quantity | None = None
    ) -> None:
        if (ratio := inner_radius.to('mm') / thickness.to('mm')) < 10:
            raise ValueError(
                "the vessel is not thin-walled: "
                f"inner-radius-to-wall-thickness ratio = {ratio} < 10.0"
            )
        self.inner_radius = self.r_i = inner_radius
        self.thickness = self.t = thickness
        self.pressure = self.P = pressure
        self.outer_radius = self.r_o = inner_radius + thickness
        self.sigma_allow = sigma_allow

    @property
    def circumferential_stress(self) -> Quantity:
        if self.P is not None:
            sigma1 = self.P * self.r_i / self.t
            return sigma1.to('MPa')

    @property
    def longitudinal_stress(self) -> Quantity:
        sigma2 = self.circumferential_stress / 2
        return sigma2

    @property
    def maximum_pressure(self) -> Quantity:
        if self.sigma_allow is not None:
            P_max = self.sigma_allow * self.t / self.r_i
            return P_max.to('Pa')


class SphericalVessel(CylindricalVessel):

    @property
    def circumferential_stress(self) -> Quantity:
        sigma1 = self.longitudinal_stress
        return sigma1

    @property
    def maximum_pressure(self) -> Quantity:
        if self.sigma_allow is not None:
            P_max = 2 * self.sigma_allow * self.t / self.r_i
            return P_max.to('Pa')

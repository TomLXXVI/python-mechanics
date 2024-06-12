from typing import Type
from mechanics import Quantity
from mechanics.geometry import create_shape, Dimensions, Shape
from .normal_stress import NormalStress
from .simple_shear import ShearStress


class Section:

    def __init__(
        self,
        shape: Type[Shape] | None = None,
        dim: Dimensions | None = None
    ) -> None:
        if isinstance(dim, Dimensions):
            self.shape = create_shape(shape, dim)
        elif shape is not None and issubclass(shape, Shape):
            self.shape = shape
        else:
            self.shape = None
        self.N: Quantity | None = None
        self.V: Quantity | None = None
        self.M: Quantity | None = None
        self.T: Quantity | None = None

    def set_resultant_internal_loadings(self, d: dict[str, Quantity]) -> None:
        self.N = d.get('F_x')
        self.V = d.get('F_y')
        self.M = d.get('M_z')
        self.T = d.get('M_x')

    def set_normal_force(self, N: Quantity) -> None:
        self.N = N

    def set_shear_force(self, V: Quantity) -> None:
        self.V = V

    def set_bending_moment(self, M: Quantity) -> None:
        self.M = M

    def set_torque(self, T: Quantity) -> None:
        self.T = T

    @property
    def sigma_avg(self) -> Quantity:
        if isinstance(self.shape, Shape):
            ns = NormalStress(self.N, A=self.shape.area())
            return ns.sigma_avg

    @property
    def tau_avg(self) -> Quantity:
        if isinstance(self.shape, Shape):
            ss = ShearStress(self.V, A=self.shape.area())
            return ss.tau_avg

    def design(
        self,
        sigma_allow: Quantity | None = None,
        sigma_fail: Quantity | None = None,
        tau_allow: Quantity | None = None,
        tau_fail: Quantity | None = None,
        safety_factor: float | None = None
    ) -> tuple[Quantity | None, Quantity | None]:
        A_sigma, A_tau = None, None
        if self.N is not None:
            ns = NormalStress(self.N, sigma_allow, sigma_fail, safety_factor)
            A_sigma = ns.design()
        if self.V is not None:
            ss = ShearStress(self.V, tau_allow, tau_fail, safety_factor)
            A_tau = ss.design()
        return A_sigma, A_tau

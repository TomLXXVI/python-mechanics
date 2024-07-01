from typing import Type
from mechanics import Quantity
from .geometry import create_geometry, Dimensions, Shape
from .stress import AxialLoading, TransverseShear, Torsion, Bending


class Section:

    def __init__(
        self,
        shape: Type[Shape],
        dim: Dimensions
    ) -> None:
        self.shape = create_geometry(shape, dim)
        self.N: Quantity | None = None
        self.V: Quantity | None = None
        self.M_z: Quantity | None = None
        self.T: Quantity | None = None

    def set_resultant_internal_loadings(self, d: dict[str, Quantity]) -> None:
        self.N = d.get('F_x')
        self.V = d.get('F_y')
        self.M_z = d.get('M_z')
        self.T = d.get('M_x')

    def set_normal_force(self, N: Quantity) -> None:
        self.N = N

    def set_shear_force(self, V: Quantity) -> None:
        self.V = V

    def set_bending_moment(self, M: Quantity) -> None:
        self.M_z = M

    def set_torque(self, T: Quantity) -> None:
        self.T = T

    @property
    def sigma_N(self) -> Quantity:
        if self.N is not None:
            return AxialLoading(self.N, self.shape).sigma
        else:
            raise ValueError('no normal force is acting at the section')

    @property
    def tau_V(self) -> Quantity:
        if self.V is not None:
            return TransverseShear(self.V, self.shape).tau
        else:
            raise ValueError('no shear force is acting at the section')

    @property
    def tau_T(self) -> Quantity:
        if self.T is not None:
            return Torsion(self.T, self.shape).tau_max
        else:
            raise ValueError('no torque is acting at the section')

    @property
    def sigma_M(self) -> Quantity:
        if self.M_z is not None:
            return Bending(self.M_z, self.shape).sigma_max
        else:
            raise ValueError('no bending moment is acting at the section')

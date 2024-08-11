from typing import Type
from mechanics import Quantity
from mechanics.statics import Force, Moment
from mechanics.geometry.shapes import *
from .stress import *


Q_ = Quantity


class Section:
    """Represents the cross-section of a beam."""

    def __init__(
        self,
        shape: Type[Shape],
        dim: Dimensions
    ) -> None:
        """Creates a `Section` object.

        Parameters
        ----------
        shape:
            Shape of the cross-section (any class derived from abstract base
            class `Shape`, see section.geometry.shapes.py).
        dim:
            Dimensions of the cross-section (instance of dataclass `Dimensions`,
            see section.geometry.shapes.py).
        """
        self.shape = self.__create_geometry(shape, dim)
        self.N: Quantity | None = None
        self.V: Quantity | None = None
        self.M_z: Quantity | None = None
        self.M_y: Quantity | None = None
        self.T: Quantity | None = None
        self.axial: AxialLoading | None = None
        self.simple_shear: SimpleShear | None = None
        self.torsion: Torsion | None = None
        self.transverse_shear: TransverseShear | None = None
        self.bending: Bending | None = None

    def set_internal_loadings(self, F_i: Force, M_i: Moment) -> None:
        """Sets the resultant internal force and moment at the cross-section."""
        self.set_normal_force(F_i.x)
        self.set_shear_force(F_i.y)
        self.set_bending_moment(M_i.z, M_i.y)
        self.set_torque(M_i.x)

    def set_normal_force(self, N: Quantity) -> None:
        """Sets the resultant internal normal force at the cross-section."""
        if N.magnitude != 0.0:
            self.N = N
            self.axial = AxialLoading(self.N, self.shape)

    def set_shear_force(self, V: Quantity) -> None:
        """Sets the resultant internal shear force at the cross-section."""
        if V.magnitude != 0.0:
            self.V = V
            self.simple_shear = SimpleShear(self.V, self.shape)
            self.transverse_shear = TransverseShear(self.V, self.shape)

    def set_bending_moment(self, M_z: Quantity, M_y: Quantity | None = None) -> None:
        """Sets the resultant internal bending moment at the cross-section."""
        if M_z.magnitude != 0.0:
            self.M_z = M_z
        if isinstance(M_y, Quantity) and M_y.magnitude != 0.0:
            self.M_y = M_y
        if self.M_z is not None and self.M_y is not None:
            self.bending = Bending((self.M_z, self.M_y), self.shape)
        elif self.M_z is not None:
            self.bending = Bending(self.M_z, self.shape)
        elif self.M_y is not None:
            M_z = Q_(0.0, self.M_y.units)
            self.bending = Bending((M_z, self.M_y), self.shape)

    def set_torque(self, T: Quantity) -> None:
        """Sets the resultant internal torsion moment at the cross-section."""
        if T.magnitude != 0.0:
            self.T = T
            self.torsion = Torsion(self.T, self.shape)

    @property
    def sigma_N_avg(self) -> Quantity:
        """Returns the average normal stress in the cross-section due to the
        resultant internal normal force N at the section.
        """
        if self.N is not None:
            return self.axial.sigma
        else:
            raise ValueError('no normal force is acting at the section')

    @property
    def tau_V_avg(self) -> Quantity:
        """Returns the average shear stress in the cross-section due to the
        resultant internal shear force V at the section.
        Can be used for small elements like bolts, pins, welds, etc.
        """
        if self.V is not None:
            return self.simple_shear.tau
        else:
            raise ValueError('no shear force is acting at the section')

    @property
    def tau_V_max(self) -> Quantity:
        """Returns the maximum shear stress in the cross-section due to the
        resultant internal shear force V at the section.
        """
        if self.V is not None:
            return self.transverse_shear.tau_max
        else:
            raise ValueError('no shear force is acting at the section')

    @property
    def tau_T_max(self) -> Quantity:
        """Returns the maximum shear stress in the cross-section due to the
        resultant internal torsion moment T at the section.
        """
        if self.T is not None:
            return self.torsion.tau_max
        else:
            raise ValueError('no torque is acting at the section')

    @property
    def sigma_M_max(self) -> Quantity:
        """Returns the maximum normal stress in the cross-section due to the
        resultant internal bending moment M at the section.
        """
        if self.M_z is not None:
            return self.bending.sigma_max
        else:
            raise ValueError('no bending moment is acting at the section')

    @staticmethod
    def __create_geometry(shape: Type[Shape], dim: Dimensions) -> Shape:
        if shape is Rectangle:
            shape_obj = Rectangle(
                width=dim.width,
                height=dim.height
            )
            return shape_obj
        elif shape is Circle:
            shape_obj = Circle(
                radius=dim.radius or dim.diameter / 2
            )
            return shape_obj
        elif shape is Annulus:
            shape_obj = Annulus(
                outer_radius=dim.outer_radius,
                thickness=dim.thickness
            )
            return shape_obj
        elif shape is HShape:
            shape_obj = HShape(
                height=dim.height,
                width=dim.width,
                web_thickness=dim.web_thickness,
                flange_thickness=dim.flange_thickness
            )
            return shape_obj
        elif shape is CShape:
            shape_obj = CShape(
                height=dim.height,
                width=dim.width,
                web_thickness=dim.web_thickness,
                flange_thickness=dim.flange_thickness,
                orientation=dim.orientation
            )
            return shape_obj
        elif shape is ZShape:
            shape_obj = ZShape(
                height=dim.height,
                width=dim.width,
                web_thickness=dim.web_thickness,
                flange_thickness=dim.flange_thickness,
                orientation=dim.orientation
            )
            return shape_obj
        elif shape is HollowRectangle:
            shape_obj = HollowRectangle(
                width=dim.width,
                height=dim.height,
                thickness=dim.thickness
            )
            return shape_obj

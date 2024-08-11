from typing import Type
from abc import ABC
from dataclasses import dataclass
import numpy as np
from mechanics import Quantity
from mechanics.statics import Force, Moment, DistributedLoad1D
from mechanics import statics
from mechanics.geometry.shapes import Shape, Dimensions
from mechanics.strength.strain import AxialStrain, TorsionalStrain, ElasticCurve, BoundaryCondition
from .section import Section

Q_ = Quantity


@dataclass
class Support(ABC):
    x: Quantity

    def __post_init__(self):
        self.bc = BoundaryCondition(self.x)


@dataclass
class Roller(Support):

    def __post_init__(self):
        self.bc = BoundaryCondition(self.x, y=Q_(0, 'm'))


@dataclass
class Hinge(Support):

    def __post_init__(self):
        self.bc = BoundaryCondition(self.x, y=Q_(0, 'm'))


@dataclass
class FixedEnd(Support):

    def __post_init__(self):
        self.bc = BoundaryCondition(self.x, y=Q_(0, 'm'), theta=Q_(0, 'rad'))


@dataclass
class FreeEnd(Support):
    pass


class Beam(statics.Beam):
    """Extends class `Beam` from subpackage `statics`.
    Base class `Beam` in subpackage `statics` implements calculation of the
    resultant internal force and/or moment in a cross-section of a beam.
    The extended class `Beam` adds geometry to the cross-section by
    encapsulating an instance of class `Section` (see `strength.section.py`).
    Through this instance it is possible to determine stresses resulting from the
    internal force and/or moment in a given cross-section.
    Furthermore, this class can also be used to determine the deformation of
    the beam: elongation if the beam is subject to axial loadings, angle of
    twist if the beam is subject to torsional loadings. If the beam is subject
    to bending, then through its attribute `elastic_curve` the methods of class
    `ElasticCurve` can be used to determine the vertical displacement and
    slope of the beam (see `strength.strain.deflection.py`).
    """
    def __init__(
        self,
        length: Quantity,
        shape_type: Type[Shape],
        shape_dim: Dimensions,
        E_modulus: Quantity,
        G_modulus: Quantity,
        loadings: list[Force | Moment | DistributedLoad1D],
        supports: list[Support] | None = None,
        units: tuple[str, str] | None = None,
        num_sections: int = 50
    ) -> None:
        """Creates a `Beam` instance.

        Parameters
        ----------
        length:
            Length of the beam.
        shape_type:
            Type of shape, derived from base class `Shape`.
        shape_dim:
            Instance of class `Dimensions` holding the dimensions of the shape.
        E_modulus:
            Modulus of elasticity or Young's modulus of the material.
        G_modulus:
            Shear modulus of elasticity or modulus of rigidity of the material.
        loadings:
            List of the external loadings applied to the beam. If some loadings
            are not fully determined yet (reaction forces or moments), these
            will be solved on instantiation of this class.
        supports:
            Optional list with the supports that support the beam. If the
            elastic curve of the beam is to be determined, the supports must be
            specified (instances of class `Roller`, `Hinge`, `FixedEnd`, and/or
            `FreeEnd`). If one end of the beam is not supported, this must also
            be indicated by a `FreeEnd` instance.
        units:
            Tuple with the units of force and the units of length to be used.
        num_sections:
            The number of sections to be made for determining the normal-force
            diagram, shear diagram and moment diagram of the beam.
        """
        self.section = None  # (*)
        super().__init__(length, loadings, units, num_sections)
        self.supports = supports
        self.section = Section(shape_type, shape_dim)
        self.E = E_modulus
        self.G = G_modulus
        self._axial_strain: AxialStrain | None = None
        self._torsional_strain: TorsionalStrain | None = None
        self.elastic_curve: ElasticCurve | None = None
        if np.any(self._N_arr):
            self._axial_strain = AxialStrain(
                self.N,
                self.section.shape.area,
                self.E
            )
        if np.any(self._T_arr):
            self._torsional_strain = TorsionalStrain(
                self.T,
                self.section.shape.polar_moment_of_inertia,
                self.G
            )
        if np.any(self._M_arr) and self.supports is not None:
            self.elastic_curve = ElasticCurve(
                self.M,
                self.E,
                self.section.shape.moment_of_inertia_xx,
                [sup.bc for sup in self.supports],
                self.num_sections
            )

    # (*): On instantiation of the derived class `Beam` when the line
    # `super().__init__(length, loadings, units, num_sections)` is executed, the
    # derived method `cut` from this class (instead of the original `cut` method
    # in the base class) is called through `self` inside the `__init__()` of
    # base class `Beam`. However, as the base class has no attribute `section`,
    # this would raise an exception. By setting `self.section` to `None` and
    # check for it in the derived method `cut`, this issue has been circumvented.

    def cut(
        self,
        x: Quantity,
        view: str = 'left'
    ) -> tuple[Force, Moment]:
        """Returns the internal force and/or moment at the position `x` along
        the longitudinal axis of the beam.
        Parameter `view` indicates whether the cross-section of the left part of
        the cut beam is to be regarded (this is the default) or the
        cross-section of its right part (`view = 'right'`).

        Through attribute `section` of this class the stresses resulting from
        the internal force and/or moment can also be determined.

        Returns
        -------
        Tuple with the internal force (`Force`-object) and internal moment
        (`Moment`-object) acting at the viewed cross-section.
        """
        F_i, M_i = super().cut(x, view)
        if self.section is not None:  # (*)
            self.section.set_internal_loadings(F_i, M_i)
        return F_i, M_i

    def elongation(self, x1: Quantity, x2: Quantity) -> Quantity:
        """Returns the elongation of the beam between its longitudinal positions
        `x1` and `x2`, where `x2` must be greater than `x1`.
        """
        if self._axial_strain is not None:
            return self._axial_strain.elongation(x1, x2)
        else:
            raise ValueError("no axial loadings are acting on the beam")

    def angle_of_twist(self, x1: Quantity, x2: Quantity) -> Quantity:
        """Returns the angle of twist of the circular shaft between its
        longitudinal positions `x1` and `x2`, where `x2` must be greater than
        `x1`.
        """
        if self._torsional_strain is not None:
            return self._torsional_strain.angle_of_twist(x1, x2)
        else:
            raise ValueError("no torques are acting on the beam")


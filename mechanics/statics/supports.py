from abc import ABC
from enum import IntEnum, StrEnum
import numpy as np
import sympy as sp
from mechanics import Quantity
from ._position import Position, set_position


Q_ = Quantity


class Support(ABC):
    """Abstract base class that represents a generic beam support that applies
    unknown reaction forces to a beam.
    """

    def __init__(self, name: str, position: Position):
        """Creates a `Support` instance.

        Parameters
        ----------
        name:
            Identifies the support of the beam.
        position:
            Position of the support along the x-axis of the beam.
        """
        self.name = name
        self.position = set_position(position)
        self.__create_symbols()
        self.F: sp.Expr | None = None
        self.theta: Quantity | None = None
        self.F_x: sp.Expr | None = None
        self.F_y: sp.Expr | None = None
        self.F_z: sp.Expr | None = None
        self.M_x: sp.Expr | None = None
        self.M_y: sp.Expr | None = None
        self.M_z: sp.Expr | None = None

    def __create_symbols(self) -> None:
        self._F = sp.Symbol(f'{self.name}.F')
        self._F_x = sp.Symbol(f'{self.name}.F_x')
        self._F_y = sp.Symbol(f'{self.name}.F_y')
        self._F_z = sp.Symbol(f'{self.name}.F_z')
        self._M_x = sp.Symbol(f'{self.name}.M_x')
        self._M_y = sp.Symbol(f'{self.name}.M_y')
        self._M_z = sp.Symbol(f'{self.name}.M_z')


class Roller2D(Support):
    """Represents a planar roller. This support has one unknown reaction force
    that is always perpendicular to the x-axis of the beam (i.e. directed along
    the y-axis).
    """
    def __init__(self, name: str, position: Position):
        super().__init__(name, position)
        self.F_y = self._F_y


class Hinge2D(Support):
    """Represents a planar hinge. This support has two unknown reaction force
    components. One component is directed along the x-axis of the beam. The
    other component is perpendicular to the x-axis (i.e. directed along the
    y-axis).
    """
    def __init__(self, name: str, position: Position):
        super().__init__(name, position)
        self.F_x = self._F_x
        self.F_y = self._F_y


class FixedEnd2D(Support):
    """Represents a planar fixed end of a beam. This support has two unknown
    reaction force components, just like a `Hinge2D` object, and also a
    unknown reaction bending moment of which the vector is perpendicular to
    the plane of the beam (i.e. the xy-plane).
    """
    def __init__(self, name: str, position: Position):
        super().__init__(name, position)
        self.F_x = self._F_x
        self.F_y = self._F_y
        self.M_z = self._M_z


class TwoForceMember2D(Support):
    """Represents a planar two-force-member, e.g. a slender rod or a cable. The
    reaction force exerted by this support has only one unknown: its magnitude.
    The direction of the reaction force coincides with the longitudinal axis of
    the two-force-member.
    """
    def __init__(self, name: str, position: Position, theta: Quantity):
        super().__init__(name, position)
        self.theta = theta
        cos_theta = round(np.cos(theta.to('rad').m), 12)
        sin_theta = round(np.sin(theta.to('rad').m), 12)
        self.F_x = self._F * cos_theta
        self.F_y = self._F * sin_theta

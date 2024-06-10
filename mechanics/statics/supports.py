from abc import ABC
from enum import IntEnum, StrEnum
import sympy as sp
from mechanics import Quantity
from .position import Position, set_position


Q_ = Quantity


class Sign(IntEnum):
    NEGATIVE = -1
    POSITIVE = 1
    ZERO = 0


class Component(StrEnum):
    F_x = 'F_x'
    F_y = 'F_y'
    F_z = 'F_z'
    M_x = 'M_x'
    M_y = 'M_y'
    M_z = 'M_z'


class Support(ABC):
    """Abstract base class that represents a generic beam support that applies
    unknown reaction forces to a beam.
    """

    def __init__(self, name: str, position: Position, **kwargs):
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
        self.__create_expressions(**kwargs)
        self.F_x: sp.Expr | None = None
        self.F_y: sp.Expr | None = None
        self.F_z: sp.Expr | None = None
        self.M_x: sp.Expr | None = None
        self.M_y: sp.Expr | None = None
        self.M_z: sp.Expr | None = None

    def __create_symbols(self) -> None:
        self._sym_F_x = sp.Symbol(f'{self.name}.F_x')
        self._sym_F_y = sp.Symbol(f'{self.name}.F_y')
        self._sym_F_z = sp.Symbol(f'{self.name}.F_z')
        self._sym_M_x = sp.Symbol(f'{self.name}.M_x')
        self._sym_M_y = sp.Symbol(f'{self.name}.M_y')
        self._sym_M_z = sp.Symbol(f'{self.name}.M_z')

    def __create_expressions(self, **kwargs) -> None:
        self._F_x = kwargs.get('F_x', Sign.POSITIVE) * self._sym_F_x
        self._F_y = kwargs.get('F_y', Sign.POSITIVE) * self._sym_F_y
        self._F_z = kwargs.get('F_z', Sign.POSITIVE) * self._sym_F_z
        self._M_x = kwargs.get('M_x', Sign.POSITIVE) * self._sym_M_x
        self._M_y = kwargs.get('M_y', Sign.POSITIVE) * self._sym_M_y
        self._M_z = kwargs.get('M_z', Sign.POSITIVE) * self._sym_M_z


class Roller2D(Support):
    """Represents a planar roller. This support has one unknown reaction force
    that is always perpendicular to the x-axis of the beam (i.e. directed along
    the y-axis).
    """
    def __init__(self, name: str, position: Position, **kwargs):
        super().__init__(name, position, **kwargs)
        self.F_y = self._F_y


class Hinge2D(Support):
    """Represents a planar hinge. This support has two unknown reaction force
    components. One component is directed along the x-axis of the beam. The
    other component is perpendicular to the x-axis (i.e. directed along the
    y-axis).
    """
    def __init__(self, name: str, position: Position, **kwargs):
        super().__init__(name, position, **kwargs)
        self.F_x = self._F_x
        self.F_y = self._F_y


class FixedEnd2D(Support):
    """Represents a planar fixed end of a beam. This support has two unknown
    reaction force components, just like a `Hinge2D` object, and also a
    unknown reaction bending moment of which the vector is perpendicular to
    the plane of the beam (i.e. the xy-plane).
    """
    def __init__(self, name: str, position: Position, **kwargs):
        super().__init__(name, position, **kwargs)
        self.F_x = self._F_x
        self.F_y = self._F_y
        self.M_z = self._M_z

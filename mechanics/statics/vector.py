from __future__ import annotations
from dataclasses import dataclass
from typing import Any
import numpy as np
import sympy as sp
from scipy.interpolate import interp1d
from scipy.integrate import quad
from mechanics import Quantity, UNITS
from mechanics.geometry import AxesRotation2D


Q_ = Quantity


class Position:
    """Represents a position in 3D space.

    Generally a position is defined by its x, y, and z-coordinate. The
    coordinates have units of length (default is meter).
    """
    num_decimals: int = 9

    def __init__(
        self,
        x: float = 0.0,
        y: float = 0.0,
        z: float = 0.0,
        units: str = 'm'
    ) -> None:
        self._x = x
        self._y = y
        self._z = z
        self.units = units
        self._pos = Q_([self._x, self._y, self._z], self.units)

    def to(self, units: str) -> Position:
        """Converts the position to the given units of length and then returns
        it."""
        self.units = units
        self._pos = self._pos.to(units)
        self._x = self._pos[0].m
        self._y = self._pos[1].m
        self._z = self._pos[2].m
        return self

    def __repr__(self) -> str:
        x, y, z = tuple(map(
            lambda coord: round(coord, self.num_decimals),
            [coord for coord in (self._x, self._y, self._z)]
        ))
        s = f"({x}; {y}; {z}) {self.units}"
        return s

    @property
    def x(self) -> Quantity:
        """Returns the x-coordinate of the position as a `Quantity` object."""
        return Q_(self._x, self.units)

    @property
    def y(self) -> Quantity:
        """Returns the y-coordinate of the position as a `Quantity` object."""
        return Q_(self._y, self.units)

    @property
    def z(self) -> Quantity:
        """Returns the z-coordinate of the position as a `Quantity` object."""
        return Q_(self._z, self.units)

    @property
    def as_quantity(self) -> Quantity:
        """Returns the position as a `Quantity` object of the form
        `Quantity([x, y, z], <unit_str>).
        """
        return Q_([self._x, self._y, self._z], self.units)


ORIGIN = Position()
"""Defines the origin (0, 0, 0) of the coordinate system."""


@dataclass
class Angle:
    """Defines an angle. An angle is defined by its magnitude and units.
    Default units are decimal degrees. Counterclockwise angles should have a
    positive magnitude, while clockwise angles should be indicated by a negative
    value.
    """
    magnitude: float = 0.0
    units: str = 'deg'

    def __post_init__(self):
        self._angle = Q_(self.magnitude, self.units)

    @classmethod
    def create(
        cls,
        v: float | Quantity,
        h: float | Quantity,
        quadrant: int = 1
    ) -> Angle:
        """Creates an angle given the slope specified by a vertical distance `v`
        and horizontal distance `h`. By specifying the quadrant, the slope angle
        can be positioned in a right-handed coordinate system. The first
        quadrant is situated between 0° and 90° (where both the x- and
        y-coordinates have positive values). The second quadrant is situated
        between 90° and 180° (where the x-coordinates have a negative and the
        y-coordinates have a positive value). The third quadrant is situated
        between 180° and 270° (where both the x- and y-coordinates have negative
        values). The fourth quadrant is situated between 270° and 360° (where
        the x-coordinates have a positive value and the y-coordinates have a
        negative value).
        """
        if isinstance(v, Quantity) and isinstance(h, Quantity):
            v = v.m
            h = h.to(v.units).m
        a = np.arctan2(v, h)
        match quadrant:
            case 1:
                return cls(a, 'rad')
            case 2:
                return cls(np.pi - a, 'rad')
            case 3:
                return cls(a + np.pi, 'rad')
            case 4:
                return cls(-a, 'rad')

    @property
    def degrees(self) -> float:
        """Returns the magnitude of the angle in degrees."""
        return self._angle.to('deg').m

    @property
    def radians(self) -> float:
        """Returns the magnitude of the angle in radians."""
        return self._angle.to('rad').m

    @property
    def deg(self) -> float:
        """Returns the magnitude of the angle in degrees."""
        return self.degrees

    @property
    def rad(self) -> float:
        """Returns the magnitude of the angle in radians."""
        return self.radians

    def __add__(self, other: Angle) -> Angle:
        """Returns the addition of angle `self` and angle `other`."""
        angle = Angle(self.degrees + other.degrees)
        return angle

    @property
    def as_quantity(self) -> Quantity:
        """Returns the angle as a `Quantity` object."""
        return self._angle


class Vector:
    """Defines a 3D (or 2D) vector, of which the magnitude and/or direction
    angles can be still undetermined, i.e. a symbolic vector.
    """
    num_decimals: int = 9
    default_units: str = UNITS.dimensionless

    def __init__(
        self,
        magnitude: float | str = 0.0,
        theta: Angle | str = Angle(0),
        gamma: Angle | str = Angle(0),
        position: Position = ORIGIN,
        units: str | None = None,
        name: str | None = None
    ) -> None:
        """Creates a `Vector` object.

        Parameters
        ----------
        magnitude:
            Magnitude of the vector, which is always a positive value.
            If the magnitude is unknown, a symbolic name (string) to indicate the
            magnitude of the vector.
        theta:
            Smallest angle measured from the positive x-axis towards the
            projection of the vector on the xy-plane. Counterclockwise angles
            are positive, while clockwise angles are negative.
            If the angle `theta` is unknown, a symbolic name (string) to
            indicate this angle.
        gamma:
            Smallest angle measured from the projection of the vector on the
            xy-plane towards the action line of the vector. Counterclockwise
            angles are positive, while clockwise angles are negative. In case
            of a 2D vector, angle `gamma` is always zero (i.e. the default
            value).
            If the angle `gamma` is unknown, a symbolic name (string) to
            indicate this angle.
        position:
            Start point of the vector. Default is the origin of the coordinate
            system.
        units:
            In case of a vector quantity, the units must be specified, e.g.
            a force in units of Newton. If `None`, the vector is dimensionless.
        name:
            Name to identify the vector.
            If no name is given and the magnitude is a string, the name is set
            to this string. Otherwise, the name is set to the id of the `Vector`
            object (returned by Python's built-in function `id()`).
        """
        self._magnitude = magnitude
        self.theta = theta
        self.gamma = gamma
        self.position = position
        self.units = units if units is not None else self.default_units
        self.name = self.__set_name(name)

        self._vector = self.__create_vector()

    def __set_name(self, name: str | None):
        if name is None and isinstance(self._magnitude, str):
            return self._magnitude
        elif isinstance(name, str):
            return name
        else:
            return str(id(self))

    def is_symbolic(self) -> bool:
        """Returns `True` if magnitude and or direction angles are unknown."""
        if isinstance(self._magnitude, str):
            return True
        elif isinstance(self.theta, str):
            return True
        elif isinstance(self.gamma, str):
            return True
        else:
            return False

    def __create_vector(self) -> Quantity | tuple[sp.Expr, ...]:
        if not self.is_symbolic():
            return self.__quantity_vector()
        else:
            return self.__symbolic_vector()

    def __quantity_vector(self) -> Quantity:
        # Decomposes the vector into its cartesian components x, y, z.
        # Internally the vector is represented as a `Quantity` object.
        vec_xy = abs(self._magnitude * np.cos(self.gamma.rad))
        vec_x = round(vec_xy * np.cos(self.theta.rad), self.num_decimals)
        vec_y = round(vec_xy * np.sin(self.theta.rad), self.num_decimals)
        vec_z = round(self._magnitude * np.sin(self.gamma.rad), self.num_decimals)
        vec = Q_([vec_x, vec_y, vec_z], self.units)
        return vec

    def __symbolic_vector(self) -> tuple[sp.Expr, ...]:
        # Decomposes the vector into its cartesian components x, y, z.
        # Internally the vector is represented as a tuple of Sympy expressions.
        magnitude = sp.Symbol(self._magnitude) if isinstance(self._magnitude, str) else self._magnitude
        gamma = sp.Symbol(self.gamma) if isinstance(self.gamma, str) else self.gamma.rad
        theta = sp.Symbol(self.theta) if isinstance(self.theta, str) else self.theta.rad
        vec_xy = None
        vec_x = sp.Symbol(f"{self.name}.x")
        vec_y = sp.Symbol(f"{self.name}.y")
        vec_z = sp.Symbol(f"{self.name}.z")
        if isinstance(gamma, float):
            vec_xy = magnitude * round(np.cos(gamma), self.num_decimals)
            vec_z = magnitude * round(np.sin(gamma), self.num_decimals)
        if isinstance(theta, float):
            if vec_xy is None:
                vec_x = vec_x * round(np.cos(theta), self.num_decimals)
                vec_y = vec_y * round(np.sin(theta), self.num_decimals)
            else:
                vec_x = vec_xy * round(np.cos(theta), self.num_decimals)
                vec_y = vec_xy * round(np.sin(theta), self.num_decimals)
        return vec_x, vec_y, vec_z

    def to(self, units: str) -> Vector:
        """Converts the vector to the given units and then returns it."""
        self.units = units
        if not self.is_symbolic():
            self._vector = self._vector.to(units)
        return self

    @property
    def x(self) -> Quantity | sp.Expr:
        """Returns the x-component of the vector as a `Quantity` object or Sympy
        expression.
        """
        if isinstance(self._vector[0], Quantity):
            return self._vector[0]
        else:
            return self._vector[0]

    @property
    def y(self) -> Quantity | sp.Expr:
        """Returns the y-component of the vector as a `Quantity` object or Sympy
        expression."""
        if isinstance(self._vector[1], Quantity):
            return self._vector[1]
        else:
            return self._vector[1]

    @property
    def z(self) -> Quantity | sp.Expr:
        """Returns the z-component of the vector as a `Quantity` object or Sympy
        expression."""
        if isinstance(self._vector[2], Quantity):
            return self._vector[2]
        else:
            return self._vector[2]

    @property
    def component_values(self) -> tuple[float | sp.Expr, ...]:
        """Returns the three cartesian components x, y, z of the vector as
        floats or as a Sympy expressions."""
        x = self.x.m if isinstance(self.x, Quantity) else self.x
        y = self.y.m if isinstance(self.y, Quantity) else self.y
        z = self.z.m if isinstance(self.z, Quantity) else self.z
        return x, y, z

    @property
    def components(self) -> tuple[Quantity | sp.Expr, ...]:
        """Returns the three cartesian components x, y, z of the vector as
        `Quantity` objects or as a Sympy expressions."""
        x = self.x if isinstance(self.x, Quantity) else self.x
        y = self.y if isinstance(self.y, Quantity) else self.y
        z = self.z if isinstance(self.z, Quantity) else self.z
        return x, y, z

    def __repr__(self) -> str:
        x, y, z = self.component_values
        s = f"<"
        l = [('x', x), ('y', y), ('z', z)]
        for i, tup in enumerate(l):
            if isinstance(tup[1], float):
                s += f"{tup[0]}: {round(tup[1], self.num_decimals)}"
            else:
                s += f"{tup[0]}: {str(tup[1])}"
            if i == 2:
                s += ">"
            else:
                s += "; "
        s += f" {self.units}"
        return s

    @classmethod
    def reverse(cls, vector: Vector, name: str | None = None) -> Vector:
        """Returns a new vector by reversing the given vector.

        Another name can be given to the new vector. If `None`, the name of the
        original vector is kept.

        Notes
        -----
        Only fully determined (i.e. non symbolic) vectors can be reversed. A
        `ValueError` exception is raised in case `vector` should be symbolic.
        """
        if isinstance(vector.theta, Angle) and isinstance(vector.gamma, Angle):
            theta = vector.theta + Angle(180)
            if vector.gamma.magnitude == 0.0:  # vector only in xy-plane
                gamma = vector.gamma
            else:
                gamma = vector.gamma + Angle(180)
            reversed_vector = cls(
                magnitude=vector._magnitude,
                theta=theta,
                gamma=gamma,
                position=vector.position,
                units=vector.units,
                name=name if name is not None else vector.name
            )
            return reversed_vector
        else:
            raise ValueError("a symbolic vector cannot be reversed")

    @classmethod
    def rotate(
        cls,
        vector: Vector,
        theta: Angle,
        gamma: Angle = Angle(0),
        name: str | None = None
    ) -> Vector:
        """Returns a new vector by rotating the given vector.

        Another name can be given to the new vector. If `None`, the name of the
        original vector is kept.

        Notes
        -----
        Only fully determined (i.e. non symbolic) vectors can be reversed. A
        `ValueError` exception is raised in case `vector` should be symbolic.
        """
        if isinstance(vector.theta, Angle) and isinstance(vector.gamma, Angle):
            theta = vector.theta + theta
            if vector.gamma.magnitude == 0.0:  # vector only in xy-plane
                gamma = vector.gamma
            else:
                gamma = vector.gamma + gamma
            rotated_vector = cls(
                magnitude=vector._magnitude,
                theta=theta,
                gamma=gamma,
                position=vector.position,
                units=vector.units,
                name=name if name is not None else vector.name
            )
            return rotated_vector
        else:
            raise ValueError("a symbolic vector cannot be rotated")

    @classmethod
    def create_from_components(
        cls,
        vec_x: float | None,
        vec_y: float | None,
        vec_z: float | None,
        position: Position = ORIGIN,
        units: str | None = None,
        name: str | None = None
    ) -> Vector:
        """Creates a `Vector` object from the cartesian components of the
        vector. If any of the components is unknown, set this component to
        `None`.

        Parameters
        ----------
        vec_x:
            x-component of the vector.
        vec_y:
            y-component of the vector.
        vec_z:
            z-component of the vector.
        position:
            Start point of the vector. Default is the origin of the coordinate
            system.
        units:
            In case of a vector quantity, the units must be specified, e.g.
            a force in units of Newton. If `None`, the vector is dimensionless.
        name:
            Name to identify the vector.
            If no name is given and any of the components is unknown, the name
            of the vector is set to this string. Otherwise, the name is set to
            the id of the `Vector` object (returned by Python's built-in
            function `id()`).
        """
        if all([isinstance(c, float) for c in (vec_x, vec_y, vec_z)]):
            mag = np.sqrt(vec_x**2 + vec_y**2 + vec_z**2)
            mag_xy = np.sqrt(vec_x**2 + vec_y**2)
            theta = Angle(np.arctan2(vec_y, vec_x), 'rad')
            gamma = Angle(np.arctan2(vec_z, mag_xy), 'rad')
            vec = cls(mag, theta, gamma, position, units, name)
            return vec
        else:
            vec = cls('unknown', 'theta', 'gamma', position, units, name)
            if vec_x is None: vec_x = sp.Symbol(f"{vec.name}.x")
            if vec_y is None: vec_y = sp.Symbol(f"{vec.name}.y")
            if vec_z is None: vec_z = sp.Symbol(f"{vec.name}.z")
            vec._vector = (vec_x, vec_y, vec_z)
            return vec

    @property
    def magnitude(self) -> Quantity:
        return Q_(self._magnitude, self.units)


class Force(Vector):
    """
    Defines a force as a derived class from class `Vector`.
    """
    default_units = 'N'

    def moment(self, point: Position = ORIGIN) -> Moment:
        """Returns the moment of the force about the given point."""
        point = point.to(self.position.units)
        x = self.position.x - point.x
        y = self.position.y - point.y
        z = self.position.z - point.z
        F_x, F_y, F_z = self.component_values
        M_x = -z * F_y + y * F_z
        M_y = z * F_x - x * F_z
        M_z = -y * F_x + x * F_y
        units = f"{self.units} * {self.position.units}"
        moment = Moment.create_from_components(
            M_x, M_y, M_z,
            point, units, f"moment of {self.name}"
        )
        return moment


class Moment(Vector):
    """
    Defines a moment or torque as a derived class from class `Vector`.
    """
    default_units = 'N * m'


class DistributedLoad1D:
    """
    Defines a distributed linear (1D) load.
    """
    def __init__(
        self,
        x_coords: Quantity,
        loads: Quantity,
        name: str | None = None,
        slope: Angle = Angle(0)
    ) -> None:
        """Creates a `DistributedLoad1D` object.

        Parameters
        ----------
        x_coords:
            Positions in ascending order along a straight axis (x-axis) where
            the specific loads (see below) are given. These x-coordinates are
            passed as a `Quantity` object. E.g. `Quantity([1, 2, 3], 'm')` are
            the x-coordinates at 1, 2, and 3 meter from the origin along the
            x-axis.
        loads:
            Specific loads (in units of force per units of length) at each of the
            given positions in the same order. These specific loads are also
            passed as a `Quantity` object. E.g. `Quantity([50, 100, 60], 'N/m')`
            means that at x = 1 m the specific load is 50 N/m,
            at x = 2 m 100 N/m, and at x = 3 m 60 N/m.
        name:
            Name to identify the distributed load.
            If no name is given the name is set to the id of the
            `DistributedLoad1D` object (returned by Python's built-in function
            `id()`).
        slope:
            Slope angle of the straight axis where the first position in
            `x_coords` is taken as the pivot point.

        Notes
        -----
        A distributed linear load must be fully determined; it cannot have any
        elements that are still unknown.
        """
        self.x_coords = x_coords
        self.loads = loads
        self.name = self.__set_name(name)
        self.slope = slope

        self._units_of_length = f"{self.x_coords.units:~P}"
        self._units_of_load = f"{self.loads.units:~P}"
        self._units_of_force = self.__get_units_of_force()

        self._q: Any = interp1d(self.x_coords.m, self.loads.m)

    def __set_name(self, name: str | None):
        if isinstance(name, str):
            return name
        else:
            return str(id(self))

    def __get_units_of_force(self) -> str:
        zero_force = (
            Q_(0, self._units_of_load)
            * Q_(0, self._units_of_length)
        )
        return f"{zero_force.units:~P.0f}"

    def to(self, units_of_load: str) -> DistributedLoad1D:
        """Converts the distributed linear load to the given units and then
        returns it.
        """
        self.loads = self.loads.to(units_of_load)
        self._units_of_load = f"{self.loads.units:~P}"
        self._units_of_force = self.__get_units_of_force()
        self._q: Any = interp1d(self.x_coords.m, self.loads.m)
        return self

    def positions(self, units_of_length: str) -> Quantity:
        """Converts the positions of the distributed linear load to the given
        units and then returns them.
        """
        self.x_coords = self.x_coords.to(units_of_length)
        self._units_of_length = f"{self.x_coords.units:~P}"
        self._units_of_force = self.__get_units_of_force()
        self._q: Any = interp1d(self.x_coords.m, self.loads.m)
        return self.x_coords

    def resultant(
        self,
        x1: Quantity | float | None = None,
        x2: Quantity | float | None = None
    ) -> Force:
        """Returns the resultant force between the positions `x1` and `x2` where
        `x2` > `x1`. If `x1` is not specified the first position in `x_coords`
        will be taken. If `x2` is not specified the last position in `x_coords`
        will be taken.
        """
        x1 = x1 if x1 is not None else self.x_coords[0]
        x2 = x2 if x2 is not None else self.x_coords[-1]

        if isinstance(x1, Quantity):
            x1 = x1.to(self._units_of_length).m
        if isinstance(x2, Quantity):
            x2 = x2.to(self._units_of_length).m

        x1 = min(x1, x2)
        x2 = max(x1, x2)
        if x2 > x1:
            Q_mag = quad(self._q, x1, x2)[0]
            x_c = quad(lambda x: x * self._q(x), x1, x2)[0] / Q_mag
            theta = Angle(90) if Q_mag >= 0 else Angle(-90)
            y_c = 0.0
            if self.slope.magnitude != 0.0:
                theta += self.slope
                x_0 = self.x_coords[0].m
                p_c = Q_([x_c - x_0, 0], self._units_of_length)
                rotator = AxesRotation2D(-self.slope.as_quantity)
                p_c = rotator(p_c)
                x_c = x_0 + p_c[0].to(self._units_of_length).m
                y_c = p_c[1].to(self._units_of_length).m
            Q = Force(
                magnitude=abs(Q_mag),
                theta=theta,
                position=Position(x=x_c, y=y_c, units=self._units_of_length),
                units=self._units_of_force,
                name=f"resultant of {self.name}"
            )
            return Q
        else:
            raise ValueError("position `x2` must be further than position `x1`")

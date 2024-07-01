"""
Some useful functions and classes...
"""
import warnings
import numpy as np
from mechanics import Quantity

warnings.filterwarnings("error", category=RuntimeWarning)

Q_ = Quantity


def create_angle(
    v: float | Quantity,
    h: float | Quantity,
    quadrant: int = 1
) -> Quantity:
    """Creates an angle (instance of class `Quantity`) given the slope
    specified by a vertical distance `v` and horizontal distance `h`. By
    specifying the quadrant, the slope angle can be positioned in a right-handed
    (x,y) coordinate system (the x-axis runs horizontally pointing to the right
    and the y-axis runs vertically pointing upward). The first quadrant is
    situated between 0° and 90° (where both the x- and y-coordinates have
    positive values). The second quadrant is situated between 90° and 180°
    (where the x-coordinates have a negative value and the y-coordinates have a
    positive value). The third quadrant is situated between 180° and 270° (where
    both the x- and y-coordinates have negative values). The fourth quadrant is
    situated between 270° and 360° (where the x-coordinates have a positive
    value and the y-coordinates have a negative value).
    """
    if isinstance(v, Quantity) and isinstance(h, Quantity):
        v = v.m
        h = h.to(v.units).m
    a = Q_(np.arctan2(v, h), 'rad')
    match quadrant:
        case 1:
            return a
        case 2:
            return np.pi - a
        case 3:
            return a + np.pi
        case 4:
            return -a


class AxesRotation2D:
    """
    Implements the transformation of coordinates due to rotation of a coordinate
    system around the same origin as the original coordinate system.
    """
    def __init__(self, angle: Quantity) -> None:
        """Creates a `AxesRotation2D` object.

        Parameters
        ----------
        angle:
            Angle between original coordinate system and rotated coordinate
            system.
        """
        self.angle = angle

    def __call__(
        self,
        vector2D: Quantity
    ) -> Quantity:
        """Returns the components of a 2D vector in another coordinate system
        which is rotated by a given angle with respect to the original
        coordinate system.
        """
        x = vector2D[0]
        y = vector2D[1]
        x_rot = x * np.cos(self.angle) + y * np.sin(self.angle)
        y_rot = -x * np.sin(self.angle) + y * np.cos(self.angle)
        vector2D_rot = Quantity.from_list([x_rot, y_rot])
        return vector2D_rot


class Line:

    def __init__(self):
        self.m = None  # slope
        self.q = None  # intersection with the y-axis (0, q)

    @classmethod
    def from_point_and_angle(
        cls,
        p: tuple[float, float],
        angle: float
    ) -> 'Line':
        """Create a straight line (`Line` object) with slope `angle`
        (in radians) through point `p` (x,y).
        """
        line = cls()
        x1, y1 = p[0], p[1]
        m = np.tan(angle)
        q = y1 - m * x1
        line.m = m
        line.q = q
        return line

    @classmethod
    def from_uvw(
        cls,
        u: float,
        v: float,
        w: float
    ) -> 'Line':
        """Create a straight line (`Line` object) from the equation
        $u.x + v.y + w = 0$.
        """
        line = cls()
        line.m = -u / v
        line.q = -w / v
        return line

    @classmethod
    def from_two_points(
        cls,
        p1: tuple[float, float],
        p2: tuple[float, float]
    ) -> 'Line':
        """Create a straight line (`Line` object) through two points `p1` and
        `p2`.
        """
        line = cls()
        x1, x2 = p1[0], p2[0]
        y1, y2 = p1[1], p2[1]
        try:
            line.m = (y2 - y1) / (x2 - x1)
        except ZeroDivisionError:
            line.m = float('inf') if y2 > y1 else -float('inf')
        line.q = y1 - line.m * x1
        return line

    def distance(self, p: tuple[float, float]) -> float:
        """Returns the distance between the line and a given point `p`(x,y)."""
        u = -self.m
        v = 1.0
        w = -self.q
        x_p, y_p = p[0], p[1]
        num = abs(u * x_p + v * y_p + w)
        den = np.sqrt(u ** 2 + v ** 2)
        d = num / den
        return d

    def y(self, x: float) -> float:
        """Returns the y-coordinate that corresponds with x-coordinate `x`.
        If the line is vertical (x = constant), the y-coordinate is undetermined
        and float NaN is returned.
        """
        y = self.m * x + self.q
        return y

    def x(self, y: float) -> float:
        """Returns the x-coordinate that corresponds with y-coordinate `y`.
        If the line is horizontal (y = constant), the x-coordinate is
        undetermined and float NaN is returned.
        """
        try:
            x = (y - self.q) / self.m
        except ZeroDivisionError:
            x = float('nan')
        return x


class LineSegment:

    def __init__(
        self,
        p1: tuple[float, float],
        p2: tuple[float, float]
    ) -> None:
        """Creates a `LineSegment` object.

        Parameters
        ----------
        p1:
            Start point of the line segment.
        p2:
            End point of the line segment.
        """
        self._line = Line.from_two_points(p1, p2)
        self.p1 = p1
        self.p2 = p2
        self.x1, self.y1 = p1[0], p1[1]
        self.x2, self.y2 = p2[0], p2[1]
        self.x_min = min(self.x1, self.x2)
        self.x_max = max(self.x1, self.x2)
        self.y_min = min(self.y1, self.y2)
        self.y_max = max(self.y1, self.y2)

    def y(self, x: float) -> float | None:
        """Returns the y-coordinate that corresponds with x-coordinate `x` or
        `None` if no point with x-coordinate `x` lies on the line segment.

        Notes
        -----
        If the line segment runs vertically (x = constant), the y-coordinate
        of a point on this line segment cannot be determined and NaN-float is
        returned.
        """
        if self.x_min <= x <= self.x_max:
            y = self._line.y(x)
            if np.isnan(y):
                return self.y1
            return y
        return None

    def x(self, y: float) -> float | None:
        """Returns the x-coordinate that corresponds with y-coordinate `y` or
        `None` if no point with y-coordinate `y` lies on the line segment.
        """
        if self.y_min <= y <= self.y_max:
            x = self._line.x(y)
            if np.isnan(x):
                return self.x1
            return x
        return None

    def is_vertical(self) -> bool:
        """Returns `True` if the line segment is vertical."""
        if np.isinf(self._line.m):
            return True
        return False

    def is_horizontal(self) -> bool:
        """Returns `True` if the line segment is horizontal."""
        if self._line.m == 0.0:
            return True
        return False

    def contains(self, point: tuple[float, float]) -> bool:
        """Returns `True` if the point lies on the line segment, otherwise
        returns `False`.
        """
        x = self.x(y=point[1])  # x-coord on line segment that corresponds with y-coord of point
        y = self.y(x=point[0])  # y-coord on line segment that corresponds with x-coord of point
        cond1 = x is not None and x == point[0]
        cond2 = y is not None and y == point[1]
        if not self.is_vertical() and cond1 and cond2:
            return True
        elif self.is_vertical() and cond1:
            return True
        elif self.is_horizontal() and cond2:
            return True
        return False

    @property
    def slope(self) -> float:
        """Returns the slope of the line segment."""
        return self._line.m

    @property
    def length(self) -> float:
        """Returns the length of the line segment."""
        if self.slope == 0:
            l = abs(self.x2 - self.x1)
            return l
        elif np.isinf(self.slope):
            l = abs(self.y2 - self.y1)
        else:
            dx = abs(self.x2 - self.x1)
            dy = abs(self.y2 - self.y1)
            l = np.sqrt(dx ** 2 + dy ** 2)
        return l

    def coincide(self, other: 'LineSegment') -> bool:
        """Returns `True` if `LineSegment self` coincides with `LineSegment
        other`.
        """
        cond1 = self._line.m == other._line.m
        cond2 = self._line.q == other._line.q
        cond3 = other.x_min <= self.x_min <= other.x_max
        # cond4 = other.y_min <= self.y_min <= other.y_max
        cond5 = self.x_min <= other.x_min <= self.x_max
        # cond6 = self.y_min <= other.y_min <= self.y_max
        if (cond1 and cond2) and (cond3 or cond5):
            return True
        return False

    def reverse(self) -> 'LineSegment':
        """Returns a new `LineSegment` with the start and end point being
        reversed.
        """
        return self.__class__(self.p2, self.p1)

    def __str__(self):
        return f"({self.x1}, {self.y1}) -> ({self.x2}, {self.y2})"

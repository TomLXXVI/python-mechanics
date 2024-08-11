from abc import ABC, abstractmethod
from dataclasses import dataclass
from collections.abc import Sequence
from math import sqrt, atan2, cos, sin, pi
from mechanics import Quantity
from mechanics.charts import LineChart


Q_ = Quantity


@dataclass
class Dimensions:
    """Groups the possible dimensions of any shape."""
    radius: Quantity | None = None
    diameter: Quantity | None = None
    width: Quantity | None = None
    height: Quantity | None = None
    inner_radius: Quantity | None = None
    outer_radius: Quantity | None = None
    thickness: Quantity | None = None
    web_thickness: Quantity | None = None
    flange_thickness: Quantity | None = None
    orientation: Quantity | None = None


class Shape(ABC):
    """Abstract base class from which concrete shapes are derived."""

    def __init__(self, dim: Dimensions) -> None:
        self.dim = dim

    @property
    @abstractmethod
    def area(self) -> Quantity:
        """Returns the area of the shape."""
        ...

    @property
    @abstractmethod
    def moment_of_inertia_xx(self) -> Quantity:
        """Returns the moment of inertia about the horizontal axis through
        the centroid of the shape.
        """
        ...

    @property
    @abstractmethod
    def moment_of_inertia_yy(self) -> Quantity:
        """Returns the moment of inertia about the vertical axis through the
        centroid of the shape.
        """
        ...

    @property
    @abstractmethod
    def product_of_inertia(self) -> Quantity:
        """Returns the product of inertia with respect to the horizontal and the
        vertical axis through the centroid of that shape.
        """
        ...

    @property
    @abstractmethod
    def polar_moment_of_inertia(self) -> Quantity:
        """Returns the polar moment of inertia of the shape about the centroid."""
        ...

    @property
    @abstractmethod
    def principal_moments_of_inertia(self) -> tuple[Quantity, Quantity, Quantity]:
        """Returns the moments of inertia about the major principal axis and the
        minor principal axis, and also the angle between the horizontal axis and
        the major principal axis.
        """
        ...


TPoints = Sequence[Quantity] | Sequence[tuple[Quantity, Quantity]]


class Polygon(Shape):

    def __init__(
        self,
        vertices: TPoints,
        dim: Dimensions | None = None
    ) -> None:
        """Creates a `Polygon` object.

        Parameters
        ----------
        vertices:
            A sequence of the vertices of the polygon in consecutive
            counter-clockwise order. The vertices are `Quantity` objects like
            `Q_([x, y], 'm')` with x and y the horizontal and vertical
            coordinate of a vertex referred to a orthogonal, right-handed
            coordinate system.
        dim:
            Main dimensions of the polygon.
        """
        super().__init__(dim)
        self._vertices = vertices
        self._units = self._vertices[0][0].units
        self._pts = self.__convert_to_float(vertices)
        if self._pts[0] != self._pts[-1]:
            self._pts = self._pts + self._pts[:1]
        self._x = [p[0] for p in self._pts]
        self._y = [p[1] for p in self._pts]
        self._A = self.__area()
        self._x_c, self._y_c = self.__centroid()
        self._I_xx, self._I_yy, self._I_xy = self.__inertia()
        self._J = self._I_xx + self._I_yy
        self._I_max, self._I_min, self._theta = self.__principal()

    @property
    def area(self) -> Quantity:
        """Area of the polygon."""
        return Q_(self._A, self._units ** 2)

    @property
    def centroid(self) -> tuple[Quantity, Quantity]:
        """Location of the centroid of the polygon. The coordinates are relative
        to the origin of the coordinate system that was used to specify the
        coordinates of the vertices of the polygon.
        """
        return Q_(self._x_c, self._units), Q_(self._y_c, self._units)

    @property
    def moment_of_inertia_xx(self) -> Quantity:
        """Moment of inertia about the horizontal x-axis of which the origin is
        in the centroid of the polygon (and the positive direction is pointing
        to the right).
        """
        return Q_(self._I_xx, self._units ** 4)

    @property
    def moment_of_inertia_yy(self) -> Quantity:
        """Moment of inertia about the vertical y-axis of which the origin is in
        the centroid of the polygon (and the positive direction is pointing
        upward).
        """
        return Q_(self._I_yy, self._units ** 4)

    @property
    def product_of_inertia(self) -> Quantity:
        """Product of inertia about the horizontal x-axis and the vertical
        y-axis of which the origin is in the centroid of the polygon.
        """
        return Q_(self._I_xy, self._units ** 4)

    @property
    def polar_moment_of_inertia(self) -> Quantity:
        """Polar moment of inertia about the centroid of the polygon."""
        return Q_(self._J, self._units ** 4)

    @property
    def principal_moments_of_inertia(self) -> tuple[Quantity, Quantity, Quantity]:
        """Tuple with the major principal moment of inertia, the minor principal
        moment of inertia, and the rotation angle of the principal axes with
        respect to the orthogonal right-handed coordinate system of which the
        origin is in the centroid of the polygon. The rotation angle is measured
        positive in the counter-clockwise sense from the horizontal x-axis to
        the major principal axis.
        """
        I_max = Q_(self._I_max, self._units ** 4)
        I_min = Q_(self._I_min, self._units ** 4)
        theta = Q_(self._theta, 'rad')
        return I_max, I_min, theta

    def plot(self) -> LineChart:
        """Returns a `LineChart` object with a drawing of the shape.

        Notes
        -----
        The cross marks the centroid of the polygon. The thick blue line
        indicates the orientation of the major principal axis, while the thin
        red line indicates the orientation of the minor principal axis.
        """
        x = [x - self._x_c for x in self._x]
        y = [y - self._y_c for y in self._y]
        x_min = min(x)
        x_max = max(x)
        y_min = min(y)
        y_max = max(y)
        b = .05 * max(x_max - x_min, y_max - y_min)
        l = min(x_max - x_min, y_max - y_min) / 10
        x_a1 = [
            -l * cos(self._theta),
            l * cos(self._theta)
        ]
        y_a1 = [
            -l * sin(self._theta),
            l * sin(self._theta)
        ]
        x_a2 = [
            -l * cos(self._theta + pi / 2),
            l * cos(self._theta + pi / 2)
        ]
        y_a2 = [
            -l * sin(self._theta + pi / 2),
            l * sin(self._theta + pi / 2)
        ]
        chart = LineChart()
        chart.add_xy_data(
            label='polygon',
            x1_values=x,
            y1_values=y,
            style_props={'color': 'black', 'lw': 2}
        )
        chart.add_xy_data(
            label='major_principal_axis',
            x1_values=x_a1,
            y1_values=y_a1,
            style_props={'color': '#0072B2', 'lw': 2}
        )
        chart.add_xy_data(
            label='minor_principal_axis',
            x1_values=x_a2,
            y1_values=y_a2,
            style_props={'color': '#D55E00'}
        )
        chart.axes.set_aspect('equal')
        chart.axes.set_xlim(xmin=x_min - b, xmax=x_max + b)
        chart.axes.set_ylim(ymin=y_min - b, ymax=y_max + b)
        return chart

    @property
    def vertices(self) -> list[Quantity]:
        """Returns the vertices of the polygon with the coordinates referred to
        the centroid of the polygon.
        """
        pts = [
            Q_([x - self._x_c, y - self._y_c], self._units)
            for x, y in zip(self._x[:-1], self._y[:-1])
        ]
        return pts

    def shift(self, x: Quantity, y: Quantity) -> None:
        x = x.to(self._units).magnitude
        y = y.to(self._units).magnitude
        dx = x - self._x_c
        dy = y - self._y_c
        self._x = [x + dx for x in self._x]
        self._y = [y + dy for y in self._y]
        self._pts = [(x, y) for x, y in zip(self._x, self._y)]
        self._x_c = self._x_c + dx
        self._y_c = self._y_c + dy

    @staticmethod
    def __convert_to_float(pts: TPoints) -> list[tuple[float, float]]:
        pts = [tuple(p) if not isinstance(p, tuple) else p for p in pts]
        pts = [(round(float(p[0].m), 12), round(float(p[1].m), 12)) for p in pts]
        return pts

    def __get_distance(self, i: int) -> float:
        x, y = self._x[i], self._y[i]
        dx = x - self._x_c
        dy = y - self._y_c
        d = sqrt(dx ** 2 + dy ** 2)
        return d

    def __area(self) -> float:
        s = 0
        for i in range(len(self._pts) - 1):
            s += self._x[i] * self._y[i+1] - self._x[i+1] * self._y[i]
        A = s / 2
        return A

    def __centroid(self) -> tuple[float, float]:
        sx = sy = 0
        for i in range(len(self._pts) - 1):
            k = self._x[i] * self._y[i + 1] - self._x[i + 1] * self._y[i]
            sx += k * (self._x[i] + self._x[i + 1])
            sy += k * (self._y[i] + self._y[i + 1])
        x_c = sx / (6 * self._A)
        y_c = sy / (6 * self._A)
        return x_c, y_c

    def __inertia(self) -> tuple[float, float, float]:
        sxx = syy = sxy = 0
        for i in range(len(self._pts) - 1):
            k = self._x[i] * self._y[i + 1] - self._x[i + 1] * self._y[i]
            sxx += k * (
                self._y[i] ** 2
                + self._y[i] * self._y[i + 1]
                + self._y[i + 1] ** 2
            )
            syy += k * (
                self._x[i] ** 2
                + self._x[i] * self._x[i + 1]
                + self._x[i + 1] ** 2
            )
            sxy += k * (
                self._x[i] * self._y[i + 1]
                + 2 * self._x[i] * self._y[i]
                + 2 * self._x[i + 1] * self._y[i + 1]
                + self._x[i + 1] * self._y[i]
            )
        I_xx = sxx / 12 - self._A * self._y_c ** 2
        I_yy = syy / 12 - self._A * self._x_c ** 2
        I_xy = sxy / 24 - self._A * self._x_c * self._y_c
        return I_xx, I_yy, I_xy

    def __principal(self) -> tuple[float, float, float]:
        avg = (self._I_xx + self._I_yy) / 2
        diff = (self._I_xx - self._I_yy) / 2
        k = sqrt(diff**2 + self._I_xy**2)
        I_max = avg + k
        I_min = avg - k
        theta = atan2(-self._I_xy, diff) / 2
        return I_max, I_min, theta


class HollowPolygon(Shape):

    def __init__(
        self,
        outer_polygon: Polygon,
        inner_polygon: Polygon
    ) -> None:
        super().__init__(Dimensions())
        self.outer_polygon = outer_polygon
        self.inner_polygon = inner_polygon
        # Shift inner polygon to center of outer polygon:
        x_c_o, y_c_o = self.outer_polygon.centroid
        self.inner_polygon.shift(x_c_o, y_c_o)

    @property
    def area(self) -> Quantity:
        A_o = self.outer_polygon.area
        A_i = self.inner_polygon.area
        A = A_o - A_i
        return A

    @property
    def centroid(self) -> tuple[Quantity, Quantity]:
        x_c_o, y_c_o = self.outer_polygon.centroid
        return x_c_o, y_c_o

    @property
    def moment_of_inertia_xx(self) -> Quantity:
        I_xx_o = self.outer_polygon.moment_of_inertia_xx
        I_xx_i = self.inner_polygon.moment_of_inertia_xx
        I_xx = I_xx_o - I_xx_i
        return I_xx

    @property
    def moment_of_inertia_yy(self) -> Quantity:
        I_yy_o = self.outer_polygon.moment_of_inertia_yy
        I_yy_i = self.inner_polygon.moment_of_inertia_yy
        I_yy = I_yy_o - I_yy_i
        return I_yy

    @property
    def product_of_inertia(self) -> Quantity:
        I_xy_o = self.outer_polygon.product_of_inertia
        I_xy_i = self.outer_polygon.product_of_inertia
        I_xy = I_xy_o - I_xy_i
        return I_xy

    @property
    def polar_moment_of_inertia(self) -> Quantity:
        I_xx = self.moment_of_inertia_xx
        I_yy = self.moment_of_inertia_yy
        J = I_xx + I_yy
        return J

    @property
    def principal_moments_of_inertia(self) -> tuple[Quantity, Quantity, Quantity]:
        I_xx = self.moment_of_inertia_xx
        u = I_xx.units
        I_xx = I_xx.magnitude
        I_yy = self.moment_of_inertia_yy.magnitude
        I_xy = self.product_of_inertia.magnitude
        avg = (I_xx + I_yy) / 2
        diff = (I_xx - I_yy) / 2
        k = sqrt(diff ** 2 + I_xy ** 2)
        I_max = avg + k
        I_min = avg - k
        theta = atan2(-I_xy, diff) / 2
        return Q_(I_max, u), Q_(I_min, u), Q_(theta, 'rad')

    # noinspection PyProtectedMember
    def plot(self) -> LineChart:
        x_c, y_c = self.centroid
        *_, theta = self.principal_moments_of_inertia
        theta = theta.to('rad').magnitude
        x_o = [x - x_c.m for x in self.outer_polygon._x]
        y_o = [y - y_c.m for y in self.outer_polygon._y]
        x_i = [x - x_c.m for x in self.inner_polygon._x]
        y_i = [y - y_c.m for y in self.inner_polygon._y]
        x_min = min(x_o)
        x_max = max(x_o)
        y_min = min(y_o)
        y_max = max(y_o)
        b = .05 * max(x_max - x_min, y_max - y_min)
        l = min(x_max - x_min, y_max - y_min) / 10
        x_a1 = [
            -l * cos(theta),
            l * cos(theta)
        ]
        y_a1 = [
            -l * sin(theta),
            l * sin(theta)
        ]
        x_a2 = [
            -l * cos(theta + pi / 2),
            l * cos(theta + pi / 2)
        ]
        y_a2 = [
            -l * sin(theta + pi / 2),
            l * sin(theta + pi / 2)
        ]
        chart = LineChart()
        chart.add_xy_data(
            label='outer_polygon',
            x1_values=x_o,
            y1_values=y_o,
            style_props={'color': 'black', 'lw': 2}
        )
        chart.add_xy_data(
            label='inner_polygon',
            x1_values=x_i,
            y1_values=y_i,
            style_props={'color': 'black', 'lw': 2}
        )
        chart.add_xy_data(
            label='major_principal_axis',
            x1_values=x_a1,
            y1_values=y_a1,
            style_props={'color': '#0072B2', 'lw': 2}
        )
        chart.add_xy_data(
            label='minor_principal_axis',
            x1_values=x_a2,
            y1_values=y_a2,
            style_props={'color': '#D55E00'}
        )
        chart.axes.set_aspect('equal')
        chart.axes.set_xlim(xmin=x_min - b, xmax=x_max + b)
        chart.axes.set_ylim(ymin=y_min - b, ymax=y_max + b)
        return chart

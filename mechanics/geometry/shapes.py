from typing import Type
from dataclasses import dataclass
from enum import Enum
from math import pi
import sympy as sp
from abc import ABC, abstractmethod
from mechanics import Quantity

Q_ = Quantity


class Shape(ABC):

    def __init__(self, *args, **kwargs) -> None:
        pass

    @abstractmethod
    def area(self) -> Quantity:
        pass

    @abstractmethod
    def second_moment_of_area(self) -> Quantity:
        pass

    @abstractmethod
    def product_moment_of_area(self) -> Quantity:
        pass


class Rectangle(Shape):

    def __init__(self, width: Quantity, height: Quantity) -> None:
        super().__init__()
        self.width = width.to('m')
        self.height = height.to('m')

    def area(self) -> Quantity:
        A = self.width * self.height
        return A.to('m ** 2')

    def second_moment_of_area(self) -> Quantity:
        I_xx = self.width * self.height ** 3 / 12
        return I_xx.to('m ** 4')

    def product_moment_of_area(self) -> Quantity:
        return Quantity(0.0, 'm ** 4')


class Triangle(Shape):

    def __init__(self, base: Quantity, height: Quantity, top_offset: Quantity) -> None:
        """Creates a `Triangle` shape.

        Parameters
        ----------
        base:
            Base width of the triangle (horizontal dimension).
        height:
            Height of the triangle (vertical dimension).
        top_offset:
            Horizontal distance between the top of the triangle and its lower
            right vertex.
        """
        super().__init__()
        self.base = base.to('m')
        self.height = height.to('m')
        self.top_offset = top_offset.to('m')

    def area(self) -> Quantity:
        A = self.base * self.height / 2
        return A.to('m ** 2')

    def second_moment_of_area(self) -> Quantity:
        I_xx = self.base * self.height ** 3 / 36
        return I_xx.to('m ** 4')

    def product_moment_of_area(self) -> Quantity:
        I_xy = self.base * (self.base - 2 * self.top_offset) * self.height ** 2 / 72
        return I_xy.to('m ** 4')

    def center_of_gravity(self) -> tuple[Quantity, Quantity]:
        """Returns the position of the center of gravity with respect to the
        lower left vertex of the triangle.
        """
        b = self.base.m
        h = self.height.m
        s = self.top_offset.m
        # Determine the expression of the median that connects the top and the
        # midpoint of the triangle's base. The general expression of this
        # straight line is `y = k * (x - a)`.
        x1, y1 = b - s, h
        x2, y2 = b / 2, 0
        k = (y2 - y1) / (x2 - x1)  # slope of the median
        a = sp.Symbol('a')
        eq = k * (b / 2 - a)
        sol = sp.solve(eq)
        a = float(sol[0])
        # The y-coordinate of the center of gravity lies at 1/3 of the
        # triangle's height with respect to the base, i.e. `y_cg = h/3` in our
        # reference frame
        x_cg, y_cg = sp.Symbol('x_cg'), h / 3
        eq = k * (x_cg - a) - y_cg
        sol = sp.solve(eq)
        x_cg = float(sol[0])
        x_cg = Q_(x_cg, 'm')
        y_cg = Q_(y_cg, 'm')
        return x_cg, y_cg


class Circle(Shape):

    def __init__(self, radius: Quantity) -> None:
        super().__init__()
        self.radius = radius

    def area(self) -> Quantity:
        A = pi * self.radius ** 2
        return A.to('m ** 2')

    def second_moment_of_area(self) -> Quantity:
        I_xx = pi * self.radius ** 4 / 4
        return I_xx.to('m ** 4')

    def product_moment_of_area(self) -> Quantity:
        return Quantity(0.0, 'm ** 4')


class SemiCircle(Shape):

    def __init__(self, radius: Quantity) -> None:
        super().__init__()
        self.radius = radius

    def area(self) -> Quantity:
        A = pi * self.radius ** 2 / 2
        return A

    def second_moment_of_area(self) -> Quantity:
        I_xx = 0.10976 * self.radius ** 4
        return I_xx.to('m ** 4')

    def product_moment_of_area(self) -> Quantity:
        return Quantity(0.0, 'm ** 4')


@dataclass
class Dimensions:
    width: Quantity | None = None
    height: Quantity | None = None
    base: Quantity | None = None
    top_offset: Quantity | None = None
    radius: Quantity | None = None


class ShapeEnum(Enum):
    RECTANGLE = Rectangle
    TRIANGLE = Triangle
    CIRCLE = Circle
    SEMICIRCLE = SemiCircle


def create_shape(shape: Type[Shape], dim: Dimensions) -> Shape:
    match shape:
        case ShapeEnum.RECTANGLE.value:
            return Rectangle(dim.width, dim.height)
        case ShapeEnum.TRIANGLE.value:
            return Triangle(dim.base, dim.height, dim.top_offset)
        case ShapeEnum.CIRCLE.value:
            return Circle(dim.radius)
        case ShapeEnum.SEMICIRCLE.value:
            return SemiCircle(dim.radius)
        case _:
            raise ValueError(f"'{shape.__name__}' is not defined")

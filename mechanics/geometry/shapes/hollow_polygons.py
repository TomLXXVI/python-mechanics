from mechanics import Quantity
from .shapes import Dimensions, HollowPolygon
from .solid_polygons import Rectangle, Circle


Q_ = Quantity


class HollowRectangle(HollowPolygon):
    """
    Represents a hollow rectangle shape.
    """
    def __init__(
        self,
        width: Quantity,
        height: Quantity,
        thickness: Quantity
    ) -> None:
        """Creates a `HollowRectangle` object.

        Parameters
        ----------
        width:
            Horizontal dimension of the outer rectangle.
        height:
            Vertical dimension of the outer rectangle.
        thickness:
            Space between outer and inner rectangle.
        """
        outer_rectangle = Rectangle(width, height)
        inner_rectangle = Rectangle(width - 2 * thickness, height - 2 * thickness)
        super().__init__(outer_rectangle, inner_rectangle)
        self.dim = Dimensions(width=width, height=height, thickness=thickness)


class Annulus(HollowPolygon):

    def __init__(
        self,
        outer_radius: Quantity,
        thickness: Quantity
    ) -> None:
        outer_circle = Circle(outer_radius)
        inner_circle = Circle(outer_radius - thickness)
        super().__init__(outer_circle, inner_circle)
        self.dim = Dimensions(outer_radius=outer_radius, thickness=thickness)

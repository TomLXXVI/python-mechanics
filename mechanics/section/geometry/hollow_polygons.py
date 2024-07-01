from mechanics import Quantity
from .shapes import HollowPolygon
from .polygons import Rectangle


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

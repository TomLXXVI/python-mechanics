# Demo based on example 7.1 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.

# Part 1: Derive a new polygon class `TShape` from base class `Polygon`.

from mechanics import Quantity
from mechanics.geometry.shapes import Polygon, Dimensions


Q_ = Quantity


class TShape(Polygon):

    def __init__(
        self,
        width: Quantity,
        height: Quantity,
        web_thickness: Quantity,
        flange_thickness: Quantity
    ) -> None:
        """Creates a `TShape` object."""
        super().__init__(
            dim=Dimensions(
                width=width,
                height=height,
                web_thickness=web_thickness,
                flange_thickness=flange_thickness
            ),
            vertices=self.__create_vertices(
                width,
                height,
                web_thickness,
                flange_thickness
            )
        )

    @staticmethod
    def __create_vertices(
        width: Quantity,
        height: Quantity,
        web_thickness: Quantity,
        flange_thickness: Quantity
    ) -> list[Quantity]:
        """Creates the vertices of the T-shape based on the main dimensions
        of the shape. The vertices must be ordered counterclockwise.
        """
        u = width.units
        w = width.m
        h = height.m
        t_w = web_thickness.m
        t_f = flange_thickness.m
        pts = [Q_([0, 0], u) for _ in range(8)]
        pts[1] = Q_([t_w, 0], u)
        pts[2] = Q_([t_w, h - t_f], u)
        pts[3] = Q_([(t_w + w) / 2, h - t_f], u)
        pts[4] = Q_([(t_w + w) / 2, h], u)
        pts[5] = Q_([(t_w - w) / 2, h], u)
        pts[6] = Q_([(t_w - w) / 2, h - t_f], u)
        pts[7] = Q_([0, h - t_f], u)
        return pts

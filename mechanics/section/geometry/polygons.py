from mechanics import Quantity
from mechanics.utils import AxesRotation2D
from .shapes import Polygon, Dimensions, TPoints

Q_ = Quantity


class Rectangle(Polygon):
    """
    Represents a rectangle shape.
    """
    def __init__(
        self,
        width: Quantity,
        height: Quantity
    ) -> None:
        """Creates a `Rectangle` object.

        Parameters
        ----------
        width:
            Horizontal dimension of the rectangle.
        height:
            Vertical dimension of the rectangle.
        """
        super().__init__(
            vertices=self.__create_vertices(
                width,
                height
            ),
            dim=Dimensions(
                width=width,
                height=height
            )
        )

    @staticmethod
    def __create_vertices(
        width: Quantity,
        height: Quantity
    ) -> TPoints:
        u = width.units
        w = width.m
        h = height.m
        p1 = Q_([0, 0], u)
        p2 = Q_([0, h], u)
        p3 = Q_([w, h], u)
        p4 = Q_([w, 0], u)
        return [p1, p4, p3, p2]


class HShape(Polygon):
    """
    Represents a wide-flange section or H-shaped section.
    """
    def __init__(
        self,
        height: Quantity,
        width: Quantity,
        web_thickness: Quantity,
        flange_thickness: Quantity
    ) -> None:
        """Creates a `WideFlange` object.

        Parameters
        ----------
        height:
            Height of the section.
        width:
            Width of the section.
        web_thickness:
            Thickness of the web.
        flange_thickness:
            Thickness of the flanges.
        """
        super().__init__(
            vertices=self.__create_vertices(
                height,
                width,
                web_thickness,
                flange_thickness
            ),
            dim=Dimensions(
                height=height,
                width=width,
                web_thickness=web_thickness,
                flange_thickness=flange_thickness
            )
        )

    @staticmethod
    def __create_vertices(
        height: Quantity,
        width: Quantity,
        web_thickness: Quantity,
        flange_thickness: Quantity
    ) -> TPoints:
        u = width.units
        h = height.m
        w = width.m
        t_w = web_thickness.m
        t_f = flange_thickness.m
        pts = [Q_([0, 0], u) for _ in range(12)]
        pts[1] = Q_([0, t_f], u)
        pts[2] = Q_([(w - t_w) / 2, t_f], u)
        pts[3] = Q_([(w - t_w) / 2, h - t_f], u)
        pts[4] = Q_([0, h - t_f], u)
        pts[5] = Q_([0, h], u)
        pts[6] = Q_([w, h], u)
        pts[7] = Q_([w, h - t_f], u)
        pts[8] = Q_([(w + t_w) / 2, h - t_f], u)
        pts[9] = Q_([(w + t_w) / 2, t_f], u)
        pts[10] = Q_([w, t_f], u)
        pts[11] = Q_([w, 0], u)
        pts.reverse()
        return pts


class CShape(Polygon):
    """
    Represents a channel section or C-shaped section.
    """
    def __init__(
        self,
        height: Quantity,
        width: Quantity,
        web_thickness: Quantity,
        flange_thickness: Quantity,
        orientation: Quantity = Q_(0, 'deg')
    ) -> None:
        """Creates a `Channel` object.

        Parameters
        ----------
        height:
            Height of the section.
        width:
            Width of the section.
        web_thickness:
            Thickness of the web.
        flange_thickness:
            Thickness of the flanges.
        orientation:
            Rotation angle of the section.

        Notes
        -----
        The dimensions of the C-shaped section are valid when the web is
        vertical. The section can be tilted by adjusting the orientation angle
        of the section.
        """
        super().__init__(
            vertices=self.__create_vertices(
                height,
                width,
                web_thickness,
                flange_thickness,
                orientation
            ),
            dim=Dimensions(
                height=height,
                width=width,
                web_thickness=web_thickness,
                flange_thickness=flange_thickness
            )
        )

    @staticmethod
    def __create_vertices(
        height: Quantity,
        width: Quantity,
        web_thickness: Quantity,
        flange_thickness: Quantity,
        orientation: Quantity
    ) -> TPoints:
        u = width.units
        h = height.m
        w = width.m
        t_w = web_thickness.m
        t_f = flange_thickness.m
        pts = [Q_([0, 0], u) for _ in range(8)]
        pts[1] = Q_([0, h], u)
        pts[2] = Q_([w, h], u)
        pts[3] = Q_([w, h - t_f], u)
        pts[4] = Q_([t_w, h - t_f], u)
        pts[5] = Q_([t_w, t_f], u)
        pts[6] = Q_([w, t_f], u)
        pts[7] = Q_([w, 0], u)
        pts.reverse()
        if orientation.m != 0.0:
            rot = AxesRotation2D(-orientation)
            pts = [rot(p) for p in pts]
        return pts


class ZShape(Polygon):

    def __init__(
        self,
        width: Quantity,
        height: Quantity,
        web_thickness: Quantity,
        flange_thickness: Quantity,
        orientation: Quantity = Q_(0, 'deg')
    ) -> None:
        super().__init__(
            vertices=self.__create_vertices(
                width,
                height,
                web_thickness,
                flange_thickness,
                orientation
            ),
            dim=Dimensions(
                width=width,
                height=height,
                web_thickness=web_thickness,
                flange_thickness=flange_thickness
            )
        )

    @staticmethod
    def __create_vertices(
        width: Quantity,
        height: Quantity,
        web_thickness: Quantity,
        flange_thickness: Quantity,
        orientation: Quantity
    ) -> TPoints:
        u = width.units
        w = width.m
        h = height.m
        t_w = web_thickness.m
        t_f = flange_thickness.m
        pts = [Q_([0, 0], u) for _ in range(8)]
        pts[1] = Q_([w - t_f, 0], u)
        pts[2] = Q_([w - t_f, -(h - t_w)], u)
        pts[3] = Q_([w, -(h - t_w)], u)
        pts[4] = Q_([w, t_w], u)
        pts[5] = Q_([t_f, t_w], u)
        pts[6] = Q_([t_f, h], u)
        pts[7] = Q_([0, h], u)
        if orientation.m != 0.0:
            rot = AxesRotation2D(-orientation)
            pts = [rot(p) for p in pts]
        return pts

# Demo based on example 7.4 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.

from mechanics import Quantity
from mechanics.strength import TransverseShear
from mechanics.geometry import Polygon

Q_ = Quantity


class PiShape(Polygon):

    def __init__(
        self,
        top_width: Quantity,
        bottom_width: Quantity,     # including the thickness of the legs
        height: Quantity,           # including the thickness of the top board
        web_thickness: Quantity,    # thickness of the top board
        flange_thickness: Quantity  # thickness of a single leg
    ) -> None:
        vertices = self.__create_vertices(
            top_width,
            bottom_width,
            height,
            web_thickness,
            flange_thickness
        )
        super().__init__(vertices)

    @staticmethod
    def __create_vertices(
        top_width: Quantity,
        bottom_width: Quantity,
        height: Quantity,
        web_thickness: Quantity,
        flange_thickness: Quantity
    ) -> list[Quantity]:
        u = top_width.units
        w_top = top_width.m
        w_bot = bottom_width.m
        h = height.m
        t_w = web_thickness.m
        t_f = flange_thickness.m
        w_w = (w_top - w_bot) / 2
        # points in counter-clockwise order:
        pts = [Q_([0, 0], u) for _ in range(12)]
        pts[1] = Q_([0, -t_f], u)
        pts[2] = Q_([w_w, -t_f], u)
        pts[3] = Q_([w_w, -h], u)
        pts[4] = Q_([w_w + t_w, -h], u)
        pts[5] = Q_([w_w + t_w, -t_f], u)
        pts[6] = Q_([w_w + w_bot - t_w, -t_f], u)
        pts[7] = Q_([w_w + w_bot - t_w, -h], u)
        pts[8] = Q_([w_w + w_bot, -h], u)
        pts[9] = Q_([w_w + w_bot, -t_f], u)
        pts[10] = Q_([w_top, -t_f], u)
        pts[11] = Q_([w_top, 0], u)
        return pts


if __name__ == '__main__':

    section = PiShape(
        top_width=Q_(250, 'mm'),
        bottom_width=Q_(145, 'mm'),
        height=Q_(310, 'mm'),
        web_thickness=Q_(10, 'mm'),
        flange_thickness=Q_(10, 'mm')
    )
    section.plot().show()

    for vertex in section.vertices: print(vertex)

    y_B = section.vertices[1][1]

    shear = TransverseShear(Q_(850, 'kN'), section)
    q = shear.flow(y=y_B)
    print(q.to('MN / m'))

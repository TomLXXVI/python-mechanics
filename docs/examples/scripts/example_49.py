# Demo based on example 12.10 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.
import pandas as pd
from mechanics import Quantity
from mechanics.statics import Position, Angle, Force
from mechanics.strength import Beam, Hinge, Roller
from mechanics.geometry import Rectangle, Dimensions
from mechanics.charts import LineChart


Q_ = Quantity


def create_beam():
    R_A = Force(
        magnitude='R_A',
        theta=Angle(90),
        position=Position(0, units='m'),
        units='kN'
    )
    R_B = Force(
        magnitude='R_B',
        theta=Angle(90),
        position=Position(8, units='m'),
        units='kN'
    )
    F = Force(
        magnitude=16,
        theta=Angle(-90),
        position=Position(6, units='m'),
        units='kN',
        name='F'
    )
    beam = Beam(
        length=Q_(8, 'm'),
        shape_type=Rectangle,
        shape_dim=Dimensions(width=Q_(100, 'mm'), height=Q_(300, 'mm')),
        E_modulus=Q_(200, 'GPa'),
        G_modulus=Q_(75, 'GPa'),
        loadings=[R_A, R_B, F],
        supports=[Hinge(x=Q_(0, 'm')), Roller(x=Q_(8, 'm'))],
        units=('kN', 'm'),
        num_sections=100
    )
    return beam


def main():
    beam = create_beam()
    beam.moment_diagram.show()

    data = {
        'x': beam.elastic_curve.x.to('m').magnitude,
        'y': beam.elastic_curve.y.to('mm').magnitude,
        'theta': beam.elastic_curve.theta.to('rad').magnitude
    }
    df = pd.DataFrame(data)
    with pd.option_context('display.max_rows', None):
        print(df)

    beam.elastic_curve.diagram().show()

    diagram = LineChart()
    diagram.add_xy_data(
        label='elastic curve - slope',
        x1_values=beam.elastic_curve.x.to('m').magnitude,
        y1_values=beam.elastic_curve.theta.to('rad').magnitude,
        style_props={'drawstyle': 'steps-post'}
    )
    diagram.x1.add_title('x, m')
    diagram.y1.add_title('theta, rad')
    diagram.show()


if __name__ == '__main__':
    main()

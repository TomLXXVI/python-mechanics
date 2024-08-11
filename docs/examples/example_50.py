# Demo based on example 12.12 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.
from mechanics import Quantity
from mechanics.statics import Position, Angle, Force
from mechanics.strength import Beam, Hinge, Roller, FreeEnd
from mechanics.geometry import HShape, Dimensions
from mechanics.charts import LineChart


Q_ = Quantity


def create_beam():
    R_A = Force(
        magnitude='R_A',
        theta='theta_A',
        position=Position(0, units='m'),
        units='kN'
    )
    R_B = Force(
        magnitude='R_B',
        theta=Angle(90),
        position=Position(4, units='m'),
        units='kN'
    )
    F_C = Force(
        magnitude=25,
        theta=Angle(-90),
        position=Position(8, units='m'),
        units='kN',
        name='F_C'
    )
    beam = Beam(
        length=Q_(8, 'm'),
        shape_type=HShape,
        shape_dim=Dimensions(
            width=Q_(102, 'mm'),
            height=Q_(313, 'mm'),
            web_thickness=Q_(6.60, 'mm'),
            flange_thickness=Q_(10.8, 'mm')
        ),
        E_modulus=Q_(200, 'GPa'),
        G_modulus=Q_(75, 'GPa'),
        loadings=[R_A, R_B, F_C],
        supports=[
            Hinge(x=Q_(0, 'm')),
            Roller(x=Q_(4, 'm')),
            FreeEnd(x=Q_(8, 'm'))
        ],
        units=('kN', 'm'),
        num_sections=500
    )
    return beam


def main():
    beam = create_beam()

    beam.moment_diagram.show()

    beam.elastic_curve.diagram().show()

    diagram = LineChart()
    diagram.add_xy_data(
        label='elastic curve - slope',
        x1_values=beam.elastic_curve.x.to('m').magnitude,
        y1_values=beam.elastic_curve.theta.to('deg').magnitude,
        style_props={'drawstyle': 'steps-post'}
    )
    diagram.x1.add_title('position, m')
    diagram.y1.add_title('slope, deg')
    diagram.show()

    print(f"{beam.elastic_curve.displacement(x=Q_(8, 'm')).to('mm'):~P.0f}")


if __name__ == '__main__':
    main()

# Demo based on example 5.5 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.
import time
from mechanics import Quantity
from mechanics.statics import Position, Angle, Moment
from mechanics.strength import Beam
from mechanics.geometry import Circle, Dimensions


Q_ = Quantity

if __name__ == '__main__':
    # on Windows, the script needs the line above since parallel processing is
    # used when parameter `num_sections` in the constructor of class `Beam` is
    # large (> 350).
    start = time.perf_counter()

    T_A = Moment(
        magnitude=80,
        theta=Angle(180),
        position=Position(0, units='m'),
        units='kN * m',
        name='T_A'
    )

    T_B = Moment(
        magnitude=150,
        theta=Angle(0),
        position=Position(3, units='m'),
        units='kN * m',
        name='T_B'
    )

    T_C = Moment(
        magnitude=60,
        theta=Angle(180),
        position=Position(5, units='m'),
        units='kN * m',
        name='T_C'
    )

    T_D = Moment(
        magnitude=10,
        theta=Angle(180),
        position=Position(6.5, units='m'),
        units='kN * m',
        name='T_D'
    )

    beam = Beam(
        length=Q_(6.5, 'm'),
        shape_type=Circle,
        shape_dim=Dimensions(diameter=Q_(200, 'mm')),
        E_modulus=Q_(200, 'GPa'),
        G_modulus=Q_(75, 'GPa'),
        loadings=[T_A, T_B, T_C, T_D],
        units=('kN', 'm'),
        num_sections=650
    )

    beam.torque_diagram.show()

    fi_1 = beam.angle_of_twist(x1=Q_(0, 'm'), x2=Q_(3, 'm'))
    fi_2 = beam.angle_of_twist(x1=Q_(3, 'm'), x2=Q_(5, 'm'))
    fi_3 = beam.angle_of_twist(x1=Q_(5, 'm'), x2=Q_(6.5, 'm'))
    fi_A = fi_1 + fi_2 + fi_3
    print(fi_A.to('rad'))
    fi_AC = fi_1 + fi_2
    print(fi_AC.to('rad'))

    end = time.perf_counter()
    print(f"execution time {end - start} s")

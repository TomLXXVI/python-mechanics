# Demo based on example 7.5 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.

from mechanics import Quantity
from mechanics.statics import Position, Angle, Force, Moment, Beam
from mechanics.strength import Section
from mechanics.geometry import HollowRectangle, Dimensions


Q_ = Quantity


def create_beam():
    F = Force(
        position=Position(1, units='m'),
        magnitude=80,
        theta=Angle(-90, 'deg'),
        units='N',
        name='F'
    )
    R_w = Force(
        magnitude='R_w',
        theta='theta_w',
        units='N',
        position=Position(0, units='m')
    )
    M_w = Moment(
        magnitude='M_w',
        gamma=Angle(90, 'deg'),
        units='N * m',
        position=Position(0, units='m')
    )
    beam = Beam(
        length=Q_(1, 'm'),
        loadings=[F, R_w, M_w],
    )
    beam.shear_diagram.show()
    beam.moment_diagram.show()
    return beam


def create_section(F_i, M_i):
    section = Section(
        shape=HollowRectangle,
        dim=Dimensions(
            width=Q_(75, 'mm'),
            height=Q_(75, 'mm'),
            thickness=Q_(15, 'mm')
        )
    )
    section.set_internal_loadings(F_i, M_i)
    return section


def main():
    beam = create_beam()

    x1 = Q_(0.50, 'm')
    delta_x = Q_(0.001, 'mm')
    x2 = x1 + delta_x

    F_i1, M_i1 = beam.cut(x=x1, view='right')
    F_i2, M_i2 = beam.cut(x=x2, view='left')

    section1 = create_section(F_i1, M_i1)
    section1.shape.plot().show()
    section2 = create_section(F_i2, M_i2)

    sigma_1 = section1.bending.sigma(z=Q_(0, 'mm'), y=Q_(30, 'mm')).to('MPa')
    sigma_2 = section2.bending.sigma(z=Q_(0, 'mm'), y=Q_(30, 'mm')).to('MPa')

    print('section 1 (right-side): ', M_i1, sigma_1)
    print('section 2 (left-side): ', M_i2, sigma_2)

    delta_sigma = abs(sigma_1) - abs(sigma_2)
    t = Q_(15, 'mm')
    w = Q_(45, 'mm')
    A = t * w
    delta_F = delta_sigma * A

    q = delta_F / delta_x
    print('q: ', q.to('N / m'))

    tau_avg = q / 2 / t
    print('tau_avg: ', tau_avg.to('MPa'))


if __name__ == '__main__':
    main()

# Demo based on example 8.2 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.

from mechanics import Quantity
from mechanics.statics import Position, Angle, Force, Moment, Beam
from mechanics.strength import Section
from mechanics.geometry import Rectangle, Dimensions


Q_ = Quantity


def create_beam():
    R_w = Force(
        magnitude='R_w',
        theta='theta_w',
        position=Position(0, 0, 0, units='m'),
        units='kN'
    )
    M_w = Moment(
        magnitude='M_w',
        gamma=Angle(90),
        position=Position(0, 0, 0, units='m'),
        units='kN * m'
    )
    F = Force(
        magnitude=300,
        theta=Angle(180),
        position=Position(1.000, -0.150, 0, units='m'),
        units='kN',
        name='F'
    )
    beam = Beam(length=Q_(1, 'm'), loadings=[R_w, M_w, F], units=('kN', 'm'))
    return beam


def create_section():
    section = Section(
        shape=Rectangle,
        dim=Dimensions(
            width=Q_(100, 'mm'),
            height=Q_(300, 'mm')
        )
    )
    return section


def main():
    beam = create_beam()
    section = create_section()

    F_i, M_i = beam.cut(x=Q_(0.5, 'm'), view='right')

    print(f"resultant internal force: {F_i}")
    print(f"resultant internal moment: {M_i}")

    section.set_internal_loadings(F_i, M_i)

    p_B = Q_([-50, 150], 'mm')
    p_C = Q_([50, -150], 'mm')

    sigma_N = section.axial.sigma.to('MPa')
    sigma_M_B = section.bending.sigma(*p_B).to('MPa')
    sigma_M_C = section.bending.sigma(*p_C).to('MPa')

    print(f"axial stress at B and C: {sigma_N:~P.3f}")
    print(f"bending stress at B: {sigma_M_B:~P.3f}")
    print(f"bending stress at C: {sigma_M_C:~P.3f}")

    sigma_B = sigma_N + sigma_M_B
    sigma_C = sigma_N + sigma_M_C

    print(f"resultant normal stress at B: {sigma_B:~P.3f}")
    print(f"resultant normal stress at C: {sigma_C:~P.3f}")


if __name__ == '__main__':
    main()

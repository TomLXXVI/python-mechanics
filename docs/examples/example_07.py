# Demo based on example 1.5 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.
from mechanics import Quantity
from mechanics.statics import Position, Angle, Force, Beam
from mechanics.strength import Section
from mechanics.geometry import Rectangle, Dimensions


Q_ = Quantity


F_A = Force(
    magnitude=12,
    theta=Angle(180),
    position=Position(0, units='cm'),
    units='kN',
    name='F_A'
)

F_B = Force(
    magnitude=2*9,
    theta=Angle(180),
    position=Position(30, units='cm'),
    units='kN',
    name='F_B'
)

F_C = Force(
    magnitude=2*4,
    theta=Angle(0),
    position=Position(60, units='cm'),
    units='kN',
    name='F_C'
)

F_D = Force(
    magnitude=22,
    theta=Angle(0),
    position=Position(90, units='cm'),
    units='kN',
    name='F_D'
)

beam = Beam(
    length=Q_(90, 'cm'),
    loadings=[F_A, F_B, F_C, F_D],
    units=('kN', 'cm'),
    num_sections=100
)

beam.normal_force_diagram.show()

x, N_max = beam.N_max()
print(x, N_max)


S = Section(Rectangle, Dimensions(width=Q_(10, 'mm'), height=Q_(35, 'mm')))
S.set_normal_force(N_max)
print(S.sigma_N_avg.to('MPa'))

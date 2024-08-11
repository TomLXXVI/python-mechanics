# Demo based on example 4.3 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.
from mechanics import Quantity
from mechanics.statics import Position, Angle, Force, System
from mechanics.strength import Beam
from mechanics.geometry import Circle, Dimensions


Q_ = Quantity


# Reaction force exerted by vertical beam (column) AC on horizontal beam:
F_A = Force(
    magnitude='F_A',
    theta=Angle(90),
    position=Position(0, units='mm'),
    units='kN'
)

# Reaction force exerted by vertical beam (column) BD on horizontal beam:
F_B = Force(
    magnitude='F_B',
    theta=Angle(90),
    position=Position(600, units='mm'),
    units='kN'
)

F_l = Force(
    magnitude=90,
    theta=Angle(-90),
    position=Position(200, units='mm'),
    units='kN',
    name='F_l'
)

# Determine the reaction forces acting on horizontal beam at A and B:
horizontal_beam = System(loadings=[F_A, F_B, F_l], units=('kN', 'm'))
sol_dict = horizontal_beam.solve()
for name, force in sol_dict.items():
    print(name, force)


# Reaction force exerted by horizontal beam on vertical beam (column) AC:
# The reaction force exerted by vertical beam AC on horizontal beam points
# vertically upward. To get the reaction force on vertical beam (column) AC
# exerted by the horizontal beam, we would normally rotate (or reverse) this
# reaction force. However, in the attached coordinate system of a beam, the
# beam's longitudinal axis coincides with the x-axis and the direction of
# normal forces N. Therefore, we need to rotate the reaction force by 90° here.
F_A = Force.rotate(sol_dict['F_A'], theta=Angle(90), name='F_A')
F_A.position = Position(300, units='mm')
print(F_A.name, F_A, F_A.theta)


# Reaction force exerted by floor on vertical beam (column) AC:
F_C = Force.reverse(F_A, name='F_C')
F_C.position = Position(0, units='mm')
print(F_C.name, F_C, F_C.theta)


beam_AC = Beam(
    length=Q_(300, 'mm'),
    shape_type=Circle,
    shape_dim=Dimensions(diameter=Q_(20, 'mm')),
    E_modulus=Q_(200, 'GPa'),
    G_modulus=Q_(75, 'GPa'),
    loadings=[F_A, F_C],
    units=('kN', 'mm'),
    num_sections=50
)

delta_AC = beam_AC.elongation(x1=Q_(0, 'mm'), x2=Q_(300, 'mm'))
print(delta_AC.to('mm'))


# Reaction force exerted by horizontal beam on vertical beam (column) BD:
F_B = Force.rotate(sol_dict['F_B'], theta=Angle(90), name='F_B')
F_B.position = Position(300, units='mm')
print(F_B.name, F_B, F_B.theta)


# Reaction force exerted by floor on vertical beam (column) BD:
F_D = Force.reverse(F_B, name='F_D')
F_D.position = Position(0, units='mm')
print(F_D.name, F_D, F_D.theta)


beam_BD = Beam(
    length=Q_(300, 'mm'),
    shape_type=Circle,
    shape_dim=Dimensions(diameter=Q_(40, 'mm')),
    E_modulus=Q_(70, 'GPa'),
    G_modulus=Q_(27, 'GPa'),
    loadings=[F_B, F_D],
    units=('kN', 'mm'),
    num_sections=50
)

delta_BD = beam_BD.elongation(x1=Q_(0, 'mm'), x2=Q_(300, 'mm'))
print(delta_BD.to('mm'))

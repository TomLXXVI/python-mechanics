# Demo based on example 8.5 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.
from mechanics import Quantity
from mechanics.statics import (
    Position,
    Angle,
    Force,
    Moment,
    System
)
from mechanics.strength import Section
from mechanics.geometry import Circle, Dimensions


Q_ = Quantity


# Coordinate system
# See docs/notes/right_handed_coordinate_system.pdf
# All forces, moments, and positions (both external and internal) are referred
# to this one coordinate system.


# Determine resultant internal force and moment in cross-section through point A
F_A = Force(
    magnitude='F_A',
    theta='theta_A',
    position=Position(0, 0, 0, units='mm'),
    units='kN'
)
M_A = Moment(
    magnitude='M_A',
    theta=Angle(90),
    position=Position(0, 0, 0, units='mm'),
    units='kN * m'
)

F_B = Force(
    magnitude=2,
    theta=Angle(0),
    position=Position(150, 0, 200, units='mm'),
    units='kN'
)

rod = System(
    loadings=[F_A, M_A, F_B],
    units=('kN', 'm')
)
rod.solve()

F_A = rod.external_forces['F_A']
M_A = rod.external_moments['M_A']
print(F_A)
print(M_A)


# Determine the state of stress at point A
S_A = Section(
    shape=Circle,
    dim=Dimensions(radius=Q_(20, 'mm'))
)
S_A.set_internal_loadings(F_A, M_A)

sigma_N = S_A.axial.sigma
sigma_M = S_A.bending.sigma(z=Q_(20, 'mm'), y=Q_(0, 'mm'))
print(sigma_N.to('MPa'))
print(sigma_M.to('MPa'))

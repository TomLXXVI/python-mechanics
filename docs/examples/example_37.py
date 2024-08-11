# Demo based on example 8.6 from Hibbeler, R. C. (2017). Mechanics of Materials
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


# Determine resultant internal force and moment in cross-section at point A

F_A = Force(
    magnitude='F_A',
    theta=Angle(90),
    position=Position(0, 0, 0, units='mm'),
    units='kN'
)
M_A = Moment(
    magnitude='M_A',
    theta=Angle(0),
    gamma='gamma_A',
    position=Position(0, 0, 0, units='mm'),
    units='kN * m'
)

F_B = Force(
    magnitude=3,
    theta=Angle(-90),
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
print(f"resultant internal force: {F_A}")
print(f"resultant internal moment: {M_A}")


# Determine the state of stress at point A

S_A = Section(
    shape=Circle,
    dim=Dimensions(radius=Q_(20, 'mm'))
)
S_A.set_internal_loadings(F_A, M_A)

if S_A.axial:
    sigma_N = S_A.axial.sigma
    print(f"axial stress: {sigma_N.to('MPa')}")

if S_A.bending:
    sigma_M = S_A.bending.sigma(z=Q_(20, 'mm'), y=Q_(0, 'mm'))
    print(f"bending stress: {sigma_M.to('MPa'):~P.2f}")

if S_A.transverse_shear:
    tau_V = S_A.transverse_shear.tau(y=Q_(0, 'mm'))
    print(f"shear stress: {tau_V.to('MPa'):~P.2f}")

if S_A.torsion:
    tau_T = S_A.torsion.tau_max
    print(f"torsion stress: {tau_T.to('MPa'):~P.2f}")

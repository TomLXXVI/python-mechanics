# Demo based on example 1.2 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.
from mechanics import Quantity
from mechanics.statics import (
    Position, Angle,
    Force,
    Beam
)

Q_ = Quantity


F_C = Force(
    magnitude='F_C',
    theta=Angle.create(3, 4, quadrant=2),
    position=Position(2, units='m'),
)


F_A = Force(
    magnitude='F_A',
    theta='theta_A',
    position=Position(0, units='m')
)


F_B = Force(
    magnitude=9.81 * 500,
    theta=Angle(-90),
    position=Position(3, units='m'),
    name='F_B'
)


boom = Beam(
    length=Q_(3, 'm'),
    loadings=[F_C, F_A, F_B]
)

for name, force in boom.external_forces.items():
    print(name, force)

for name, torque in boom.external_moments.items():
    print(name, torque)

iF_E, iM_E = boom.cut(x=Q_(1, 'm'), view='left')
print(iF_E)
print(iM_E)

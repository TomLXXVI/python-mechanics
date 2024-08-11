# Demo based on example 1.12 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.
import numpy as np
from mechanics import Quantity
from mechanics.statics import (
    Position, Angle, Force,
    System
)
from mechanics.strength.stress import SimpleShear

Q_ = Quantity


R_AB = Force(
    position=Position(0, 200, units='mm'),
    magnitude='R_AB',
    theta=Angle(0, 'deg')
)


R_C = Force(
    position=Position(0, 0, units='mm'),
    magnitude='R_C',
    theta='theta_C'
)

F_D = Force(
    position=Position(75, 0, units='mm'),
    magnitude=15,
    theta=Angle(-90),
    units='kN',
    name='F_D'
)

F_E = Force(
    position=Position(125, 0, units='mm'),
    magnitude=25,
    theta=Angle.create(3, 4, quadrant=4),
    units='kN',
    name='F_E'
)

sys = System([R_AB, R_C, F_D, F_E])
sol_dict = sys.solve()
for name, force in sol_dict.items(): print(f"{name}: {force}")


R_AB = sol_dict['R_AB']
R_C = sol_dict['R_C']


# Determine the required cross-section area of the steel pin at A
A_AB = SimpleShear.design(
    V=R_AB.magnitude,
    tau_fail=Q_(82.5, 'MPa'),  # failure shear stress for the steel
    safety_factor=1.5   # factor of safety for shear
)
d_AB = np.sqrt(4 * A_AB / np.pi)
print(d_AB.to('mm'))


# Determine the required cross-section area of the steel pin at C
A_C = SimpleShear.design(
    V=R_C.magnitude / 2,  # double shear
    tau_fail=Q_(82.5, 'MPa'),
    safety_factor=1.5
)
d_C = np.sqrt(4 * A_C / np.pi)
print(d_C.to('mm'))

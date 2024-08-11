# Demo based on example 6.2 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.
from mechanics import Quantity
from mechanics.statics import (
    Position, Angle, Force, Moment,
    DistributedLoad1D,
    Beam
)

Q_ = Quantity


R_A = Force(
    magnitude='R_A',
    theta='theta_A',
    position=Position(0, units='m'),
)

M_A = Moment(
    magnitude='M_A',
    gamma=Angle(90),
    position=Position(0, units='m')
)

q = DistributedLoad1D(
    x_coords=Q_([0, 3], 'm'),
    loads=Q_([0, -2], 'kN / m')
)

beam = Beam(
    length=Q_(3, 'm'),
    loadings=[R_A, M_A, q]
)
print(beam.external_forces['R_A'])
print(beam.external_moments['M_A'])

beam.shear_diagram.show()
beam.moment_diagram.show()

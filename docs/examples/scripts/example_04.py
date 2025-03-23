# Demo based on example 6.1 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.
from mechanics import Quantity
from mechanics.statics import (
    Position, Angle, Force,
    DistributedLoad1D,
    Beam
)

Q_ = Quantity


R_A = Force(
    magnitude='R_A',
    theta=Angle(90),
    position=Position(0, units='m')
)

R_B = Force(
    magnitude='R_B',
    theta=Angle(90),
    position=Position(4, units='m')
)

q = DistributedLoad1D(
    x_coords=Q_([0, 4], 'm'),
    loads=Q_([-3, -3], 'kN / m')
)

beam = Beam(
    length=Q_(4, 'm'),
    loadings=[R_A, R_B, q]
)

beam.shear_diagram.show()
beam.moment_diagram.show()

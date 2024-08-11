# Demo based on example 6.3 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.
from mechanics import Quantity
from mechanics.statics import (
    Position, Angle, Force,
    DistributedLoad1D,
    Beam
)

Q_ = Quantity

# roller (1 unknown)
R_A = Force(
    magnitude='R_A',
    theta=Angle(90),
    position=Position(0, units='m')
)

# hinge (2 unknowns)
R_B = Force(
    magnitude='R_B',
    theta='theta_B',
    position=Position(6, units='m')
)

q = DistributedLoad1D(
    x_coords=Q_([0, 6], 'm'),
    loads=Q_([10, 30], 'kN / m')
)


beam = Beam(
    length=Q_(6, 'm'),
    loadings=[R_A, R_B, q],
    units=('kN', 'm'),
    num_sections=100
)


beam.shear_diagram.show()
beam.moment_diagram.show()


x, V_max = beam.V_max()
print(x, V_max)


x, M_max = beam.M_max()
print(x, M_max)

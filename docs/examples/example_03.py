# Demo based on example 1.3 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.
from mechanics import Quantity
from mechanics.statics import (
    Position, Angle,
    Force,
    DistributedLoad1D,
    System,
    Beam
)

Q_ = Quantity


F_A = Force(
    magnitude=1500,
    theta=Angle(-90),
    position=Position(0, 0, units='m'),
    units='N',
    name='F_A'
)

q = DistributedLoad1D(
    x_coords=Q_([2, 5], 'm'),
    loads=Q_([-600, 0], 'N / m'),
    name='q'
)

F_C = Force(
    magnitude='F_C',
    theta=Angle(0),
    position=Position(5, 1.5, units='m')
)

F_E = Force(
    magnitude='F_E',
    theta='theta_E',
    position=Position(5.0, 0.0, units='m')
)


sys = System([F_A, F_C, F_E, q])

solutions = sys.solve()

for name, solution in solutions.items():
    print(f"{name}: {solution}")

# ------------------------------------------------------------------------------

R_A = Force(
    magnitude='R_A',
    theta=Angle.create(3, 4, quadrant=1),
    position=Position(0, units='m')
)

R_D = Force(
    magnitude='R_D',
    position=Position(2, units='m'),
    theta=Angle(90)
)

F_E = solutions['F_E']

beam = Beam(
    length=Q_(5, 'm'),
    loadings=[F_A, q, F_E, R_A, R_D]
)

iF_G, iM_G = beam.cut(x=Q_(1, 'm'), view='left')

print(iF_G)
print(iM_G)

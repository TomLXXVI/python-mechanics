# Demo: Determine the resultant internal force and moment in a cross-section
# of an externally loaded shaft.
from mechanics import Quantity
from mechanics.statics import (
    Position, Angle,
    Force,
    DistributedLoad1D,
    Beam
)

Q_ = Quantity


# Coordinate system
# See docs/notes/right_handed_coordinate_system.pdf
# All forces, moments, and positions (both external and internal) are referred
# to this one coordinate system.


# External loadings acting on the shaft:

F_D = Force(
    magnitude=225,
    theta=Angle(-90, 'deg'),
    position=Position(500, units='mm'),
)

q = DistributedLoad1D(
    x_coords=Q_([200, 350], 'mm'),
    loads=Q_([-800, -800], 'N / m')
)

# Unknown reaction force exerted on the shaft by a roller support:
R_A = Force(
    magnitude='R_A',
    theta=Angle(90, 'deg'),
    position=Position(0, units='mm'),
)

# Unknown reaction force exerted on the shaft by a roller support:
R_B = Force(
    magnitude='R_B',
    theta=Angle(90),
    position=Position(400, units='mm')
)

# Define the shaft (any unknown reaction forces will be calculated on
# instantiation)
shaft = Beam(
    length=Q_(500, 'mm'),
    loadings=[F_D, q, R_A, R_B]
)

# Get the internal force and moment in the cross-section at x = 250 mm.
iF_C, iM_C = shaft.cut(x=Q_(250, 'mm'), view='left')

print(iF_C)
print(iM_C)

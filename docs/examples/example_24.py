# Demo based on example 6.15 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.

# The rectangular cross section shown in [Fig. 6–33a] is subjected to a bending
# moment of M = 12 kN*m. Determine the normal stress developed at each corner of
# the section, and specify the orientation of the neutral axis.

from mechanics import Quantity
from mechanics.statics import Angle, Moment
from mechanics.strength import Bending
from mechanics.geometry import Rectangle

Q_ = Quantity


# Define the resultant internal bending moment at the cross-section (see
# example_24.pdf in docs/notes).
M = Moment(
    magnitude=12,
    theta=Angle(-90, 'deg'),
    gamma=Angle.create(4, 3, quadrant=3),
    units='kN * m'
)
print(M)


# Define the shape of the cross-section
rect = Rectangle(width=Q_(0.4, 'm'), height=Q_(0.2, 'm'))


# Apply the bending moment to the cross-section:
bending = Bending(
    M=(M.z, M.y),
    shape=rect
)

# Get the maximum bending stress in the cross-section:
print(f"sigma_max = {bending.sigma_max.to('MPa'):~P.2f}")


# Get the bending stresses at the corners of the rectangle:
for p in rect.vertices:
    print(
        f"point ({p[0]:~P.1f}, {p[1]:~P.1f}): "
        f"sigma = {bending.sigma(p[0], p[1]).to('MPa'):~P.2f}"
    )

# Demo based on example 6.13 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.

# The beam shown in [Fig. 6–27a] has a cross-sectional area in the shape of a
# channel, [Fig. 6–27b]. Determine the maximum bending stress that occurs in the
# beam at section a–a.

from mechanics import Quantity
from mechanics.statics import Position, Angle, Force, Moment, Beam
from mechanics.strength import Section
from mechanics.geometry import CShape, Dimensions


Q_ = Quantity

# ------------------------------------------------------------------------------
# BEAM WITH EXTERNAL LOADINGS:

F = Force(
    position=Position(0, 0.05909, units='m'),  # y = vertical distance from centroid of C-shape
    magnitude=2.6,
    theta=Angle.create(12, 5, quadrant=3),
    units='kN',
    name='F'
)

# fixed end: reaction force
R_A = Force(
    position=Position(3, units='m'),
    magnitude='R_A',
    theta='theta_A',
    units='kN'
)

# fixed end: reaction moment about z-axis:
M_A = Moment(
    position=Position(3, units='m'),
    magnitude='M_A',
    gamma=Angle(90),
    units='kN * m'
)


beam = Beam(
    length=Q_(3, 'm'),
    loadings=[F, R_A, M_A],
    units=('kN', 'm')
)

# M-diagram
beam.moment_diagram.show()

# Resultant internal force and moment at section a-a:
F_i, M_i = beam.cut(x=Q_(2, 'm'))

# ------------------------------------------------------------------------------
# CROSS-SECTION

# Shape and dimensions of cross-section:
dim = Dimensions(
    height=Q_(280, 'mm'),
    width=Q_(200, 'mm'),
    web_thickness=Q_(20, 'mm'),
    flange_thickness=Q_(15, 'mm'),
    orientation=Q_(-90, 'deg')  # turns the C-shape 90° degrees clockwise
)
section = Section(CShape, dim)

# Set resultant internal bending moment at section:
section.set_bending_moment(M_i.z)

# Get maximum bending stress at section:
print(section.bending.sigma_max.to('MPa'))

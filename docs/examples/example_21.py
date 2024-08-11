# Demo based on example 6.12 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.

# The simply supported beam in [Fig. 6–26a] has the cross-sectional area shown in
# [Fig. 6–26b]. Determine the absolute maximum bending stress in the beam and draw
# the stress distribution over the cross section at this location. Also, what is
# the stress at point B?

from mechanics import Quantity
from mechanics.statics import *
from mechanics.strength import Section
from mechanics.geometry import HShape, Dimensions


Q_ = Quantity

# %%
# Create beam with external loadings

# Hinge at position A
R_A = Force(
    magnitude='R_A',
    theta='theta_A',
    position=Position(0, units='m')
)

# Roller support at position B
R_B = Force(
    magnitude='R_B',
    theta=Angle(90),
    position=Position(6, units='m')
)

# Distributed rectangular load
q = DistributedLoad1D(
    x_coords=Q_([0, 6], 'm'),
    loads=Q_([-5, -5], 'kN / m')
)

beam = Beam(
    length=Q_(6, 'm'),
    loadings=[R_A, R_B, q],
    units=('kN', 'm'),
    num_sections=100
)

# Plot M-diagram
beam.moment_diagram.show()

# %%
# Maximum resultant internal bending moment (the accuracy depends on the number
# of sections that were set to draw the M-diagram).
print(beam.M_max())

# From the M-diagram it can be seen that the maximum lies at position x = 3 m.
F_i, M_i = beam.cut(x=Q_(3, 'm'))
print(M_i)

# %%
section = Section(
    HShape,
    dim=Dimensions(
        height=Q_(340, 'mm'),
        width=Q_(250, 'mm'),
        web_thickness=Q_(20, 'mm'),
        flange_thickness=Q_(20, 'mm')
    )
)
section.set_internal_loadings(F_i, M_i)
print(section.bending.sigma_max.to('MPa'))
print(section.bending.sigma(z=Q_(0, 'mm'), y=Q_(150, 'mm')).to('MPa'))

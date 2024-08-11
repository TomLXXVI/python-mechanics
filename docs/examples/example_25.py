# Demo based on example 6.16 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.

from mechanics import Quantity
from mechanics.statics import Moment
from mechanics.strength import Bending
from mechanics.geometry import ZShape

Q_ = Quantity


# Define bending moment about the y-axis:
M = Moment.create_from_components(
    vec_x=0.0,
    vec_y=20.0,
    vec_z=0.0,
    units='kN * m'
)


# Define the cross-section:
z_section = ZShape(
    width=Q_(600, 'mm'),
    height=Q_(400, 'mm'),
    web_thickness=Q_(100, 'mm'),
    flange_thickness=Q_(100, 'mm'),
    orientation=Q_(90, 'deg')
)
z_section.plot().show()


# Apply the bending moment to the cross-section:
bending = Bending(M=(M.z, M.y), shape=z_section)


# Maximum bending stress in cross-section:
print(bending.sigma_max.to('MPa'))


# Bending stresses at the vertices:
for i, vertex in enumerate(z_section.vertices):
    print(
        f"vertex {i+1} @ "
        f"({vertex[0].to('mm'):~P.0f}, {vertex[1].to('mm'):~P.0f}): "
        f"{bending.sigma(vertex[0], vertex[1]).to('MPa'):~P.2f}"
    )

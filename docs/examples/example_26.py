from mechanics import Quantity
from mechanics.statics import Force, Moment, DistributedLoad1D, Beam
from mechanics.strength import Section
from mechanics.geometry import Rectangle, Dimensions


Q_ = Quantity

# ------------------------------------------------------------------------------
# BEAM

R_A = Force(
    magnitude='R_A',
    theta='theta_A',
)

M_A = Moment.create_from_components(0.0, 0.0, None)

q = DistributedLoad1D(
    x_coords=Q_([0, 1], 'm'),
    loads=Q_([-1, -1], 'kN / m'),
    name='q'
)

beam = Beam(
    length=Q_(1, 'm'),
    loadings=[R_A, M_A, q],
    units=('kN', 'm'),
    num_sections=100
)

beam.moment_diagram.show()

F_i, M_i = beam.cut(x=Q_(0.5, 'm'), view='right')
print(M_i)

# ------------------------------------------------------------------------------
# CROSS-SECTION

section = Section(
    shape=Rectangle,
    dim=Dimensions(width=Q_(200, 'mm'), height=Q_(400, 'mm'))
)
section.set_internal_loadings(F_i, M_i)

section.shape.plot().show()

for corner in section.shape.vertices:
    print(
        f"@ ({corner[0]:~P.0f}, {corner[1]:~P.0f}): "
        f"{section.bending.sigma(corner[0], corner[1]).to('MPa'):~P.4f}"
    )

# Demo based on example 5.3 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.
from mechanics import Quantity
from mechanics.statics import Position, Angle, Force, Moment, System
from mechanics.geometry.shapes import Annulus, Dimensions
from mechanics.strength import Section


Q_ = Quantity


# Force parallel to the xy-plane (gamma = 0)
F1 = Force(
    position=Position(0, 0, -200, units='mm'),
    magnitude=80,
    theta=Angle(90, 'deg'),  # pointing in the positive y-direction (upward)
    units='N'
)
print(F1)


# Force parallel to the xy-plane (gamma = 0)
F2 = Force(
    position=Position(0, 0, 300, units='mm'),
    magnitude=80,
    theta=Angle(-90, 'deg'),  # pointing in the negative y-direction (downward)
    units='N'
)
print(F2)


# Torsion moment at point C about the x-axis (theta and gamma 0).
T_C = Moment(
    position=Position(300, units='mm'),
    magnitude='T_C',
    units='N * m'
)
print(T_C)


# Solve for the torsion moment at C
sys = System(loadings=[F1, F2, T_C])
sol_dict = sys.solve()
for name, sol in sol_dict.items():
    print(name, sol)


# Define cross-section at C (circular tube, i.e. an annulus)
sect_C = Section(
    shape=Annulus,
    dim=Dimensions(
        inner_radius=Q_(40, 'mm'),
        outer_radius=Q_(50, 'mm')
    )
)

# Get torsional shear stress at inner and outer tube wall:
sect_C.set_torque(sol_dict['T_C'].magnitude)
print(f"{sect_C.torsion.tau(r=sect_C.shape.dim.inner_radius).to('MPa'):~P.3f}")
print(f"{sect_C.torsion.tau(r=sect_C.shape.dim.outer_radius).to('MPa'):~P.3f}")

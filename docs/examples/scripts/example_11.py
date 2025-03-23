# Demo based on example 1.9 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.
from mechanics import Quantity
from mechanics.statics import (
    Position, Angle, Force,
    System
)
from mechanics.geometry.shapes import Circle, Dimensions
from mechanics.strength import Section

Q_ = Quantity

# %%
F_A = Force(
    magnitude='F_A',
    theta='theta_A',
    position=Position(0, units='m')
)

F_B = Force(
    magnitude='F_B',
    theta=Angle.create(4, 3, quadrant=1),
    position=Position(6, units='m')
)

F = Force(
    magnitude=30,
    theta=Angle(-90),
    position=Position(2, units='m'),
    units='kN',
    name='F'
)

sys = System(loadings=[F_A, F_B, F], units=('kN', 'm'))
sol = sys.solve()
for n, s in sol.items():
    print(f"{n}: {s}")

# %%
F_A = sol['F_A']

sec_A = Section(Circle, Dimensions(radius=Q_(20, 'mm') / 2))
sec_A.set_shear_force(F_A.magnitude / 2)  # double shear

print(sec_A.tau_V_avg.to('MPa'))

# %%
F_B = sol['F_B']

sec_B = Section(Circle, dim=Dimensions(radius=Q_(30, 'mm') / 2))
sec_B.set_shear_force(F_B.magnitude)  # single shear

print(sec_B.tau_V_avg.to('MPa'))

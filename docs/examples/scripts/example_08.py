# Demo based on example 1.6 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.
from mechanics import Quantity
from mechanics.statics import (
    Position, Angle, Force,
    System
)
from mechanics.strength import Section
from mechanics.geometry import Circle, Dimensions


Q_ = Quantity


F_BA = Force(
    magnitude='F_BA',
    theta=Angle(180 - 60),
    position=Position(0, 0),
)


F_BC = Force(
    magnitude='F_BC',
    theta=Angle.create(3, 4, quadrant=1),
    position=Position(0, 0),
)


W_lamp = Force(
    magnitude=9.81 * 80,
    theta=Angle(-90),
    position=Position(0, 0),
)

sys = System([F_BA, F_BC, W_lamp])
solutions = sys.solve()
for name, force in solutions.items():
    print(f"{name}: {force}")


F_BA = solutions['F_BA']
sec_BA = Section(Circle, Dimensions(radius=Q_(10 / 2, 'mm')))
sec_BA.set_normal_force(F_BA.magnitude)
sigma_BA = sec_BA.sigma_N_avg
print(sigma_BA.to('MPa'))


F_BC = solutions['F_BC']
sec_BC = Section(Circle, Dimensions(radius=Q_(8 / 2, 'mm')))
sec_BC.set_normal_force(F_BC.magnitude)
sigma_BC = sec_BC.sigma_N_avg
print(sigma_BC.to('MPa'))

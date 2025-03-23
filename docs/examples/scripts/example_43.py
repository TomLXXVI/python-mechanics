# Demo based on example 9.11 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.
from mechanics import Quantity
from mechanics.strength import PlaneStressTransformation


Q_ = Quantity


pst = PlaneStressTransformation(
    sigma_x=Q_(-20, 'MPa'),
    sigma_y=Q_(0, 'MPa'),
    tau_xy=Q_(-40, 'MPa')
)


sigma_1, sigma_2, theta_1 = pst.principal_normal_stresses()

print(f"{sigma_1.to('MPa'):~P.2f}")
print(f"{sigma_2.to('MPa'):~P.2f}")
print(f"{theta_1.to('deg'):~P.2f}")
print()

tau_abs_max, sigma_avg = pst.absolute_maximum_shear_stress()

print(f"{tau_abs_max.to('MPa'):~P.2f}")
print(f"{sigma_avg.to('MPa'):~P.2f}")

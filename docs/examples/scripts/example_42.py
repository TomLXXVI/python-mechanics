# Demo based on example 9.10 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.
from mechanics import Quantity
from mechanics.strength import PlaneStressTransformation


Q_ = Quantity


pst = PlaneStressTransformation(
    sigma_x=Q_(16, 'MPa'),
    sigma_y=Q_(32, 'MPa'),
    tau_xy=Q_(0, 'MPa')
)

sigma_1, sigma_2, _ = pst.principal_normal_stresses()

print(f"{sigma_1.to('MPa'):~P.2f}")
print(f"{sigma_2.to('MPa'):~P.2f}")


tau_abs_max, sigma_avg = pst.absolute_maximum_shear_stress()

print(f"{tau_abs_max.to('MPa'):~P.2f}")
print(f"{sigma_avg.to('MPa'):~P.2f}")


tau_max, sigma_avg, _ = pst.maximum_shear_stress()

print(f"{tau_max.to('MPa'):~P.2f}")
print(f"{sigma_avg.to('MPa'):~P.2f}")
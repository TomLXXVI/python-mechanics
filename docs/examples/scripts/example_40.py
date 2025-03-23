# Demo based on example 9.4 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.
from mechanics import Quantity
from mechanics.strength import PlaneStressTransformation


Q_ = Quantity


pst = PlaneStressTransformation(
    sigma_x=Q_(-20, 'MPa'),
    sigma_y=Q_(90, 'MPa'),
    tau_xy=Q_(60, 'MPa')
)

tau_max, sigma_avg, theta = pst.maximum_shear_stress()

print(f"{tau_max.to('MPa'):~P.2f}")
print(f"{sigma_avg.to('MPa'):~P.2f}")
print(f"{theta.to('deg'):~P.2f}")

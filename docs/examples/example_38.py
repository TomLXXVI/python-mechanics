# Demo based on example 9.2 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.
from mechanics import Quantity
from mechanics.strength import PlaneStressTransformation


Q_ = Quantity


pst = PlaneStressTransformation(
    sigma_x=Q_(-80, 'MPa'),
    sigma_y=Q_(50, 'MPa'),
    tau_xy=Q_(-25, 'MPa')
)

sigma_x, sigma_y, tau_xy = pst.transform(theta=Q_(-30, 'deg'))

print(f"{sigma_x.to('MPa'):~P.2f}")
print(f"{sigma_y.to('MPa'):~P.2f}")
print(f"{tau_xy.to('MPa'):~P.2f}")
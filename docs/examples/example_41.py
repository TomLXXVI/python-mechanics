# Demo based on example 9.7 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.
from mechanics import Quantity
from mechanics.strength import PlaneStressTransformation


Q_ = Quantity


pst = PlaneStressTransformation(
    sigma_x=Q_(-12, 'MPa'),
    sigma_y=Q_(0, 'MPa'),
    tau_xy=Q_(-6, 'MPa')
)

mohr_circle = pst.draw_mohr_circle()
mohr_circle.show()

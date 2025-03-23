# Demo based on example 5.3 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.

# A solid steel shaft AB is to be used to transmit 5 hp from the motor M to
# which it is attached. If the shaft rotates at v = 175 rpm and the steel has an
# allowable shear stress of tallow = 100 MPa determine the required diameter of
# the shaft to the nearest mm.

from mechanics import Quantity
from mechanics.geometry.shapes import Circle
from mechanics.strength.stress import Torsion


Q_ = Quantity

P = Q_(5, 'hp')
n = Q_(175, 'rpm')
T = (P / n).to('N * m')

r = Torsion.design(Circle, T, tau_allow=Q_(100, 'MPa'))
print(f"{r.to('mm'):~P.0f}")

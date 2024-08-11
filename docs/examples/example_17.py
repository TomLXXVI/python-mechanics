# Demo based on example 5.1 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.

# The solid shaft and tube shown in [Fig. 5–8] are made of a material having an
# allowable shear stress of 75 MPa. Determine the maximum torque that can be
# applied to each cross section, and show the stress acting on a small element
# of material at point A of the shaft, and points B and C of the tube.

from mechanics import Quantity
from mechanics.geometry.shapes import Circle, Annulus
from mechanics.strength.stress import Torsion


Q_ = Quantity

# Maximum allowable torque of shaft
circle = Circle(radius=Q_(100, 'mm'))
T_max_shaft = Torsion.allowable_torque(circle, tau_allow=Q_(75, 'MPa'))
print(f"{T_max_shaft.to('kN * m'):~P.1f}")

# The inverse: torque is given, determine required radius of the shaft.
r = Torsion.design(shape=Circle, T=T_max_shaft, tau_allow=Q_(75, 'MPa'))
print(f"{r.to('mm'):~P.0f}")


# Maximum allowable torque of tube
tube = Annulus(radius_in=Q_(75, 'mm'), radius_out=Q_(100, 'mm'))
T_max_tube = Torsion.allowable_torque(tube, tau_allow=Q_(75, 'MPa'))
print(f"{T_max_tube.to('kN * m'):~P.1f}")


# The inverse: torque is given, determine required outer radius of the tube
# with a given wall thickness.
r_o = Torsion.design(
    shape=Annulus,
    T=T_max_tube,
    t=Q_(25, 'mm'),
    tau_allow=Q_(75, 'MPa')
)
print(f"{r_o.to('mm'):~P.0f}")


# Shear stress at the inner radius of the tube.
torsion = Torsion(T=T_max_tube, shape=tube)
print(torsion.tau(r=Q_(75, 'mm')).to('MPa'))

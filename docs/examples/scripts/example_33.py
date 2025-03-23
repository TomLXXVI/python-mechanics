# Demo based on example 8.1 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.
from mechanics import Quantity
from mechanics.strength.thin_wall_vessel import CylindricalVessel, SphericalVessel

Q_ = Quantity


cyl_vessel = CylindricalVessel(
    inner_radius=Q_(1.2 / 2, 'm'),
    thickness=Q_(12, 'mm'),
    sigma_allow=Q_(140, 'MPa')
)
print(cyl_vessel.maximum_pressure.to('MPa'))


sph_vessel = SphericalVessel(
    inner_radius=Q_(1.2 / 2, 'm'),
    thickness=Q_(12, 'mm'),
    sigma_allow=Q_(140, 'MPa')
)
print(sph_vessel.maximum_pressure.to('MPa'))

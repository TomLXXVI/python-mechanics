# Demo based on example 8.3 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.
from math import pi
from mechanics import Quantity
from mechanics.strength.thin_wall_vessel import CylindricalVessel

Q_ = Quantity


vessel = CylindricalVessel(
    inner_radius=Q_(600, 'mm'),
    thickness=Q_(10, 'mm'),
    pressure=Q_(450, 'kPa')
)

sigma_cir = vessel.circumferential_stress
sigma_lon = vessel.longitudinal_stress

print(f"circumferential stress: {sigma_cir.to('MPa'):~P.1f}")
print(f"longitudinal stress: {sigma_lon.to('MPa'):~P.1f}")

# Axial load at the top of the vessel exerts a compressive force, while the
# longitudinal stress due to the gas pressure is a tensile stress.
load = Q_(200, 'kN')
A_sect = pi * (vessel.r_o ** 2 - vessel.r_i ** 2)
sigma_load = load / A_sect
print(f"axial stress due to load: {sigma_load.to('MPa'):~P.2f}")

sigma_res = sigma_lon - sigma_load
print(f"resultant axial stress: {sigma_res.to('MPa'):~P.2f}")

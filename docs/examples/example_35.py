# Demo based on example 8.4 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.

from mechanics import Quantity
from mechanics.statics import (
    Position,
    Angle,
    Force,
    DistributedLoad1D,
    System,
    Beam
)
from mechanics.strength import Section
from mechanics.geometry import Rectangle, Dimensions


Q_ = Quantity


# External loadings

R_A = Force(
    magnitude='R_A',
    theta='theta_A',
    position=Position(0, 0, units='m'),
    units='kN'
)

R_B = Force(
    magnitude='R_B',
    theta=Angle.create(1.5, 2, quadrant=4) + Angle(90),
    position=Position(6, -1.5, units='m'),
    units='kN'
)

q = DistributedLoad1D(
    x_coords=Q_([4, 6.5], 'm'),
    loads=Q_([-50, -50], 'kN / m'),
    slope=Angle.create(1.5, 2, quadrant=4),
    name='q'
)

system = System(loadings=[R_A, R_B, q], units=('kN', 'm'))
reaction_forces = system.solve()
for name, force in reaction_forces.items():
    print(f"{name} {force}, magnitude = {force.magnitude.to('kN'):~P}")


# Internal loadings

beam = Beam(
    length=Q_(4, 'm'),
    loadings=[reaction_forces['R_A']],
    units=('kN', 'm')
)

F_i, M_i = beam.cut(x=Q_(1.5, 'm'), view='left')

F_i.name = 'F_i at pos C'
M_i.name = 'M_i at pos C'

print(f"{F_i.name}: {F_i}")
print(f"{M_i.name}: {M_i}")


# Stresses in cross-section at C

section = Section(
    shape=Rectangle,
    dim=Dimensions(width=Q_(50, 'mm'), height=Q_(250, 'mm'))
)

section.set_internal_loadings(F_i, M_i)

sigma_N = section.axial.sigma
tau_V = section.transverse_shear.tau(y=Q_(125, 'mm'))
sigma_M = section.bending.sigma(z=Q_(-25, 'mm'), y=Q_(125, 'mm'))

print(
    f"normal stress at C: {sigma_N.to('MPa'):~P.2f}",
    f"shear stress at C: {tau_V.to('MPa'):~P.2f}",
    f"bending stress at C: {sigma_M.to('MPa'):~P.2f}",
    sep='\n'
)

# Stresses in cross-section at D

sigma_N = section.axial.sigma
tau_V = section.transverse_shear.tau(y=Q_(0, 'mm'))
sigma_M = section.bending.sigma(z=Q_(-25, 'mm'), y=Q_(0, 'mm'))

print(
    f"normal stress at D: {sigma_N.to('MPa'):~P.2f}",
    f"shear stress at D: {tau_V.to('MPa'):~P.2f}",
    f"bending stress at D: {sigma_M.to('MPa'):~P.2f}",
    sep='\n'
)

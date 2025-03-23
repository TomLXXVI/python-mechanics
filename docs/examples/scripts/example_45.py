# Demo based on example 10.12 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.
from mechanics import Quantity
from mechanics.statics import Position, Angle, Force, Moment, Beam
from mechanics.strength import (
    Section,
    PlaneStressTransformation,
    VonMisesYieldCriterion
)
from mechanics.geometry import Circle, Dimensions


Q_ = Quantity

# Determine resultant internal force and moment in the cross-section:

F = Force(
    magnitude=65,
    theta=Angle(0),
    position=Position(0, units='m'),
    units='kN',
    name='F'
)

T = Moment(
    magnitude=350,
    theta=Angle(0),
    position=Position(0, units='m'),
    units='N * m',
    name='T'
)

beam = Beam(
    length=Q_(1, 'm'),
    loadings=[F, T],
    units=('kN', 'm')
)

F_i, M_i = beam.cut(x=Q_(0.5, 'm'), view='left')
print(f"resultant internal force = {F_i}")
print(f"resultant internal moment = {M_i}")

# Determine maximum stresses in the cross-section:

section = Section(shape=Circle, dim=Dimensions(radius=Q_(12.5, 'mm')))
section.set_internal_loadings(F_i, M_i)
if section.axial:
    print(f"axial stress: {section.axial.sigma.to('MPa'):~P.2f}")
if section.torsion:
    print(f"torsion shear stress: {section.torsion.tau_max.to('MPa'):~P.2f}")

# Get principal stresses:

pst = PlaneStressTransformation(
    sigma_x=section.axial.sigma,
    sigma_y=Q_(0, 'MPa'),
    tau_xy=section.torsion.tau_max
)
sigma1, sigma2, theta1 = pst.principal_normal_stresses()
print(f"maximum normal stress = {sigma1.to('MPa'):~P.2f}")
print(f"minimum normal stress = {sigma2.to('MPa'):~P.2f}")

# Apply von Mises yield criterion (maximum distortion energy theory):

vmc = VonMisesYieldCriterion(sigma_yield=Q_(250, 'MPa'))
vmc.draw_chart(sigma1, sigma2, units='MPa').show()
print(f"von Mises stress = {vmc.von_mises_stress(sigma1, sigma2).to('MPa'):~P.2f}")
if vmc.check(sigma1, sigma2) is True:
    print(f"shear failure will not occur.")
else:
    print(f"shear failure will occur.")

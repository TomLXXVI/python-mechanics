# Demo based on example 5.2 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.
from mechanics import Quantity
from mechanics.statics import Position, Angle, Moment, Beam
from mechanics.strength import Section
from mechanics.geometry import Circle, Dimensions


Q_ = Quantity

# External loadings:

T1 = Moment(
    position=Position(0.5, units='m'),
    magnitude=4.25,
    theta=Angle(0),
    units='kN * m',
    name='T1'
)
print(T1)


T2 = Moment(
    position=Position(1, units='m'),
    magnitude=3.0,
    theta=Angle(180),
    units='kN * m',
    name='T2'
)
print(T2)


T3 = Moment(
    position=Position(1.5, units='m'),
    magnitude=1.25,
    theta=Angle(180),
    units='kN * m',
    name='T3'
)
print(T3)


# Definition of the shaft:
shaft = Beam(length=Q_(2, 'm'), loadings=[T1, T2, T3])


# Resultant internal torque at section a-a:
F_i, M_i = shaft.cut(Q_(1.25, 'm'))
print(M_i.x.to('kN * m'))


# Torsional shear stress at outer radius (maximum):
S_aa = Section(shape=Circle, dim=Dimensions(radius=Q_(75, 'mm')))
S_aa.set_internal_loadings(F_i, M_i)
print(f"{S_aa.torsion.tau_max.to('MPa'):~P.2f}")


# Torsional shear stress at intermediate radius:
print(f"{S_aa.torsion.tau(r=Q_(15, 'mm')).to('MPa'):~P.2f}")

# Demo based on example 4.1 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.

# Note: From here on, the `Beam` class in subpackage `strength` is used.

from mechanics import Quantity
from mechanics.statics import Position, Angle, Force
from mechanics.strength import Beam
from mechanics.geometry import Circle, Dimensions


Q_ = Quantity


F_A = Force(
    magnitude='F_A',
    theta=Angle(0),
    position=Position(0, units='m'),
    units='kN'
)

F_B = Force(
    magnitude=80,
    theta=Angle(0),
    position=Position(2, units='m'),
    units='kN',
    name='F_B'
)

F_C = Force(
    magnitude=40,
    theta=Angle(0),
    position=Position(3, units='m'),
    units='kN',
    name='F_C'
)

F_D = Force(
    magnitude=70,
    theta=Angle(180),
    position=Position(x=4.5, units='m'),
    units='kN',
    name='F_D'
)

beam = Beam(
    length=Q_(4.5, 'm'),
    shape_type=Circle,
    shape_dim=Dimensions(radius=Q_(25, 'mm')),
    E_modulus=Q_(200, 'GPa'),
    G_modulus=Q_(75, 'GPa'),
    loadings=[F_A, F_B, F_C, F_D],
    units=('kN', 'm'),
    num_sections=50
)

delta = beam.elongation(x1=Q_(0, 'm'), x2=Q_(4.5, 'm'))
print(delta.to('mm'))

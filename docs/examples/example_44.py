# Demo based on example 11.1 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.
from mechanics import Quantity
from mechanics.statics import Position, Angle, Force, Beam, DistributedLoad1D
from mechanics.strength import TransverseShear
from mechanics.geometry import HShape


Q_ = Quantity


# External loadings

F_1 = Force(
    magnitude=120,
    theta=Angle(-90),
    position=Position(2, units='m'),
    units='kN',
    name='F_1'
)

F_2 = Force(
    magnitude=60,
    theta=Angle(-90),
    position=Position(6, units='m'),
    units='kN',
    name='F_2'
)

R_A = Force(
    magnitude='R_A',
    theta='theta_A',
    position=Position(0, units='m'),
    units='kN'
)

R_B = Force(
    magnitude='R_B',
    theta=Angle(90),
    position=Position(4, units='m')
)

# Beam with shear and moment diagrams

beam1 = Beam(
    length=Q_(6, 'm'),
    loadings=[F_1, F_2, R_A, R_B],
    units=('kN', 'm'),
    num_sections=100
)

beam1.shear_diagram.show()
beam1.moment_diagram.show()
x_V_max, V_max = beam1.V_max()
x_M_max, M_max = beam1.M_max()
print(f"V_max = {V_max:~P.0f} at x = {x_V_max:~P.3f}")
print(f"M_max = {M_max:~P.0f} at x = {x_M_max:~P.3f}")

sigma_allow = Q_(165, 'MPa')
S_req = abs(M_max) / sigma_allow
print(f"required section modulus = {S_req.to('mm ** 3'):~P.3e}")
print()

# Select actual beam...(W410 x 46 with S = 774e3 mm**3)

# Check actual beam including its weight

w = DistributedLoad1D(
    x_coords=Q_([0, 6], 'm'),
    loads=Q_([-46, -46], 'kg / m') * Q_(9.81, 'N / kg')
)

beam2 = Beam(
    length=Q_(6, 'm'),
    loadings=[F_1, F_2, R_A, R_B, w],
    units=('kN', 'm'),
    num_sections=100
)

x_V_max, V_max = beam2.V_max()
x_M_max, M_max = beam2.M_max()
print(f"V_max = {V_max:~P.0f} at x = {x_V_max:~P.3f}")
print(f"M_max = {M_max:~P.0f} at x = {x_M_max:~P.3f}")
S_req = abs(M_max) / sigma_allow
print(f"required section modulus including beam weight = {S_req.to('mm ** 3'):~P.3e}")


W410x46 = HShape(
    height=Q_(403, 'mm'),
    width=Q_(140, 'mm'),
    web_thickness=Q_(6.99, 'mm'),
    flange_thickness=Q_(11.2)
)

shear_obj = TransverseShear(V_max, W410x46)
print(f"maximum shear stress = {shear_obj.tau_max.to('MPa'):~P.2f}")

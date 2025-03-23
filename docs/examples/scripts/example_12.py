# Demo based on example 1.11 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.
from mechanics import Quantity
from mechanics.statics import Angle, Force
from mechanics.strength import Section
from mechanics.geometry import Rectangle, Dimensions


Q_ = Quantity


# %%
F = Force(
    magnitude=3000,
    theta=Angle.create(4, 3, quadrant=3),
    units='N'
)
print(F.x)
print(F.y)

# %%
# Average compressive stress in section AB
section_AB = Section(
    Rectangle,
    Dimensions(width=Q_(40, 'mm'), height=Q_(25, 'mm'))
)
section_AB.set_normal_force(F.x)
print(f"sigma_AB: {section_AB.sigma_N_avg.to('MPa')}")

# %%
# Average compressive stress in section BC
section_BC = Section(
    Rectangle,
    Dimensions(width=Q_(50, 'mm'), height=Q_(40, 'mm'))
)
section_BC.set_normal_force(F.y)
print(f"sigma_BC: {section_BC.sigma_N_avg.to('MPa')}")

# %%
# Average shear stress in section DB
section_DB = Section(
    Rectangle,
    Dimensions(width=Q_(75, 'mm'), height=Q_(40, 'mm'))
)
section_DB.set_shear_force(F.x)
print(f"tau_DB: {section_DB.tau_V_avg.to('MPa')}")

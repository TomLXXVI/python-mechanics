# Demo of `Rectangle` class
from mechanics import Quantity
from mechanics.geometry import Rectangle

Q_ = Quantity

rect = Rectangle(width=Q_(300, 'mm'), height=Q_(500, 'mm'))

print(f"{rect.area.to('m**2'):~P}")

print(f"{rect.moment_of_inertia_xx:~P.2e}")
print(f"{rect.moment_of_inertia_yy:~P.2e}")
print(f"{rect.product_of_inertia:~P.2e}")
print(f"{rect.polar_moment_of_inertia:~P.2e}")

I_max, I_min, theta = rect.principal_moments_of_inertia
print(f"{I_max:~P.2e}, {I_min:~P.2e}, {theta.to('deg'):~P}")

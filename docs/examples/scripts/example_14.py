# Demo of `Polygon` class
from mechanics import Quantity
from mechanics.geometry.shapes import Polygon


Q_ = Quantity

# Specify the vertices of the polygon in a consecutive, counterclockwise order.
poly = Polygon([
    Q_([0, 0], 'm'),
    Q_([5, 0], 'm'),
    Q_([5, 1], 'm'),
    Q_([3.125, 1], 'm'),
    Q_([2.125, 3], 'm'),
    Q_([0.875, 3], 'm'),
    Q_([1.875, 1], 'm'),
    Q_([0, 1], 'm')
])

print(poly.area)
print(poly.moment_of_inertia_xx)
print(poly.moment_of_inertia_yy)
print(poly.product_of_inertia)
print(poly.polar_moment_of_inertia)

I_max, I_min, theta = poly.principal_moments_of_inertia
print(I_max)
print(I_min)
print(theta.to('deg'))


diagram = poly.plot()
diagram.show()

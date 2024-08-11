# Demo of `Polygon` class
from mechanics import Quantity
from mechanics.geometry.shapes import Polygon


Q_ = Quantity


poly = Polygon([
    Q_([0, 0], 'mm'),
    Q_([0, 400], 'mm'),
    Q_([100, 400], 'mm'),
    Q_([100, 100], 'mm'),
    Q_([600, 100], 'mm'),
    Q_([600, -300], 'mm'),
    Q_([500, -300], 'mm'),
    Q_([500, 0], 'mm'),
    Q_([0, 0], 'mm')
])

I_max, I_min, theta = poly.principal_moments_of_inertia
print(f"{theta.to('deg'):~P.1f}")
print(f"{I_max:~P.2e}")
print(f"{I_min:~P.2e}")

diagram = poly.plot()
diagram.show()

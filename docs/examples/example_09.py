# Demo based on example 1.7 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.
from mechanics import Quantity
from mechanics.geometry.shapes import Circle
from mechanics.strength.stress import AxialLoading


Q_ = Quantity


d = Q_(200, 'mm')
circle = Circle(radius=d/2)


h1 = Q_(800, 'mm')
V1 = circle.area * h1


rho = Q_(7850, 'kg / m ** 3')
m1 = rho * V1.to('m ** 3')
G1 = Q_(9.81, 'm / s ** 2') * m1


sigma = AxialLoading(N=G1, shape=circle).sigma
print(sigma.to('kN / m ** 2'))

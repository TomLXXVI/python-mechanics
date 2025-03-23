# Demo based on example 1.8 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.
import numpy as np
from scipy.optimize import root_scalar
from mechanics import Quantity

from mechanics.statics import Position, Angle, Force, System
from mechanics.charts import LineChart
from mechanics.strength import AxialLoading

Q_ = Quantity


def get_stresses(x: Quantity) -> tuple[Quantity, Quantity]:
    # Get stresses in tie rod AB and smooth support C when the vertical force
    # is acting at position `x`.
    F_AB = Force(
        position=Position(0, 0, units='mm'),
        magnitude='F_AB',
        theta=Angle(90)
    )
    F_C = Force(
        position=Position(200, 0, units='mm'),
        magnitude='F_C',
        theta=Angle(90)
    )
    F = Force(
        magnitude=3,
        theta=Angle(-90),
        units='kN',
        position=Position(x.to('mm').magnitude, 0, units='mm')
    )
    sys = System([F_AB, F_C, F])
    solutions = sys.solve()
    F_AB = solutions['F_AB']
    F_C = solutions['F_C']

    sigma_AB = AxialLoading(
        N=F_AB.magnitude,
        A=Q_(400, 'mm ** 2')
    ).sigma.to('MPa')

    sigma_C = AxialLoading(
        N=F_C.magnitude,
        A=Q_(650, 'mm ** 2')
    ).sigma.to('MPa')

    return sigma_AB, sigma_C


def find_position() -> Quantity:
    # Find the position of vertical force for which the stresses in tie rod AB
    # and at smooth support C have equal magnitude.
    x_min = Q_(0, 'mm')
    x_max = Q_(200, 'mm')

    def __fun(x: float) -> float:
        x = Q_(x, 'mm')
        sigma_AB, sigma_C = get_stresses(x)
        delta = abs(sigma_AB) - abs(sigma_C)
        return delta.m

    sol = root_scalar(__fun, bracket=(x_min.m, x_max.m))
    x_sol = Q_(sol.root, 'mm')
    return x_sol


def main():
    # Calculate the stresses in tie rod AB and at smooth support C for a range
    # of positions of the vertical force F.
    x_arr = Q_(np.linspace(0, 200), 'mm')
    sigma_AB_arr, sigma_C_arr = zip(*[get_stresses(x) for x in x_arr])
    sigma_AB_arr = Quantity.from_list(sigma_AB_arr)
    sigma_C_arr = Quantity.from_list(sigma_C_arr)

    diagram = LineChart()
    diagram.add_xy_data(
        label='sigma_AB',
        x1_values=x_arr.m,
        y1_values=sigma_AB_arr.m
    )
    diagram.add_xy_data(
        label='sigma_C',
        x1_values=x_arr.m,
        y1_values=sigma_C_arr.m
    )
    diagram.x1.add_title('position x, mm')
    diagram.y1.add_title('sigma, MPa')
    diagram.show()

    # Find the position of vertical force for which the stresses in tie rod AB
    # and at smooth support C have equal magnitude.
    x_equal_stress = find_position()
    print(x_equal_stress)


if __name__ == '__main__':
    main()

from collections.abc import Sequence
from scipy.interpolate import interp1d
from scipy.integrate import quad
from mechanics import Quantity
from .force import Force

Q_ = Quantity


class DistributedLoad2D:
    """
    A distributed piece-wise linear load acting perpendicular to a 2D-beam.
    """
    def __init__(
        self,
        points: Sequence[tuple[Quantity, Quantity]],
        name: str = ''
    ) -> None:
        """Creates a `DistributedLoad2D` instance.

        Parameters
        ----------
        points:
            List of 2-tuples that lie on the piece-wise linear load profile.
            The first element of a tuple is the position of a point along
            the x-axis (longitudinal axis) of the beam. The second element is
            the force per unit of length at this position.
            The first tuple in the list is where the distributed load starts on
            the beam. The last tuple is where the distributed load ends on the
            beam.
        """
        self.points = points
        self.name = name
        self.__set_units()
        self.__create_load_interpolant()

    def resultant(
        self,
        x1: Quantity | None = None,
        x2: Quantity | None = None
    ) -> Force:
        """Returns the resulting force of the distributed load between positions
        `x1` and `x2` along the x-axis of the beam (with `x2` > `x1`).
        """
        x1 = (
            x1.to(self.u_pos).m
            if x1 is not None
            else self.points[0][0].to(self.u_pos).m
        )
        x2 = (
            x2.to(self.u_pos).m
            if x2 is not None
            else self.points[-1][0].to(self.u_pos).m
        )
        # noinspection PyTypeChecker
        Q_mag = quad(self.q_x, x1, x2)[0]
        x_c = quad(lambda x: x * self.q_x(x), x1, x2)[0] / Q_mag
        Q = self.__create_force(Q_mag, x_c)
        return Q

    def __set_units(self) -> None:
        self.u_pos = self.points[0][0].units
        self.u_load = self.points[0][1].units

    def __create_load_interpolant(self):
        x_data, q_data = zip(*[
            (x.to(self.u_pos).m, q.to(self.u_load).m)
            for x, q in self.points
        ])
        self.q_x = interp1d(x_data, q_data)

    def __create_force(self, Q: float, x_c: float) -> Force:
        if Q > 0:
            theta = Q_(90, 'deg')
        elif Q < 0:
            theta = Q_(-90, 'deg')
            Q = -Q
        else:
            theta = Q_(0, 'deg')
        Q = Force(
            action_point=(Q_(x_c, self.u_pos), Q_(0, self.u_pos)),
            magnitude=Q_(Q, self.u_pos * self.u_load),
            theta=theta
        )
        return Q

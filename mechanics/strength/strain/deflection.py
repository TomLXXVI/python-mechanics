from typing import Callable
from dataclasses import dataclass
from abc import ABC, abstractmethod
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
from mechanics import Quantity
from mechanics.charts import LineChart


Q_ = Quantity


@dataclass
class BoundaryCondition:
    x: Quantity
    y: Quantity | None = None
    theta: Quantity | None = None


class _Case(ABC):

    def __init__(
        self,
        M: Callable[[float], float],
        E: float,
        I: float,
        num_intervals: int
    ) -> None:
        self._M = M
        self._E = E
        self._I = I
        self._n = num_intervals
        self._h = 0.0

    def _calc_slope(self, y: np.ndarray) -> np.ndarray:
        y_0_fw = y[:-2]
        y_1_fw = y[1:-1]
        y_2_fw = y[2:]
        theta_fw = (-3 * y_0_fw + 4 * y_1_fw - 1 * y_2_fw) / (2 * self._h)
        y_0_bw = y[-1:1:-1]
        y_1_bw = y[-2:0:-1]
        y_2_bw = y[-3::-1]
        theta_bw = (3 * y_0_bw - 4 * y_1_bw + 1 * y_2_bw) / (2 * self._h)
        head = theta_fw[:2]
        tail = theta_bw[1::-1]
        body = (theta_fw[2:] + theta_bw[-3::-1]) / 2
        theta = np.concatenate((head, body, tail))
        return theta

    @abstractmethod
    def solve(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        pass


class _Case1(_Case):
    # case 1: two boundary conditions - displacements y_a and y_b given.
    def __init__(
        self,
        M: Callable[[float], float],
        E: float,
        I: float,
        bc1: BoundaryCondition,
        bc2: BoundaryCondition,
        num_intervals: int
    ) -> None:
        super().__init__(M, E, I, num_intervals)
        x_a = bc1.x.to('m').magnitude
        x_b = bc2.x.to('m').magnitude
        self._h = (x_b - x_a) / self._n
        self._x = np.arange(x_a, x_b + self._h, self._h)
        self._y_a = bc1.y.to('m').magnitude
        self._y_b = bc2.y.to('m').magnitude

    @staticmethod
    def __create_A(n: int) -> np.ndarray:
        A = np.zeros((n + 1, n + 1))
        A[0, 0] = 1
        for i in range(1, n):
            A[i, i - 1] = 1
            A[i, i] = -2
            A[i, i + 1] = 1
        A[n, n] = 1
        return A

    @staticmethod
    def __create_B(
        n: int,
        h: float,
        x: np.ndarray,
        y_a: float,
        y_b: float,
        f: Callable[[float | np.ndarray], float]
    ) -> np.ndarray:
        B = np.zeros(n + 1)
        B[0] = y_a
        for i in range(1, n):
            B[i] = h ** 2 * f(x[i])
        B[n] = y_b
        return B

    def solve(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        A = self.__create_A(self._n)
        B = self.__create_B(
            self._n, self._h, self._x,
            self._y_a, self._y_b,
            lambda x: self._M(x) / (self._E * self._I)
        )
        y = np.linalg.solve(A, B)
        theta = self._calc_slope(y)
        return self._x, y, theta


class _Case2(_Case):
    # case 2: two boundary conditions - displacement y_a and slope theta_b given
    def __init__(
        self,
        M: Callable[[float], float],
        E: Quantity,
        I: Quantity,
        bc1: BoundaryCondition,
        bc2: BoundaryCondition,
        num_intervals: int
    ) -> None:
        super().__init__(M, E, I, num_intervals)
        x_a = bc1.x.to('m').magnitude
        x_b = bc2.x.to('m').magnitude
        self._h = (x_b - x_a) / self._n
        self._x = np.arange(x_a, x_b + self._h, self._h)
        self._y_a = bc1.y.to('m').magnitude
        self._theta_b = bc2.theta.to('rad').magnitude

    @staticmethod
    def __create_A(n: int) -> np.ndarray:
        A = np.zeros((n + 1, n + 1))
        A[0, 0] = 1
        for i in range(1, n):
            A[i, i - 1] = 1
            A[i, i] = -2
            A[i, i + 1] = 1
        A[n, n - 1] = 2
        A[n, n] = -2
        return A

    @staticmethod
    def __create_B(
        n: int,
        h: float,
        x: np.ndarray,
        y_a: float,
        theta_b: float,
        f: Callable[[float | np.ndarray], float]
    ) -> np.ndarray:
        B = np.zeros(n + 1)
        B[0] = y_a
        for i in range(1, n):
            B[i] = h ** 2 * f(x[i])
        B[n] = h ** 2 * f(x[n]) - 2 * h * theta_b
        return B

    def solve(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        A = self.__create_A(self._n)
        B = self.__create_B(
            self._n, self._h, self._x,
            self._y_a, self._theta_b,
            lambda x: self._M(x) / (self._E * self._I)
        )
        y = np.linalg.solve(A, B)
        theta = self._calc_slope(y)
        return self._x, y, theta


class _Case3(_Case):
    # case 3: two boundary conditions - slope theta_a and displacement y_b given
    def __init__(
        self,
        M: Callable[[Quantity], Quantity],
        E: Quantity,
        I: Quantity,
        bc1: BoundaryCondition,
        bc2: BoundaryCondition,
        num_intervals: int
    ) -> None:
        super().__init__(M, E, I, num_intervals)
        x_a = bc1.x.to('m').magnitude
        x_b = bc2.x.to('m').magnitude
        self._h = (x_b - x_a) / self._n
        self._x = np.arange(x_a, x_b + self._h, self._h)
        self._theta_a = bc1.theta.to('rad').magnitude
        self._y_b = bc2.y.to('m').magnitude

    @staticmethod
    def __create_A(n: int) -> np.ndarray:
        A = np.zeros((n + 1, n + 1))
        A[0, 0] = -2
        A[0, 1] = 2
        for i in range(2, n):
            A[i, i - 1] = 1
            A[i, i] = -2
            A[i, i + 1] = 1
        A[n, n] = 1
        return A

    @staticmethod
    def __create_B(
        n: int,
        h: float,
        x: np.ndarray,
        theta_a: float,
        y_b: float,
        f: Callable[[float | np.ndarray], float]
    ) -> np.ndarray:
        B = np.zeros(n + 1)
        B[0] = h ** 2 * f(x[0]) + 2 * h * theta_a
        for i in range(1, n):
            B[i] = h ** 2 * f(x[i])
        B[n] = y_b
        return B

    def solve(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        A = self.__create_A(self._n)
        B = self.__create_B(
            self._n, self._h, self._x,
            self._theta_a, self._y_b,
            lambda x: self._M(x) / (self._E * self._I)
        )
        y = np.linalg.solve(A, B)
        theta = self._calc_slope(y)
        return self._x, y, theta


class _Case4(_Case):
    # case 4: single boundary condition - displacement y_a and slope theta_a given
    def __init__(
        self,
        M: Callable[[Quantity], Quantity],
        E: Quantity,
        I: Quantity,
        bc1: BoundaryCondition,
        bc2: BoundaryCondition,
        num_intervals: int
    ) -> None:
        super().__init__(M, E, I, num_intervals)
        self._x_a = bc1.x.to('m').magnitude
        self._x_b = bc2.x.to('m').magnitude
        self._h = (self._x_b - self._x_a) / self._n
        self._x = np.round(np.arange(self._x_a, self._x_b + self._h, self._h), 12)
        self._y_a = bc1.y.to('m').magnitude
        self._theta_a = bc1.theta.to('rad').magnitude

    def F(self, x: float, y: np.ndarray) -> np.ndarray:
        y_dot = np.zeros(2)
        y_dot[0] = y[1]
        y_dot[1] = self._M(x) / (self._E * self._I)
        return y_dot

    def solve(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        sol = solve_ivp(
            self.F,
            t_span=(self._x_a, self._x_b),
            y0=[self._y_a, self._theta_a],
            t_eval=self._x
        )
        y = sol.y[0]
        theta = self._calc_slope(y)
        return self._x, y, theta


class ElasticCurve:
    """
    Determine the deflection and slope of beams and shafts.
    Assumptions: the beam is initially straight, it is elastically deformed by
    the loads, such that the slope and deflection of the elastic curve are
    very small, and the deformations are only caused by bending.
    """
    def __init__(
        self,
        M: Callable[[Quantity], Quantity],
        E: Quantity,
        I: Quantity,
        boundary_conditions: list[BoundaryCondition],
        num_intervals: int = 50
    ) -> None:
        """Creates an `ElasticCurve` object.

        Parameters
        ----------
        M:
            Function that takes the position `x` along the beam and returns the
            resultant internal bending moment in the cross-section of the beam
            at this position.
        E:
            Modulus of elasticity of the material.
        I:
            Area moment of inertia of the cross-section about the neutral axis.
        boundary_conditions:
            List of boundary conditions (instances of class `BoundaryCondition`)
            where the vertical displacement and/or slope of the elastic curve
            are fixed, i.e. at the beam supports. In case of a roller or hinge
            the vertical displacement is zero, but the slope is unknown. In case
            of a fixed end the slope is also zero. In case of a free end, both
            the vertical displacement and the slope are unknown.
        num_intervals:
            Number of intervals where the differential equation of the elastic
            curve is to be evaluated.
        """
        self.M = M
        self.E = E
        self.I = I
        # make sure boundary conditions are sorted from left to right along the
        # longitudinal axis of the beam
        self.bound_conds = sorted(boundary_conditions, key=lambda bc: bc.x)
        self.num_intervals = self._n = num_intervals

        _M = self.__create_M_function()
        _E = self.E.to('Pa').magnitude
        _I = self.I.to('m**4').magnitude

        # group boundary conditions in pairs
        self.bound_conds_grouped = [
            (self.bound_conds[i], self.bound_conds[i + 1])
            for i in range(len(self.bound_conds) - 1)
        ]

        # between each pair of boundary conditions, the elastic curve is
        # calculated (vertical displacement and slope); depending on the
        # boundary conditions, 4 different cases can be distinguished that
        # require a different method to solve the differential equation of the
        # elastic curve
        x, y, theta = [], [], []
        for group in self.bound_conds_grouped:
            bc1, bc2 = group
            if bc1.y is not None and bc2.y is not None:
                case1 = _Case1(_M, _E, _I, bc1, bc2, num_intervals)
                x_1, y_1, theta_1 = case1.solve()
                x.append(x_1)
                y.append(y_1)
                theta.append(theta_1)
                bc1.theta = Q_(theta_1[0], 'rad')
                bc2.theta = Q_(theta_1[-1], 'rad')
            elif bc1.y is not None and bc2.theta is not None:
                case2 = _Case2(_M, _E, _I, bc1, bc2, num_intervals)
                x_2, y_2, theta_2 = case2.solve()
                x.append(x_2)
                y.append(y_2)
                theta.append(theta_2)
                bc1.theta = Q_(theta_2[0], 'rad')
                bc2.y = Q_(y_2[-1], 'm')
            elif bc1.theta is not None and bc2.y is not None:
                case3 = _Case3(_M, _E, _I, bc1, bc2, num_intervals)
                x_3, y_3, theta_3 = case3.solve()
                x.append(x_3)
                y.append(y_3)
                theta.append(theta_3)
                bc1.y = Q_(y_3[0], 'm')
                bc2.theta = Q_(theta_3[-1], 'rad')
            elif bc1.y is not None and bc1.theta is not None:
                case4 = _Case4(_M, _E, _I, bc1, bc2, num_intervals)
                x_4, y_4, theta_4 = case4.solve()
                x.append(x_4)
                y.append(y_4)
                theta.append(theta_4)
                bc2.y = Q_(y_4[-1], 'm')
                bc2.theta = Q_(theta_4[-1], 'rad')
        if y:
            x = np.concatenate(x)
            y = np.concatenate(y)
            theta = np.concatenate(theta)
        else:
            raise ValueError('boundary conditions are wrong...')

        self._y_interp = interp1d(x, y)
        self._theta_interp = interp1d(x, theta)
        self.x = Q_(x, 'm')  # array of positions where vertical displacement and slope are evaluated.
        self.y = Q_(y, 'm')  # array of vertical displacements
        self.theta = Q_(theta, 'rad')  # array of slopes

    def __create_M_function(self) -> Callable[[float], float]:
        def M(x: float) -> float:
            x = Q_(x, 'm')
            M = self.M(x)
            M = M.to('N * m').magnitude
            return M
        return M

    def displacement(self, x: Quantity) -> Quantity:
        """Returns the vertical displacement of the elastic curve at
        position `x`.
        """
        y = self._y_interp(x.to('m').magnitude)
        return Q_(y, 'm')

    def slope(self, x: Quantity) -> Quantity:
        """Returns the slope of the elastic curve at position `x`."""
        theta = self._theta_interp(x.to('m').magnitude)
        return Q_(theta, 'rad')

    def diagram(self, units: tuple[str, str] = ('m', 'mm')) -> LineChart:
        """Returns a `LineChart` object with a diagram of the elastic curve.
        Through parameter `units` the display units can be set for the position
        `x` along the longitudinal axis of the beam and for the vertical
        displacement of the elastic curve.
        """
        diagram = LineChart()
        diagram.add_xy_data(
            label='elastic curve',
            x1_values=self.x.to(units[0]).magnitude,
            y1_values=self.y.to(units[1]).magnitude,
            style_props={'drawstyle': 'steps-post'}
        )
        diagram.x1.add_title(f'longitudinal axis, {units[0]}')
        diagram.y1.add_title(f'vertical displacement, {units[1]}')
        return diagram

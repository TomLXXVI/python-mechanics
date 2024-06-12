from __future__ import annotations
import warnings
import numpy as np
import sympy as sp
from ._position import *
from ._units import Units


class Force:
    num_decimals: int = 3

    def __init__(
        self,
        action_point: Position = ORIGIN,
        magnitude: Quantity | str = Q_(0.0, 'N'),
        theta: Quantity | str = Q_(0.0, 'deg'),
        gamma: Quantity | str = Q_(0.0, 'deg'),
        name: str = ''
    ) -> None:
        self.action_point = set_position(action_point)
        if isinstance(magnitude, Quantity) and magnitude.m < 0:
            warnings.warn(
                "Negative magnitudes aren't allowed. The absolute value is taken.",
                category=RuntimeWarning
            )
            magnitude = np.abs(magnitude)
        self.magnitude = magnitude
        self.theta = theta
        self.gamma = gamma
        if name:
            self.name = name
        elif isinstance(magnitude, str):
            self.name = magnitude
        else:
            self.name = 'not specified'
            c1 = isinstance(self.magnitude, Quantity)
            c2 = isinstance(self.theta, Quantity)
            c3 = isinstance(self.gamma, Quantity)
            if not (c1 and c2 and c3):
                warnings.warn(
                    "An undetermined `Force` object should have a meaningful name.",
                    category=RuntimeWarning
                )
        self.__set_units()
        self._magnitude = self.__get_magnitude(magnitude)
        self._theta = self.__get_angle(theta)
        self._gamma = self.__get_angle(gamma)
        F_xy = self.__get_component_xy()
        self._F_x = self.__get_component_x(F_xy)
        self._F_y = self.__get_component_y(F_xy)
        self._F_z = self.__get_component_z()

    @classmethod
    def create_from_components(
        cls,
        action_point: Position = ORIGIN,
        F_x: Quantity | sp.Expr = Q_(0, 'N * m'),
        F_y: Quantity | sp.Expr = Q_(0, 'N * m'),
        F_z: Quantity | sp.Expr = Q_(0, 'N * m'),
        name: str = ''
    ) -> Force:
        c1 = isinstance(F_x, Quantity)
        c2 = isinstance(F_y, Quantity)
        c3 = isinstance(F_z, Quantity)
        if all([c1, c2, c3]):
            F = np.sqrt(F_x ** 2 + F_y ** 2 + F_z ** 2)
            F_xy = np.sqrt(F_x ** 2 + F_y ** 2)
            gamma = Q_(np.arctan2(F_z, F_xy), 'rad')
            theta = Q_(np.arctan2(F_y, F_x), 'rad')
            force = cls(action_point, F, theta, gamma, name)
        else:
            force = cls(action_point, name=name)
        force._F_x = F_x
        force._F_y = F_y
        force._F_z = F_z
        return force

    @property
    def components(self) -> tuple[Quantity | sp.Expr, ...]:
        return self._F_x, self._F_y, self._F_z

    @property
    def x(self) -> Quantity | sp.Expr:
        return self._F_x

    @property
    def y(self) -> Quantity | sp.Expr:
        return self._F_y

    @property
    def z(self) -> Quantity | sp.Expr:
        return self._F_z

    def moment(self, ref_point: Position | None = None) -> Moment:
        x_act = self.action_point[0].to(self.__u_pos)
        y_act = self.action_point[1].to(self.__u_pos)
        z_act = self.action_point[2].to(self.__u_pos)
        if ref_point is not None:
            ref_point = set_position(ref_point)
            x_ref = ref_point[0].to(self.__u_pos)
            y_ref = ref_point[1].to(self.__u_pos)
            z_ref = ref_point[2].to(self.__u_pos)
            x = x_act - x_ref  # distance from force action point to ref. point
            y = y_act - y_ref
            z = z_act - z_ref
        else:
            ref_point = self.action_point
            x = x_act
            y = y_act
            z = z_act
        F_x, F_y, F_z = self.components
        M_x = self.__get_moment_x(z, y, F_y, F_z)
        M_y = self.__get_moment_y(x, z, F_x, F_z)
        M_z = self.__get_moment_z(x, y, F_x, F_y)
        moment = Moment.create_from_components(ref_point, M_x, M_y, M_z)
        return moment

    @classmethod
    def reverse(cls, force: Force) -> Force:
        if not (isinstance(force.theta, str) and isinstance(force.gamma, str)):
            reversed_force = cls(
                action_point=force.action_point,
                magnitude=force.magnitude,
                theta=force.theta + Q_(180, 'deg'),
                gamma=force.gamma if force.gamma.m == 0.0 else force.gamma + Q_(180, 'deg'),
                name=f"-{force.name}"
            )
            return reversed_force
        else:
            raise ValueError("cannot reverse a symbolic force")

    def __repr__(self) -> str:
        s = "<"
        l = [('x', self._F_x), ('y', self._F_y), ('z', self._F_z)]
        for i, tup in enumerate(l):
            if isinstance(tup[1], Quantity):
                s += f"{tup[0]}: {tup[1]:~P.{self.num_decimals}f}"
            else:
                s += f"{tup[0]}: {str(tup[1])}"
            if i == 2:
                s += ">"
            else:
                s += "; "
        return s

    def __set_units(self) -> None:
        self.__u_pos = Units.u_pos
        self.__u_force = Units.u_force

    def __get_magnitude(self, magnitude) -> float | sp.Symbol:
        if isinstance(magnitude, Quantity):
            F = magnitude.to(self.__u_force).m
            return F
        elif isinstance(magnitude, str):
            F = sp.Symbol(magnitude)
            return F
        else:
            raise ValueError("neither a `Quantity`, nor `str`")

    @staticmethod
    def __get_angle(a) -> float | sp.Symbol:
        if isinstance(a, Quantity):
            a = a.to('rad').m
        elif isinstance(a, str):
            a = sp.Symbol(a)
        else:
            raise ValueError("neither a `Quantity`, nor `str`")
        return a

    def __get_component_xy(self) -> float | sp.Expr:
        F_xy = self._magnitude * sp.cos(self._gamma)
        if isinstance(F_xy, sp.Float) and F_xy < 0.0:
            F_xy = np.abs(F_xy)
        elif isinstance(self._gamma, float):
            cos_gamma = round(np.cos(self._gamma), 12)
            if cos_gamma == 0:
                F_xy = 0.0
            elif cos_gamma == 1:
                F_xy = self._magnitude
            elif cos_gamma == -1:
                F_xy = -self._magnitude
            else:
                F_xy = self._magnitude * cos_gamma
        else:
            F_xy = sp.Symbol(f"{self.name}.xy")
        return F_xy

    def __get_component_x(self, F_xy: float | sp.Expr) -> Quantity | sp.Expr:
        F_x = F_xy * sp.cos(self._theta)
        if isinstance(F_x, sp.Number):
            F_x = Q_(round(float(F_x), 12), self.__u_force)
        else:
            if isinstance(self._theta, float):
                cos_theta = round(np.cos(self._theta), 12)
                if cos_theta == 0:
                    F_x = Q_(0.0, self.__u_force)
                elif cos_theta == 1:
                    F_x = F_xy
                elif cos_theta == -1:
                    F_x = -F_xy
                else:
                    pass
            else:
                F_x = sp.Symbol(f"{self.name}.x")
        return F_x

    def __get_component_y(self, F_xy: float | sp.Expr) -> Quantity | sp.Expr:
        F_y = F_xy * sp.sin(self._theta)
        if isinstance(F_y, sp.Number):
            F_y = Q_(round(float(F_y), 12), self.__u_force)
        else:
            if isinstance(self._theta, float):
                sin_theta = round(np.sin(self._theta), 12)
                if sin_theta == 0:
                    F_y = Q_(0.0, self.__u_force)
                elif sin_theta == 1:
                    F_y = F_xy
                elif sin_theta == -1:
                    F_y = -F_xy
                else:
                    pass
            else:
                F_y = sp.Symbol(f"{self.name}.y")
        return F_y

    def __get_component_z(self) -> Quantity | sp.Expr:
        F_z = self._magnitude * sp.sin(self._gamma)
        if isinstance(F_z, sp.Number):
            F_z = Q_(round(float(F_z), 12), self.__u_force)
        else:
            if isinstance(self._gamma, float):
                sin_gamma = round(np.sin(self._gamma), 12)
                if sin_gamma == 0:
                    F_z = Q_(0.0, self.__u_force)
                elif sin_gamma == 1:
                    F_z = self._magnitude
                elif sin_gamma == -1:
                    F_z = -self._magnitude
                else:
                    F_z = self._magnitude * sin_gamma
            else:
                F_z = sp.Symbol(f"{self.name}.z")
        return F_z

    def __get_moment_x(
        self,
        z: Quantity,
        y: Quantity,
        F_y: Quantity | sp.Expr,
        F_z: Quantity | sp.Expr
    ) -> Quantity | sp.Expr:
        if not (isinstance(F_y, Quantity) and isinstance(F_z, Quantity)):
            z = z.to(self.__u_pos).m
            y = y.to(self.__u_pos).m
            if isinstance(F_y, Quantity):
                F_y = F_y.to(self.__u_force).m
            if isinstance(F_z, Quantity):
                F_z = F_z.to(self.__u_force).m
        M_x = -z * F_y + y * F_z
        return M_x

    def __get_moment_y(
        self,
        x: Quantity,
        z: Quantity,
        F_x: Quantity | sp.Expr,
        F_z: Quantity | sp.Expr
    ) -> Quantity | sp.Expr:
        if not (isinstance(F_x, Quantity) and isinstance(F_z, Quantity)):
            x = x.to(self.__u_pos).m
            z = z.to(self.__u_pos).m
            if isinstance(F_x, Quantity):
                F_x = F_x.to(self.__u_force).m
            if isinstance(F_z, Quantity):
                F_z = F_z.to(self.__u_force).m
        M_y = z * F_x - x * F_z
        return M_y

    def __get_moment_z(
        self,
        x: Quantity,
        y: Quantity,
        F_x: Quantity | sp.Expr,
        F_y: Quantity | sp.Expr
    ) -> Quantity | sp.Expr:
        if not (isinstance(F_x, Quantity) and isinstance(F_y, Quantity)):
            x = x.to(self.__u_pos).m
            y = y.to(self.__u_pos).m
            if isinstance(F_x, Quantity):
                F_x = F_x.to(self.__u_force).m
            if isinstance(F_y, Quantity):
                F_y = F_y.to(self.__u_force).m
        M_z = -y * F_x + x * F_y
        return M_z


class Moment:
    num_decimals: int = 3

    def __init__(
        self,
        action_point: Position = ORIGIN,
        magnitude: Quantity | str = Q_(0.0, 'N * m'),
        theta: Quantity | str = Q_(0.0, 'deg'),
        gamma: Quantity | str = Q_(0.0, 'deg'),
        name: str = ''
    ) -> None:
        self.action_point = set_position(action_point)
        if isinstance(magnitude, Quantity) and magnitude.m < 0:
            warnings.warn(
                "Negative magnitudes aren't allowed. The absolute value is taken.",
                category=RuntimeWarning
            )
            magnitude = np.abs(magnitude)
        self.magnitude = magnitude
        self.theta = theta
        self.gamma = gamma
        if name:
            self.name = name
        elif isinstance(magnitude, str):
            self.name = magnitude
        else:
            c1 = isinstance(self.magnitude, Quantity)
            c2 = isinstance(self.theta, Quantity)
            c3 = isinstance(self.gamma, Quantity)
            if not (c1 and c2 and c3):
                warnings.warn(
                    "An undetermined `Moment` object should have a meaningful name.",
                    category=RuntimeWarning
                )
        self.__set_units()
        self._magnitude = self.__get_magnitude(magnitude)
        self._theta = self.__get_angle(theta)
        self._gamma = self.__get_angle(gamma)
        M_xy = self.__get_component_xy()
        self._M_x = self.__get_component_x(M_xy)
        self._M_y = self.__get_component_y(M_xy)
        self._M_z = self.__get_component_z()

    @classmethod
    def create_from_components(
        cls,
        action_point: Position = ORIGIN,
        M_x: Quantity | sp.Expr = Q_(0, 'N * m'),
        M_y: Quantity | sp.Expr = Q_(0, 'N * m'),
        M_z: Quantity | sp.Expr = Q_(0, 'N * m'),
        name: str = ''
    ) -> Moment:
        c1 = isinstance(M_x, Quantity)
        c2 = isinstance(M_y, Quantity)
        c3 = isinstance(M_z, Quantity)
        if all([c1, c2, c3]):
            M = np.sqrt(M_x ** 2 + M_y ** 2 + M_z ** 2)
            M_xy = np.sqrt(M_x ** 2 + M_y ** 2)
            gamma = Q_(np.arctan2(M_z, M_xy), 'rad')
            theta = Q_(np.arctan2(M_y, M_x), 'rad')
            moment = cls(action_point, M, theta, gamma, name)
        else:
            moment = cls(action_point, name=name)
        moment._M_x = M_x
        moment._M_y = M_y
        moment._M_z = M_z
        return moment

    @property
    def components(self) -> tuple[Quantity | sp.Expr, ...]:
        return self._M_x, self._M_y, self._M_z

    @property
    def x(self) -> Quantity | sp.Expr:
        return self._M_x

    @property
    def y(self) -> Quantity | sp.Expr:
        return self._M_y

    @property
    def z(self) -> Quantity | sp.Expr:
        return self._M_z

    @classmethod
    def reverse(cls, moment: Moment) -> Moment:
        if not (isinstance(moment.theta, str) and isinstance(moment.gamma, str)):
            reversed_moment = Moment(
                action_point=moment.action_point,
                magnitude=moment.magnitude,
                theta=moment.theta + Q_(180, 'deg'),
                gamma=moment.gamma if moment.gamma.m == 0.0 else moment.gamma + Q_(180, 'deg'),
                name=f"-{moment.name}"
            )
            return reversed_moment
        else:
            raise ValueError("cannot reverse a symbolic moment")

    def __repr__(self) -> str:
        s = "<"
        l = [('x', self._M_x), ('y', self._M_y), ('z', self._M_z)]
        for i, tup in enumerate(l):
            if isinstance(tup[1], Quantity):
                s += f"{tup[0]}: {tup[1]:~P.{self.num_decimals}f}"
            else:
                s += f"{tup[0]}: {str(tup[1])}"
            if i == 2:
                s += ">"
            else:
                s += "; "
        return s

    def __set_units(self) -> None:
        self.__u_pos = Units.u_pos
        self.__u_moment = Units.u_moment

    def __get_magnitude(self, magnitude) -> float | sp.Symbol:
        if isinstance(magnitude, Quantity):
            M = magnitude.to(self.__u_moment).m
        elif isinstance(magnitude, str):
            M = sp.Symbol(magnitude)
        else:
            raise ValueError("neither a `Quantity`, nor `str`")
        return M

    @staticmethod
    def __get_angle(a) -> float | sp.Symbol:
        if isinstance(a, Quantity):
            a = a.to('rad').m
        elif isinstance(a, str):
            a = sp.Symbol(a)
        else:
            raise ValueError("neither a `Quantity`, nor `str`")
        return a

    def __get_component_xy(self) -> float | sp.Expr:
        M_xy = self._magnitude * sp.cos(self._gamma)
        if isinstance(M_xy, float) and M_xy < 0.0:
            M_xy = np.abs(M_xy)
        elif isinstance(self._gamma, float):
            cos_gamma = round(np.cos(self._gamma), 12)
            if cos_gamma == 0:
                M_xy = 0.0
            elif cos_gamma == 1:
                M_xy = self._magnitude
            else:
                M_xy = self._magnitude * cos_gamma
        else:
            M_xy = sp.Symbol(self.name)
        return M_xy

    def __get_component_x(self, M_xy: float | sp.Expr) -> Quantity | sp.Expr:
        M_x = M_xy * sp.cos(self._theta)
        if isinstance(M_x, sp.Number):
            M_x = Q_(float(M_x), self.__u_moment)
        else:
            if isinstance(self._theta, float):
                cos_theta = round(np.cos(self._theta), 12)
                if cos_theta == 0:
                    M_x = Q_(0.0, self.__u_moment)
                else:
                    pass
            else:
                M_x = sp.Symbol(f"{str(M_xy)}.x")
        return M_x

    def __get_component_y(self, M_xy: float | sp.Expr) -> Quantity | sp.Expr:
        M_y = M_xy * sp.sin(self._theta)
        if isinstance(M_y, sp.Number):
            M_y = Q_(float(M_y), self.__u_moment)
        else:
            if isinstance(self._theta, float):
                sin_theta = round(np.sin(self._theta), 12)
                if sin_theta == 0:
                    M_y = Q_(0.0, self.__u_moment)
                else:
                    pass
            else:
                M_y = sp.Symbol(f"{str(M_xy)}.y")
        return M_y

    def __get_component_z(self) -> Quantity | sp.Expr:
        M_z = self._magnitude * sp.sin(self._gamma)
        if isinstance(M_z, sp.Number):
            M_z = Q_(float(M_z), self.__u_moment)
        else:
            if isinstance(self._gamma, float):
                sin_gamma = round(np.sin(self._gamma), 12)
                if sin_gamma == 0:
                    M_z = Q_(0.0, self.__u_moment)
                else:
                    pass
            else:
                M_z = sp.Symbol(f"{str(self._magnitude)}.z")
        return M_z


class Torque(Moment):
    pass

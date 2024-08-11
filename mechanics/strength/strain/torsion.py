from typing import Callable
from scipy.integrate import quad
from mechanics import Quantity


Q_ = Quantity


class TorsionalStrain:
    """
    Determine the angle of twist along a circular shaft due to torsional
    loading.
    """
    def __init__(
        self,
        T: Quantity | Callable[[Quantity], Quantity],
        J: Quantity | Callable[[Quantity], Quantity],
        G: Quantity | Callable[[Quantity], Quantity]
    ) -> None:
        """Creates a `TorsionalStrain` object.

        Parameters
        ----------
        T:
            Resultant internal torque, being a `Quantity` object if the torque
            is constant throughout the length of the shaft. Otherwise, it must
            be a function that takes the position `x` along the longitudinal
            axis of the shaft and returns the normal force at that position.
        J:
            The shaft's polar moment of inertia, being a `Quantity` object if
            the circular shaft has a constant cross-section. Otherwise, it must
            be a function that takes the position `x` along the longitudinal
            axis of the shaft and returns the polar moment of inertia at that
            position.
        G:
            The shear modulus of elasticity for the material, being a `Quantity`
            if the modulus is constant throughout the length of the shaft.
            Otherwise, it is a function that takes the position `x` along the
            longitudinal axis of the shaft and returns the modulus at that
            position.
        """
        self.T = T
        self.J = J
        self.G = G
        if isinstance(self.T, Quantity):
            self._T = lambda x: self.T.to('N * m').magnitude
        else:
            self._T = self.__create_T_function()
        if isinstance(self.J, Quantity):
            self._J = lambda x: self.J.to('m**4').magnitude
        else:
            self._J = self.__create_J_function()
        if isinstance(self.G, Quantity):
            self._G = lambda x: self.G.to('Pa').magnitude
        else:
            self._G = self.__create_G_function()

    def __create_T_function(self) -> Callable[[float], float]:
        def T(x: float) -> float:
            x = Q_(x, 'm')
            T = self.T(x)
            T = T.to('N * m').magnitude
            return T
        return T

    def __create_J_function(self) -> Callable[[float], float]:
        def J(x: float) -> float:
            x = Q_(x, 'm')
            J = self.J(x)
            J = J.to('m**4').magnitude
            return J
        return J

    def __create_G_function(self) -> Callable[[float], float]:
        def G(x: float) -> float:
            x = Q_(x, 'm')
            G = self.G(x)
            G = G.to('Pa').magnitude
            return G
        return G

    def __angle_of_twist_with_constants(self, x1: Quantity, x2: Quantity) -> Quantity:
        L = x2 - x1
        fi = self.T * L / (self.J * self.G)
        return fi.to('rad')

    def __angle_of_twist_by_integration(self, x1: Quantity, x2: Quantity) -> Quantity:
        x1 = x1.to('m').magnitude
        x2 = x2.to('m').magnitude

        def f(x):
            f = self._T(x) / (self._J(x) * self._G(x))
            return f

        fi = quad(f, x1, x2)[0]
        fi = Q_(fi, 'rad')
        return fi

    def angle_of_twist(self, x1: Quantity, x2: Quantity) -> Quantity:
        """Returns the angle of twist of the circular shaft between its
        longitudinal positions `x1` and `x2`, where `x2` must be greater than
        `x1`.
        """
        if (
            isinstance(self.T, Quantity)
            and isinstance(self.J, Quantity)
            and isinstance(self.G, Quantity)
        ):
            fi = self.__angle_of_twist_with_constants(x1, x2)
        else:
            fi = self.__angle_of_twist_by_integration(x1, x2)
        return fi

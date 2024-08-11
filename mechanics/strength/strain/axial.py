from typing import Callable
from scipy.integrate import quad
from mechanics import Quantity


Q_ = Quantity


class AxialStrain:
    """
    Determine the axial strain or elongation of a beam due to axial loading.
    """
    def __init__(
        self,
        N: Quantity | Callable[[Quantity], Quantity],
        A: Quantity | Callable[[Quantity], Quantity],
        E: Quantity | Callable[[Quantity], Quantity]
    ) -> None:
        """Creates an `AxialStrain` object.

        Parameters
        ----------
        N:
            Resultant internal normal force, being a `Quantity` object if the
            the normal force is constant throughout the length of the beam.
            Otherwise, it must be a function that takes the position `x` along
            the longitudinal axis of the beam and returns the normal force at
            that position.
        A:  Area of the beam's cross-section, being a `Quantity` object if the
            area of the cross-section is constant throughout the length of the
            beam. Otherwise it must be a function that takes the position `x`
            along the longitudinal axis of the beam and returns the area at that
            position.
        E:
            Modulus of elasticity of the material, being a `Quantity` object if
            the modulus is constant throughout the length of the beam. Otherwise
            it must be a function that takes the position `x` along the
            longitudinal axis of the beam and returns the modulus at that
            position.
        """
        self.N = N
        self.A = A
        self.E = E
        if isinstance(self.N, Quantity):
            self._N = lambda x: self.N.to('N').magnitude
        else:
            self._N = self.__create_N_function()
        if isinstance(self.A, Quantity):
            self._A = lambda x: self.A.to('m**2').magnitude
        else:
            self._A = self.__create_A_function()
        if isinstance(self.E, Quantity):
            self._E = lambda x: self.E.to('Pa').magnitude
        else:
            self._E = self.__create_E_function()

    def __create_N_function(self) -> Callable[[float], float]:
        def N(x: float) -> float:
            x = Q_(x, 'm')
            N = self.N(x)
            N = N.to('N').magnitude
            return N
        return N

    def __create_A_function(self) -> Callable[[float], float]:
        def A(x: float) -> float:
            x = Q_(x, 'm')
            A = self.A(x)
            A = A.to('m**2').magnitude
            return A
        return A

    def __create_E_function(self) -> Callable[[float], float]:
        def E(x: float) -> float:
            x = Q_(x, 'm')
            E = self.E(x)
            E = E.to('Pa').magnitude
            return E
        return E

    def __elongation_with_constants(self, x1: Quantity, x2: Quantity) -> Quantity:
        L = x2 - x1
        delta = self.N * L / (self.A * self.E)
        return delta

    def __elongation_by_integration(self, x1: Quantity, x2: Quantity) -> Quantity:
        x1 = x1.to('m').magnitude
        x2 = x2.to('m').magnitude

        def f(x):
            f = self._N(x) / (self._A(x) * self._E(x))
            return f

        delta = quad(f, x1, x2)[0]
        delta = Q_(delta, 'm')
        return delta

    def elongation(self, x1: Quantity, x2: Quantity) -> Quantity:
        """Returns the elongation of the beam between its longitudinal positions
        `x1` and `x2`, where `x2` must be greater than `x1`.
        """
        if (
            isinstance(self.N, Quantity)
            and isinstance(self.A, Quantity)
            and isinstance(self.E, Quantity)
        ):
            delta = self.__elongation_with_constants(x1, x2)
        else:
            delta = self.__elongation_by_integration(x1, x2)
        return delta

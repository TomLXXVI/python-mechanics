from mechanics import Quantity


class Stress:

    def __init__(
        self,
        F_mag: Quantity,
        stress_allow: Quantity | None = None,
        stress_fail: Quantity | None = None,
        safety_factor: float | None = None,
        A: Quantity | None = None
    ) -> None:
        self.F_mag = F_mag
        if A is not None:
            # analysis
            self._A = A
            self._stress_avg = self.__calculate_avg_stress(self.F_mag, A)
        else:
            # design
            self._A = self.__calculate_area(self.F_mag, stress_allow, stress_fail, safety_factor)
            self._stress_avg = self.__calculate_avg_stress(self.F_mag, self._A)

    def design(self) -> Quantity:
        return self._A

    @property
    def average(self) -> Quantity:
        return self._stress_avg

    @staticmethod
    def __calculate_avg_stress(F_mag: Quantity, A: Quantity) -> Quantity:
        return F_mag / A

    @staticmethod
    def __calculate_area(
        F_mag: Quantity,
        stress_allow: Quantity | None,
        stress_fail: Quantity | None,
        safety_factor: float | None
    ) -> Quantity:
        if (safety_factor is not None) and (stress_fail is not None):
            stress_allow = stress_fail / safety_factor
        A = F_mag / stress_allow
        return A

from mechanics import Quantity
from ._stress import Stress


class ShearStress(Stress):

    def __init__(
        self,
        V: Quantity,
        tau_allow: Quantity | None = None,
        tau_fail: Quantity | None = None,
        safety_factor: float | None = None,
        A: Quantity | None = None
    ) -> None:
        super().__init__(V, tau_allow, tau_fail, safety_factor, A)

    @property
    def tau_avg(self) -> Quantity:
        return self._stress_avg

from mechanics import Quantity
from ._stress import Stress


class NormalStress(Stress):

    def __init__(
        self,
        N: Quantity,
        sigma_allow: Quantity | None = None,
        sigma_fail: Quantity | None = None,
        safety_factor: float | None = None,
        A: Quantity | None = None
    ) -> None:
        super().__init__(N, sigma_allow, sigma_fail, safety_factor, A)

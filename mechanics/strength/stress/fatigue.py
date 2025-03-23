"""CALCULATION OF FATIGUE STRENGTH.

References
----------
[1] Matek, W., Muhs, D., Wittel, H., Becker, M. Roloff/Matek Machine-onderdelen
    (1996). Academic Service.
"""
import numpy as np
from mechanics import Quantity

Q_ = Quantity


class Fatigue:

    def __init__(
        self,
        ts: Quantity,
        ys: Quantity,
        K1: float,
        K2: float,
        diameter: Quantity | None = None,
        r_notch: Quantity | None = None,
        f_notch: float | None = None,
        Rz: Quantity = Q_(1, 'µm'),
        is_torsional_load: bool = False
    ) -> None:
        """Creates a `FatigueStrength` object.

        Parameters
        ----------
        ts:
            Static tensile strength of material.
        ys:
            Static yield strength of material.
        K1:
            Ratio of the alternating fatigue strength to the static tensile
            strength of the material. This value depends on the material and the
            type of load (tensile, bending, or torsional load), e.g. see ref
            [1], table 3.1.
        K2:
            Ratio of the dynamic yield strength to the static yield strength of
            the material. This value depends on the material and the type of
            load (tensile, bending, or torsional load), e.g. see ref [1],
            table 3.1.

        Other Parameters
        ----------------
        These parameters can be used to calculate the design strength of a round
        machine part.
        diameter: optional
            Cross-section diameter in case of a round spindle or axle.
        r_notch: optional
            Notch radius.
        f_notch: optional
            Notch shape factor, e.g. see ref [1], table 3.7.
        Rz: optional
            Surface roughness.
        is_torsional_load: optional
            Indicates if design strength needs to be determined for a torsional
            load.
        """
        self.ts = ts.to('N / mm**2')
        self.ys = ys.to('N / mm**2')
        self.K1 = K1
        self.K2 = K2
        self.d = diameter
        self.r_notch = r_notch
        self.f_notch = f_notch
        self.Rz = Rz.to('µm')
        self.is_torsional_load = is_torsional_load
        self.as_dyn = self.K1 * self.ts  # dynamic alternating strength
        self.ys_dyn = self.K2 * self.ys  # dynamic yield strength

    def strength(
        self,
        stress_ratio: float | None = None,
        avg_stress: Quantity | None = None
    ) -> Quantity:
        """Returns the (approximate) fatigue strength of the material which
        depends on characteristics of the applied dynamic load.

        Either parameter `stress_ratio`, or parameter `avg_stress` needs to be
        specified. If both are specified, parameter `stress_ratio` takes
        precedence over parameter `avg_stress`.

        Parameters
        ----------
        stress_ratio:
            Ratio of the minimum to maximum stress level of the dynamic load.
            The stress ratio has a value between -1 and 1.
            * 1 :  the load is static and minimum and maximum stress are equal;
            * 0 :  the load is a swelling load with the minimum stress being
                   equal to zero;
            * -1 : the load is an alternating load for which the minimum stress
                   is the negative of the maximum stress.
        avg_stress:
            Average value around which the dynamic (shear) stress fluctuates.
        """
        if stress_ratio is not None:
            n = self.as_dyn
            d = 1 - ((1 + stress_ratio) * (1 - self.K1)) / (2 - self.K1)
            fs = min(n / d, self.ys_dyn)
            return fs
        elif avg_stress is not None:
            k = (1 - self.K1) / (1 - self.K1 / 2)
            fs = k * avg_stress + self.as_dyn
            fs = min(fs, self.ys_dyn)
            return fs

    def amplitude_strength(
        self,
        stress_ratio: float
    ) -> Quantity:
        """Returns the amplitude strength of the material which depends on
        characteristics of the applied dynamic load.

        Parameters
        ----------
        stress_ratio:
            Ratio of the minimum to maximum stress level of the dynamic load.
            The stress ratio has a value between -1 and 1.
            * 1 :  the load is static and minimum and maximum stress are equal;
            * 0 :  the load is a swelling load with the minimum stress being
                   equal to zero;
            * -1 : the load is an alternating load for which the minimum stress
                   is the negative of the maximum stress.
        """
        fs = self.strength(stress_ratio)
        as_ = fs / 2 * (1 - stress_ratio)
        return as_

    def _surface_condition_factor(self) -> float:
        Rz = self.Rz.to('µm').m
        ts = self.ts.to('N / mm**2').m
        b = 1 - 0.22 * np.log10(Rz) * (np.log10(ts / 20) - 1)
        if self.is_torsional_load:
            b = 0.575 * b + 0.425
        return b

    def _size_factor(self) -> float:
        if self.d is not None:
            d = self.d.to('mm').m
            k_t = 1 - 0.25 * (np.log(d / 7.5) / np.log(20))
            # k_t should be 1 for general and tempered structural steel
            k_g = 1 - 0.2 * (np.log(d / 7.5) / np.log(20))
            # k_g should be 1 for tensile/compressive loading
            if self.f_notch is not None:
                k_a = 1 - 0.2 * np.log(self.f_notch) * (np.log(d / 7.5) / np.log(20))
            else:
                k_a = 1.0
        else:
            k_t = 1.0
            k_g = 1.0
            k_a = 1.0
        b = k_t * k_g * k_a
        return b

    def _dynamic_notch_factor(self) -> float:
        if self.r_notch is not None:
            r = self.r_notch.to('mm').m
            ys = self.ys.to('N / mm**2').m
            ts = self.ts.to('N / mm**2').m
            eta_k = 1 / (1 + (8 / r) * (1 - ys / ts) ** 3)
            beta_k = 1 + eta_k * (self.f_notch - 1)
        else:
            beta_k = 1
        return beta_k

    def design_strength(
        self,
        stress_ratio: float | None = None,
        avg_stress: Quantity | None = None
    ) -> Quantity:
        """Returns the (approximate) fatigue strength of the machine part
        depending on the characteristics of the dynamic load, while also taking
        strength reducing factors of the machine part into account, namely its
        surface roughness, its size (diameter), or the presence of a notch.

        Either parameter `stress_ratio`, or parameter `avg_stress` needs to be
        specified. If both are specified, parameter `stress_ratio` takes
        precedence over parameter `avg_stress`.

        Parameters
        ----------
        stress_ratio:
            Ratio of the minimum to maximum stress level of the dynamic load.
            The stress ratio has a value between -1 and 1.
            * 1 :  the load is static and minimum and maximum stress are equal;
            * 0 :  the load is a swelling load with the minimum stress being
                   equal to zero;
            * -1 : the load is an alternating load for which the minimum stress
                   is the negative of the maximum stress.
        avg_stress:
            Average value around which the dynamic (shear) stress fluctuates.
        """
        fs = self.strength(stress_ratio, avg_stress)
        b1 = self._surface_condition_factor()
        b2 = self._size_factor()
        b3 = self._dynamic_notch_factor()
        fs_des = fs * (b1 * b2 / b3)
        return fs_des

    def design_amplitude_strength(
        self,
        stress_ratio: float
    ) -> Quantity:
        """Returns the amplitude strength of the machine part depending on the
        characteristics of the dynamic load, while also taking strength reducing
        factors of the machine part into account, namely its surface roughness,
        its size (diameter), or the presence of a notch.

        Parameters
        ----------
        stress_ratio:
            Ratio of the minimum to maximum stress level of the dynamic load.
            The stress ratio has a value between -1 and 1.
            * 1 :  the load is static and minimum and maximum stress are equal;
            * 0 :  the load is a swelling load with the minimum stress being
                   equal to zero;
            * -1 : the load is an alternating load for which the minimum stress
                   is the negative of the maximum stress.
        """
        as_ = self.amplitude_strength(stress_ratio)
        b1 = self._surface_condition_factor()
        b2 = self._size_factor()
        b3 = self._dynamic_notch_factor()
        as_des = as_ * (b1 * b2 / b3)
        return as_des

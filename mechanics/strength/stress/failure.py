import numpy as np
from mechanics import Quantity
from mechanics.charts import LineChart


class VonMisesYieldCriterion:
    """aka maximum distortion energy theory. Applies to materials that behave in
    a ductile manner.
    """
    def __init__(self, sigma_yield: Quantity):
        """Creates a `VonMisesYieldCriterion` object.

        Parameters
        ----------
        sigma_yield:
            Yield strength of the material.
        """
        self.sigma_yield = sigma_yield.to('N / mm**2')

    def draw_chart(
        self,
        sigma1: Quantity,
        sigma2: Quantity,
        units: str = 'N / mm**2'
    ) -> LineChart:
        """Returns a `LineChart` object with a graphical representation of
        the von Mises yield criterion.

        Parameters
        ----------
        sigma1:
            First principal normal stress.
        sigma2:
            Second principal normal stress.
        units:
            Units in which the stresses must be expressed.

        The point (`sigma1`, `sigma2`) is shown by a red dot on the chart. If
        the dot lies within the boundaries of the von Misses ellipse, no failure
        of the material should happen. Otherwise, failure (yielding) is
        probable.
        """
        arr_sigma1, arr_sigma2 = self.__calculate_curve(units)
        chart = LineChart()
        chart.axes.set_aspect(1)
        chart.axes.axhline(0, color='black', lw=1.0)
        chart.axes.axvline(0, color='black', lw=1.0)
        chart.add_xy_data(
            label='curve',
            x1_values=arr_sigma1,
            y1_values=arr_sigma2,
        )
        chart.add_xy_data(
            label='point',
            x1_values=[sigma1.to(units).m],
            y1_values=[sigma2.to(units).m],
            style_props={'marker': 'o', 'color': 'red', 'linestyle': 'none'}
        )
        chart.x1.add_title(f'sigma 1, {units}')
        chart.y1.add_title(f'sigma 2, {units}')
        return chart

    @staticmethod
    def von_mises_stress(sigma1: Quantity, sigma2: Quantity) -> Quantity:
        """Returns the von Mises stress, given the principal normal stresses
        `sigma1` and `sigma2`.
        """
        sigma1 = sigma1.to('N / mm**2')
        sigma2 = sigma2.to('N / mm**2')
        sigma_vm = np.sqrt(sigma1 ** 2 - sigma1 * sigma2 + sigma2 ** 2)
        return sigma_vm

    def check(self, sigma1: Quantity, sigma2: Quantity) -> bool:
        """Checks if the von Mises yield criterion (aka maximum distortion
        energy criterion) is satisfied. For this, the calculated von Mises
        stress must be smaller than the material's yield strength.
        """
        ys = self.sigma_yield.to('N / mm**2')
        sigma_vm = self.von_mises_stress(sigma1, sigma2)
        if sigma_vm < ys:
            return True
        return False

    @staticmethod
    def __solve_sigma2(sigma1: float, ys: float) -> tuple[tuple[float, float], ...] | None:
        # Solves the quadratic equation
        # `sigma2**2 - sigma1 * sigma2 + (sigma_1**2 - ys**2) = 0`.
        a = 1.0
        b = -sigma1
        c = sigma1 ** 2 - ys ** 2
        D = b ** 2 - 4 * a * c
        if D >= 0:
            sigma2_1 = (-b + np.sqrt(D)) / (2 * a)
            sigma2_2 = (-b - np.sqrt(D)) / (2 * a)
            return (sigma1, sigma2_1), (sigma1, sigma2_2)
        return None

    def __calculate_curve(self, units: str) -> tuple[list[float], list[float]]:
        ys = self.sigma_yield.to(units)._magnitude
        # Select a range of values on the sigma1 axis:
        arr_sigma1 = np.arange(-1.20 * ys, 1.20 * ys, 1)
        # Calculate for each sigma1 in the range the corresponding sigma2.
        # The equation is quadratic, so for each sigma1 there are 2 solutions
        # for sigma2. We store the first and the second solution in two separate
        # lists.
        solutions1, solutions2 = [], []
        for sigma1 in arr_sigma1:
            sol = self.__solve_sigma2(sigma1, ys)
            # __solve_sigma2 returns two tuples if solutions exists (discriminant
            # not negative). The first element of a tuple is sigma1, the second
            # element is one of the two solutions for sigma2.
            if sol is not None:
                solutions1.append(sol[0])
                solutions2.append(sol[1])
        # The tuples (sigma1, sigma2) in the two lists with solutions are
        # unzipped into two tuples, one with only the sigma1 values, and the
        # other with only the sigma2 values:
        arr1_sigma1, arr1_sigma2 = zip(*solutions1)
        arr2_sigma1, arr2_sigma2 = zip(*solutions2)
        # The tuples are converted to lists and we add the second list with
        # sigma1/sigma2 values in reversed order to the first list with sigma1/
        # sigma2 values.
        arr_sigma1 = list(arr1_sigma1) + list(arr2_sigma1)[::-1]
        arr_sigma2 = list(arr1_sigma2) + list(arr2_sigma2)[::-1]
        # Now we have the coordinates of consecutive points on the curve. To
        # close the curve, we add the first pair of coordinates also to the end
        # of the coordinate lists:
        arr_sigma1.append(arr_sigma1[0])
        arr_sigma2.append(arr_sigma2[0])
        return arr_sigma1, arr_sigma2


class MohrsFailureCriterion:
    """Applies to materials that behave in a brittle manner and for which the
    properties under tensile load and compressive loadings are different.
    """
    def __init__(
        self,
        sigma_ult_t: Quantity,
        sigma_ult_c: Quantity | None = None
    ) -> None:
        """Creates a `MohrsFailureCriterion` object.

        Parameters
        ----------
        sigma_ult_t:
            Ultimate tensile strength of the material with a positive sign.
        sigma_ult_c:
            Ultimate compressive strength of the material with a negative sign.
            If None, it is assumed that the compressive strength equals the
            tensile strength.
        """
        self.sigma_ult_t = sigma_ult_t.to('N / mm**2')
        if sigma_ult_c is not None:
            if sigma_ult_c.m > 0: sigma_ult_c = -sigma_ult_c
            self.sigma_ult_c = sigma_ult_c.to('N / mm**2')
        else:
            self.sigma_ult_c = -self.sigma_ult_t

    def check(self, sigma1: Quantity, sigma2: Quantity) -> bool:
        """Checks if Mohr's failure criterion is satisfied.

        Parameters
        ----------
        sigma1:
            Maximum normal stress.
        sigma2:
            Minimum normal stress.
        """
        sigma1 = sigma1.to('N / mm**2')
        sigma2 = sigma2.to('N / mm**2')
        if sigma1.m >= 0 and sigma2.m >= 0:
            if sigma1 < self.sigma_ult_t and sigma2 < self.sigma_ult_t:
                return True
            else:
                return False
        elif sigma1.m <= 0 and sigma2.m <= 0:
            if sigma1 > self.sigma_ult_c and sigma2 > self.sigma_ult_c:
                return True
            else:
                return False
        elif sigma1.m > 0 > sigma2.m:
            a = sigma1 / self.sigma_ult_t
            b = sigma2 / self.sigma_ult_c
            c = a + b
            if c < 1:
                return True
            else:
                return False
        elif sigma1.m < 0 < sigma2.m:
            a = sigma1 / self.sigma_ult_c
            b = sigma2 / self.sigma_ult_t
            c = a + b
            if c < 1:
                return True
            else:
                return False
        return False

    def draw_chart(
        self,
        sigma1: Quantity,
        sigma2: Quantity,
        units: str = 'N / mm**2'
    ) -> LineChart:
        sigma_ult_t = self.sigma_ult_t.to(units).magnitude
        sigma_ult_c = self.sigma_ult_c.to(units).magnitude
        arr_sigma1_neg = np.linspace(sigma_ult_c, 0, 50, endpoint=True)
        arr_sigma1_pos = np.linspace(0, sigma_ult_t, 50, endpoint=True)
        arr_sigma2_pos = (1 - arr_sigma1_neg / sigma_ult_c) * sigma_ult_t
        arr_sigma2_neg = (1 - arr_sigma1_pos / sigma_ult_t) * sigma_ult_c
        lines = [
            ('line1', [0, sigma_ult_t], [sigma_ult_t, sigma_ult_t]),
            ('line2', [sigma_ult_t, sigma_ult_t], [0, sigma_ult_t]),
            ('line3', [sigma_ult_c, sigma_ult_c], [0, sigma_ult_c]),
            ('line4', [sigma_ult_c, 0], [sigma_ult_c, sigma_ult_c]),
            ('line5', arr_sigma1_neg, arr_sigma2_pos),
            ('line6', arr_sigma1_pos, arr_sigma2_neg)
        ]
        chart = LineChart()
        chart.axes.set_aspect(1)
        chart.axes.axhline(0, color='black', lw=1.0)
        chart.axes.axvline(0, color='black', lw=1.0)
        for line in lines:
            chart.add_xy_data(
                label=line[0],
                x1_values=line[1],
                y1_values=line[2],
                style_props={'color': 'tab:blue'}
            )
        chart.add_xy_data(
            label='point',
            x1_values=[sigma1.to(units).m],
            y1_values=[sigma2.to(units).m],
            style_props={'marker': 'o', 'color': 'red', 'linestyle': 'none'}
        )
        chart.x1.add_title(f'sigma 1, {units}')
        chart.y1.add_title(f'sigma 2, {units}')
        return chart

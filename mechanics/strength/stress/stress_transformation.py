import numpy as np
from mechanics import Quantity
from mechanics.charts import LineChart


Q_ = Quantity


class PlaneStressTransformation:
    """
    Transformation of plane stress components acting on a surface element at a
    point in a cross-section into components that act on a rotated surface
    element at the same point.
    """
    def __init__(
        self,
        sigma_x: Quantity,
        sigma_y: Quantity,
        tau_xy: Quantity
    ) -> None:
        """Creates a `PlaneStressTransformation` object.

        Parameters
        ----------
        sigma_x:
            Normal stress acting on the +x-face of the original element.
            If the direction of the normal stress is opposite to the positive
            x-direction, it must have a negative value.
        sigma_y:
            Normal stress acting on the +y-face of the original element.
            If the direction of the normal stress is opposite to the positive
            y-direction, it must have a negative value.
        tau_xy:
            Shear stress acting in the +x-face of the original element.
            If the direction of the shear stress is opposite to the positive
            y-direction, it must have a negative value.

        Notes
        -----
        The coordinate system in which the element is observed is illustrated in
        the document *right_handed_coordinate_system.pdf* in folder *docs/notes*.
        """
        self.sigma_x = sigma_x.to('N / mm**2')
        self.sigma_y = sigma_y.to('N / mm**2')
        self.tau_xy = tau_xy.to('N / mm**2')

    def transform(self, theta: Quantity) -> tuple[Quantity, ...]:
        """Returns the normal stresses and the shear stress acting on a surface
        element which is rotated an angle `theta` about the z-axis with respect
        to the original element. A positive value for `theta` corresponds with a
        counterclockwise rotation in the applied coordinate system.
        """
        theta = theta.to('rad')
        sigma_x = self._transform_sigma(theta)
        sigma_y = self._transform_sigma(theta + np.pi / 2)
        tau_xy = self._transform_tau(theta)
        return sigma_x, sigma_y, tau_xy

    def _transform_sigma(self, theta: Quantity) -> Quantity:
        a = (self.sigma_x + self.sigma_y) / 2
        b = (self.sigma_x - self.sigma_y) / 2 * np.cos(2 * theta)
        c = self.tau_xy * np.sin(2 * theta)
        sigma = a + b + c
        return sigma

    def _transform_tau(self, theta: Quantity) -> Quantity:
        a = -(self.sigma_x - self.sigma_y) / 2 * np.sin(2 * theta)
        b = self.tau_xy * np.cos(2 * theta)
        tau = a + b
        return tau

    def _principal_normal_stress_angle(self) -> Quantity:
        n = self.tau_xy
        d = (self.sigma_x - self.sigma_y) / 2
        theta = np.arctan2(n, d) / 2
        return theta

    def principal_normal_stresses(self) -> tuple[Quantity, Quantity, Quantity]:
        """Returns the maximum and minimum normal stress acting on a surface
        element at the given point in the cross-section, together with the
        rotation angle about the z-axis measured from the x-axis to the
        direction of the maximum normal stress.
        """
        theta_a = self._principal_normal_stress_angle()
        theta_b = theta_a + Q_(np.pi / 2, 'rad')
        sigma_a = self._transform_sigma(theta_a)
        sigma_b = self._transform_sigma(theta_b)
        if sigma_a > sigma_b:
            sigma_1 = sigma_a
            sigma_2 = sigma_b
            theta_1 = theta_a
        else:
            sigma_1 = sigma_b
            sigma_2 = sigma_a
            theta_1 = theta_b
        return sigma_1, sigma_2, theta_1

    def _maximum_shear_stress_angle(self) -> float:
        n = -(self.sigma_x - self.sigma_y) / 2
        d = self.tau_xy
        theta = np.arctan2(n, d) / 2
        return theta

    def maximum_shear_stress(self) -> tuple[Quantity, Quantity, Quantity]:
        """Returns the maximum in-plane shear stress and the average normal
        stress acting on a surface element at the given point in the
        cross-section, together with the rotation angle about the z-axis
        measured from the x-axis to the direction of the average normal stress.
        """
        theta = self._maximum_shear_stress_angle()
        tau_max = self._transform_tau(theta)
        sigma = self._transform_sigma(theta)
        return tau_max, sigma, theta

    def absolute_maximum_shear_stress(self) -> tuple[Quantity, Quantity]:
        """Returns the absolute maximum shear stress and the average normal
        stress acting on a surface element at the given point in the
        cross-section, together with the rotation angle measured from the
        x-axis to the direction of the average normal stress.

        Notes
        -----
        If the in-plane principal stresses both have the same sign, the absolute
        maximum shear stress will occur out of the xy-plane. This shear stress
        is also greater than the in-plane maximum shear stress.

        If the in-plane principal stresses are of opposite signs, then the
        absolute maximum shear stress will equal the maximum in-plane shear
        stress.
        """
        sigma_1, sigma_2, _ = self.principal_normal_stresses()
        if np.sign(sigma_1) == np.sign(sigma_2):
            tau_abs_max = sigma_1 / 2
            sigma = sigma_1 / 2
        else:
            tau_abs_max, sigma, _ = self.maximum_shear_stress()
        return tau_abs_max, sigma

    def mohr_circle(self):
        """Returns the center on the σ-axis and the radius of Mohr's circle."""
        sigma_avg = (self.sigma_x + self.sigma_y) / 2
        R = np.sqrt(((self.sigma_x - self.sigma_y) / 2) ** 2 + self.tau_xy ** 2)
        return sigma_avg, R

    def draw_mohr_circle(self) -> LineChart:
        """Returns a `LineChart` object with a drawing of Mohr's circle.

        Notes
        -----
        Call method `show()` on the `LineChart` object to actually see the
        drawing.
        """
        sigma_avg, R = self.mohr_circle()
        angles = np.linspace(0, 2 * np.pi, 100)
        x_values = R * np.cos(angles) + sigma_avg
        y_values = R * np.sin(angles)
        x_min = min(x_values.m)
        x_max = max(x_values.m)
        y_min = min(y_values.m)
        y_max = max(y_values.m)
        b = .05 * max(x_max - x_min, y_max - y_min)
        chart = LineChart()
        chart.axes.set_xlim(xmin=x_min - b, xmax=x_max + b)
        chart.axes.set_ylim(ymin=y_min - b, ymax=y_max + b)
        chart.axes.invert_yaxis()
        chart.axes.set_aspect(1)
        chart.axes.axhline(0, color='black', lw=1.5)
        chart.axes.axvline(0, color='black', lw=1.5)
        chart.add_xy_data(
            label='circle',
            x1_values=x_values.m,
            y1_values=y_values.m,
            style_props={'color': 'tab:blue'}
        )
        chart.add_xy_data(
            label='center',
            x1_values=[sigma_avg.m],
            y1_values=[0],
            style_props={
                'marker': 'o',
                'linestyle': 'none',
                'color': 'tab:blue'
            }
        )
        chart.add_xy_data(
            label='reference point',
            x1_values=[self.sigma_x.m],
            y1_values=[self.tau_xy.m],
            style_props={
                'marker': 'o',
                'linestyle': 'none',
                'color': 'tab:blue'
            }
        )
        chart.add_xy_data(
            label='reference line',
            x1_values=[sigma_avg.m, self.sigma_x.m],
            y1_values=[0, self.tau_xy.m],
            style_props={
                'linestyle': '--',
                'color': 'tab:blue'
            }
        )
        chart.x1.add_title('sigma, MPa')
        chart.y1.add_title('tau, MPa')
        return chart

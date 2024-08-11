import numpy as np
from mechanics import Quantity
from mechanics.geometry.shapes import Shape, Polygon, Circle, Annulus, HollowPolygon
from mechanics.geometry import AxesRotation2D, Line


Q_ = Quantity


class Bending:
    """
    Determine the (maximum) normal stress in a cross-section associated with
    the resultant internal bending moment.
    """
    def __init__(
        self,
        M: Quantity | tuple[Quantity, Quantity],
        shape: Shape
    ) -> None:
        """Creates a `Bending` instance.

        Parameters
        ----------
        M:
            Components of the resultant internal bending moment at the
            cross-section.
            If only one component is given, it is assumed to be the bending
            moment about about the horizontal axis of the cross-section.
            If two components are given, the second component must be the
            bending moment about the vertical axis of the cross-section.
        shape:
            Shape of the cross-section.
        """
        if isinstance(M, tuple):
            self.M_z, self.M_y = M[0], M[1]
        else:
            self.M_z = M
            self.M_y = Q_(0.0, M.units)
        self.shape = shape
        self.__set_principal_components()
        self.sigma_max = self.__get_sigma_max()

    def sigma(
        self,
        z: Quantity,
        y: Quantity
    ) -> Quantity:
        """
        Returns the normal stress at the cross-section's point with coordinates
        `z` and `y` where `z` is the horizontal coordinate and `y` is the
        vertical coordinate.

        Notes
        -----
        A negative sigma-value indicates a tensile stress, while a positive
        value indicates a compression stress (see notes about the right-handed
        coordinate system).
        """
        z_acc, y_acc = self._rot((z, y))
        sigma = -self._M_maj * y_acc / self._I_maj
        sigma += self._M_min * z_acc / self._I_min
        return sigma

    def section_modulus(self) -> Quantity:
        """Returns the section modulus of the beam."""
        I_xx = self.shape.moment_of_inertia_xx
        c = self.__farthest_vertex(self.shape)[1]
        S = I_xx / c
        return S

    def __set_principal_components(self) -> None:
        """
        Determines the components of the bending moment about the principal
        axes.
        """
        # Get moments of inertia along major and minor principal axis and the
        # angle between the horizontal axis and the major principal axis:
        self._I_maj, self._I_min, self._theta = self.shape.principal_moments_of_inertia
        # Define the rotation of the coordinate system to align with the
        # principal axes of the cross-section:
        self._rot = AxesRotation2D(self._theta)
        # Get components of the bending moment about the principal axes:
        self._M_maj, self._M_min = self._rot((self.M_z, self.M_y))

    def __neutral_axis_angle(self) -> Quantity:
        """
        Returns the orientation angle of the neutral axis, being the angle
        measured from the major principal axis to the neutral axis.
        """
        gamma = np.arctan2(  # angle between bending moment and major principal axis
            self._M_min.m,
            self._M_maj.m
        )
        alpha = np.arctan2(  # angle between neutral axis and major principal axis
            self._I_maj * np.sin(gamma),
            self._I_min * np.cos(gamma)
        )
        return Q_(alpha, 'rad')

    def __farthest_vertex(self, shape: Shape) -> Quantity:
        """
        Returns the vertex of the polygon at the farthest distance from the
        neutral axis.
        """
        if isinstance(shape, Polygon):
            # Create neutral axis as a `Line` object in the coordinate system
            # defined by the principal axes. The neutral axis always goes
            # through the centroid of the cross-section which is the origin of
            # the coordinate system.
            neutral_axis = Line.from_point_and_angle(
                p=(0.0, 0.0),
                angle=self.__neutral_axis_angle().m
            )
            # Transform the coordinates of the vertices to the principal axes:
            vertices_rotated = [self._rot(p) for p in shape.vertices]
            # Get the distance of each vertex to the neutral axis:
            distances = [
                neutral_axis.distance((p[0].m, p[1].m))
                for p in vertices_rotated
            ]
            # Get the index of the largest distance:
            i_farthest = np.argmax(np.array(distances))
            # This is also the index of the vertex farthest from the neutral
            # axis:
            p_farthest = shape.vertices[i_farthest]
            return p_farthest

    def __get_sigma_max(self) -> Quantity:
        """
        Returns the maximum bending stress in the cross-section.
        """
        if isinstance(self.shape, Polygon):
            p_farthest = self.__farthest_vertex(self.shape)
            sigma_max = self.sigma(p_farthest[0], p_farthest[1])
            return sigma_max
        elif isinstance(self.shape, HollowPolygon):
            p_farthest = self.__farthest_vertex(self.shape.outer_polygon)
            sigma_max = self.sigma(p_farthest[0], p_farthest[1])
            return sigma_max
        elif isinstance(self.shape, (Circle, Annulus)):
            if isinstance(self.shape, Circle):
                r = self.shape.dim.radius
            else:
                r = self.shape.dim.outer_radius
            sigma_max = self.sigma(Q_(0.0, r.units), r)
            return sigma_max

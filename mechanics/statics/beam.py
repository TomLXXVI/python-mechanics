from collections.abc import Sequence
import sympy as sp
import numpy as np
from mechanics import Quantity
from .force import Force, Moment
from .supports import Support
from .distributed_load import DistributedLoad1D
from ._units import Units

Q_ = Quantity


class Beam2D:
    """Represents a planar beam.

    The origin of the right-handed coordinate system is assumed to be always on
    the left end of the beam. The positive x-axis direction points to the right.
    along the centerline of the beam. The positive y-axis direction points
    upward.
    """
    def __init__(
        self,
        length: Quantity,
        loads: Sequence[Force | Moment | DistributedLoad1D] | None = None,
        supports: Sequence[Support] | None = None
    ):
        """Creates a `Beam2D` instance.

        Parameters
        ----------
        length:
            Length of the beam.
        loads:
            Sequence of known external loads exerted on the beam (concentrated
            forces, moments, and/or distributed loadings).
        supports:
            Sequence of supports that support the beam and of which the reaction
            forces exerted on the beam are yet unknown.
        """
        self.length = length
        self.action_forces: list[Force] = []
        self.action_torques: list[Moment] = []
        self.distributed_loads: list[DistributedLoad1D] = []
        self.supports: list[Support] = []
        self.__set_units()
        for load in loads:
            if isinstance(load, Force):
                self.add_forces(load)
            elif isinstance(load, DistributedLoad1D):
                self.add_distributed_loads(load)
            elif isinstance(load, Moment):
                self.add_torques(load)
            else:
                raise ValueError("not a valid load")
        if supports is not None:
            for support in supports:
                if isinstance(support, Support):
                    self.add_supports(support)
                else:
                    raise ValueError("not a valid support")

    def add_forces(self, *forces: Force) -> None:
        """Applies one or more known action forces to the beam."""
        self.action_forces.extend(forces)

    def add_supports(self, *supports: Support) -> None:
        """Adds one or more supports to the beam."""
        self.supports.extend(supports)

    def add_distributed_loads(self, *loads: DistributedLoad1D) -> None:
        """Places one or more distributed loads on the beam."""
        self.distributed_loads.extend(loads)

    def add_torques(self, *torques: Moment) -> None:
        """Applies one or more torques to the beam."""
        self.action_torques.extend(torques)

    def get_reactions(self) -> dict[str, Force | Moment]:
        """Returns the solutions for the unknown reaction forces applied by
        the beam supports.

        Returns
        -------
        A tuple holding two dictionaries. The keys are the names that were given
        to the supports. The first dictionary holds the forces as `Force`
        objects. The second dictionary holds the torques as `Moment` objects.
        If no reaction torques are present, the second dictionary will be empty.
        """
        # Get the action forces (point forces and from distributed loads)
        # applied to the whole beam:
        x_min = Q_(0, self.__u_pos)
        x_max = self.length
        # Get the known action forces acting on the beam:
        forces_a = self.__get_action_forces(x_min, x_max)
        forces_d = self.__get_distributed_forces(x_min, x_max)
        forces = forces_a + forces_d
        # Get the components of the known action forces:
        tup_F_x, tup_F_y, _ = zip(*[F.components for F in forces])
        # Get the moments about the z-axis of the known action forces and torques:
        *_, tup_M_z_1 = zip(*[F.moment().components for F in forces])
        torques_a = self.__get_action_torques(x_min, x_max)
        if torques_a:
            *_, tup_M_z_2 = zip(*[T.components for T in torques_a])
            lst_M_z = list(tup_M_z_1) + list(tup_M_z_2)
        else:
            lst_M_z = list(tup_M_z_1)
        # Take the sum of the components:
        F_a_x = sum(F_x.to(self.__u_force).m for F_x in tup_F_x)
        F_a_y = sum(F_y.to(self.__u_force).m for F_y in tup_F_y)
        M_a_z = sum(M_z.to(self.__u_moment).m for M_z in lst_M_z)
        # Get the components of the unknown reaction forces and torques exerted
        # by the supports on the beam:
        F_r_x = sum([sup.F_x for sup in self.supports if sup.F_x is not None])
        F_r_y = sum([sup.F_y for sup in self.supports if sup.F_y is not None])
        M_r_z = sum([sup.M_z for sup in self.supports if sup.M_z is not None])
        # Get the moments of the reaction forces about the origin:
        M_r_z += sum([
            sup.position[0].to(self.__u_pos).m * sup.F_y
            for sup in self.supports
            if sup.F_y is not None
        ])
        M_r_z += sum([
            -sup.position[1].to(self.__u_pos).m * sup.F_x
            for sup in self.supports
            if sup.F_x is not None
        ])
        # Create the static equilibrium equations and solve for the components
        # of the unknown reaction forces:
        equations = []
        if F_r_x:
            F_x_bal = sp.Eq(F_r_x + F_a_x, 0)
            equations.append(F_x_bal)
        if F_r_y:
            F_y_bal = sp.Eq(F_r_y + F_a_y, 0)
            equations.append(F_y_bal)
        if M_r_z:
            M_z_bal = sp.Eq(M_r_z + M_a_z, 0)
            equations.append(M_z_bal)
        sol = sp.solve(equations)
        # Add the solved reaction forces and/or reaction moments to the
        # dictionary of known action forces, respectively the dictionary of
        # known action moments:
        sol = self.__add_reactions(sol)
        return sol

    def cut(self, x: Quantity, side: str = 'left') -> dict[str, Quantity]:
        """Returns the internal normal force, the internal shear force, and the
        internal bending moment as `Quantity` objects, which are present in the
        cross-section at position `x` along the beam. A cut can be located
        on the left or on the right part of the beam, which is indicated by
        parameter `side`.
        """
        if side == 'left':
            x_min = Q_(0, self.__u_pos)
            x_max = x
        else:
            x_min = x
            x_max = self.length
        # Get the known action forces (point forces and/or from distributed
        # loads) between `x_min` and `x_max`:
        forces_a = self.__get_action_forces(x_min, x_max)
        forces_d = self.__get_distributed_forces(x_min, x_max)
        forces = forces_a + forces_d
        # Get the x- and y- components of the known action forces and their
        # moments (about the origin of the coordinate system) and sum them:
        if forces:
            tup_F_x, tup_F_y, _ = zip(*[F.components for F in forces])
            *_, tup_M_z = zip(*[F.moment().components for F in forces])
            F_a_x = sum(F_x.to(self.__u_force).m for F_x in tup_F_x)
            F_a_y = sum(F_y.to(self.__u_force).m for F_y in tup_F_y)
            M_a_z = sum(M_z.to(self.__u_moment).m for M_z in tup_M_z)
        else:
            F_a_x = 0.0
            F_a_y = 0.0
            M_a_z = 0.0
        # Get the known action torques between `x_min` and `x_max`:
        torques = self.__get_action_torques(x_min, x_max)
        # Get the z-components of the known action torques and add them to
        # the sum of the moments of forces:
        if torques:
            *_, tup_T_z = zip(*[T.components for T in torques])
            M_a_z += sum(T_z.to(self.__u_moment).m for T_z in tup_T_z)
        # Create the static equilibrium equations and solve for the components
        # of the unknown reaction forces:
        F_i_x, F_i_y, M_i_z = sp.symbols(['F_x', 'F_y', 'M_z'])
        F_x_bal = sp.Eq(F_i_x, -F_a_x)
        F_y_bal = sp.Eq(F_i_y, -F_a_y)
        M_z_bal = sp.Eq(M_i_z + x.to(self.__u_pos).m * F_i_y, -M_a_z)
        sol_dict = sp.solve([F_x_bal, F_y_bal, M_z_bal])
        sol = self.__create_internal_forces(sol_dict)
        return sol

    def __set_units(self) -> None:
        self.__u_pos = Units.u_pos
        self.__u_force = Units.u_force
        self.__u_moment = Units.u_moment

    def __get_action_forces(
        self,
        x_min: Quantity,
        x_max: Quantity
    ) -> list[Force]:
        """Returns the action forces which applied to the beam between positions
        `x_min` and `x_max`.
        """
        forces = []
        for force in self.action_forces:
            x = force.action_point[0]
            if x_min <= x <= x_max:
                forces.append(force)
        return forces

    def __get_distributed_forces(
        self,
        x_min: Quantity,
        x_max: Quantity
    ) -> list[Force]:
        """Returns the resultants of distribution loads which are present
        between positions `x_min` and `x_max` on the beam.
        """
        forces = []
        for load in self.distributed_loads:
            x1, _ = load.points[0]
            x2, _ = load.points[-1]
            if x1 >= x_min and x2 <= x_max:
                F = load.resultant()
                forces.append(F)
            elif x1 < x_max < x2:
                F = load.resultant(x1, x_max)
                forces.append(F)
            elif x1 < x_min < x2:
                F = load.resultant(x_min, x2)
                forces.append(F)
        return forces

    def __get_action_torques(
        self,
        x_min: Quantity,
        x_max: Quantity
    ) -> list[Moment]:
        """Returns the action torques applied to the beam between positions
        `x_min` and `x_max`.
        """
        torques = []
        for torque in self.action_torques:
            x = torque.action_point[0]
            if x_min <= x <= x_max:
                torques.append(torque)
        return torques

    def __add_reactions(
        self,
        sol_dict: dict[str, sp.Float]
    ) -> dict[str, Force | Moment]:
        """Used in method `get_reactions`. Takes the Sympy solutions of the
        reaction forces and torques, and adds these to the list of known action
        forces and/or to the list of known action torques.
        Returns two dictionaries: the first contains the solved reaction forces
        as `Force` objects; the second contains the solved reaction torques as
        `Moment` objects.
        """
        reactions = {}
        for key, value in sol_dict.items():
            s = str(key).split('.')
            name, component = s[0], s[1]
            value = float(value)
            if name in reactions.keys():
                reactions[name][component] = value
            else:
                reactions[name] = {component: value}
        for name, d in reactions.items():
            F_mag: float | None = d.get('F', None)
            F_x: float | None = d.get('F_x', None)
            F_y: float | None = d.get('F_y', None)
            if F_mag is not None:
                F = self.__create_reaction_force(name, None, None, F=F_mag)
                self.action_forces.append(F)
                reactions[name] = F
            elif F_x is not None and F_y is None:
                F = self.__create_reaction_force(name, F_x, 0.0)
                self.action_forces.append(F)
                reactions[name] = F
            elif F_x is None and F_y is not None:
                F = self.__create_reaction_force(name, 0.0, F_y)
                self.action_forces.append(F)
                reactions[name] = F
            elif F_x is not None and F_y is not None:
                F = self.__create_reaction_force(name, F_x, F_y)
                self.action_forces.append(F)
                reactions[name] = F
            else:
                pass
            M_z: float | None = d.get('M_z', None)
            if M_z is not None:
                M = self.__create_reaction_torque(M_z)
                self.action_torques.append(M)
                reactions[name] = M
        return reactions

    def __create_reaction_force(
        self,
        name: str,
        F_x: float | None,
        F_y: float | None,
        F: float = None
    ) -> Force:
        """Creates a `Force` object from the Sympy solution of a reaction
        force, when using the `get_reactions` method.
        """
        d = {support.name: support for support in self.supports}
        support = d[name]
        if F is None:
            F = Q_(np.sqrt(F_x ** 2 + F_y ** 2), self.__u_force)
            theta = Q_(np.arctan2(F_y, F_x), 'rad')
            force = Force(support.position, F, theta)
        else:
            if F < 0:
                F = Q_(abs(F), self.__u_force)
                theta = support.theta + Q_(np.pi, 'rad')
            else:
                F = Q_(F, self.__u_force)
                theta = support.theta
            force = Force(support.position, F, theta)
        return force

    def __create_reaction_torque(self, M_z: float) -> Moment:
        """Creates a `Moment` object from the Sympy solution of a reaction
        torque, when using the `get_reactions` method.
        """
        M = Moment.create_from_components(M_z=Q_(M_z, self.__u_moment))
        return M

    def __create_internal_forces(self, sol_dict: dict[str, sp.Float]) -> dict[str, Quantity]:
        """Creates a dictionary with the Sympy solutions of the internal forces
        as `Quantity` objects, when using the `cut` method.
        """
        d = {}
        for key, value in sol_dict.items():
            name = str(key)
            value = round(float(value), 12)
            if name == 'M_z':
                d[name] = Q_(value, self.__u_moment)
            else:
                d[name] = Q_(value, self.__u_force)
        return d

    def diagrams_data(self, num: int = 50) -> tuple[Quantity, ...]:
        """Returns `Quantity`-arrays `x_arr`, `N_arr`, `V_arr`, and `M_arr` with
        `x_arr` an array of equally spaced positions along the longitudinal axis
        of the beam, `N_arr` an array with the internal normal force at each
        position, `V_arr` an array with the internal shear force at each
        position, and `M_arr` an array with the internal bending moment at each
        position.
        """
        x_min = 0.0
        x_max = self.length.m
        x_arr = Q_(np.linspace(x_min, x_max, num, endpoint=True), self.length.units)
        irl_list = [self.cut(x) for x in x_arr]
        N_arr = Quantity.from_list([irl['F_x'] for irl in irl_list])
        V_arr = Quantity.from_list([irl['F_y'] for irl in irl_list])
        M_arr = Quantity.from_list([irl['M_z'] for irl in irl_list])
        return x_arr, N_arr, V_arr, M_arr

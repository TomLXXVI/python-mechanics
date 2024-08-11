from concurrent.futures import ProcessPoolExecutor
import sympy as sp
import numpy as np
from scipy.interpolate import interp1d
from mechanics import Quantity
from mechanics.charts import LineChart
from .vector import Force, Moment, DistributedLoad1D, Angle, Position


Q_ = Quantity


class System:
    """Solves a system of forces and/or moments acting on a arbitrary body for
    the unknown forces and/or moments such that static equilibrium of the body
    is preserved.

    Notes
    -----
    In case of a 3D system, there are 6 equations available, so the system cannot
    have more than 6 unknowns (magnitude and/or direction).
    In case of a 2D system, there are only 3 equations available, so the system
    cannot have more 3 unknowns.
    """
    _units_of_force = 'N'
    _units_of_moment = 'N * m'
    _units_of_length = 'm'
    _units_of_load = 'N / m'

    def __init__(
        self,
        loadings: list[Force | Moment | DistributedLoad1D],
        units: tuple[str, str] | None = None
    ) -> None:
        """Creates a `System` object of which the loadings can be forces,
        moments (moments), or distributed linear loads.

        Parameters
        ----------
        loadings:
            List with the forces, moments (moments) and distributed loads acting
            on the body.
        units:
            Tuple of two strings. The first string contains the units of force
            to be used (default units of force are 'N', Newtons). The second
            string contains the units of length (default units of length are
            'm', meters). If `units` is left to `None`, the default units are
            used.
            From the given units of force and of length, the units of
            distributed loads and of moments (moments) are derived.
            Note that the units are set at class level, i.e. when the units are
            set on one `System` object, they will be inherited as default units
            when instantiating other subsequent `System` objects.

        Notes
        -----
        Distributed linear loads cannot have unknowns; these should always be
        fully determined. Distributed loads will be replaced by their equivalent
        resultant force.
        """
        if isinstance(units, tuple):
            self.__class__.units(units[0], units[1])
        self.external_forces: dict[str, Force] = self.__get_external_forces(*loadings)
        self.external_moments: dict[str, Moment] = self.__get_external_moments(*loadings)
        self.external_distributed_loads: dict[str, DistributedLoad1D] = self.__get_external_distributed_loads(*loadings)

    def __get_external_forces(self, *loadings) -> dict[str, Force]:
        forces = [
            f.to(self._units_of_force)  # all forces must have the same units
            for f in loadings
            if isinstance(f, Force)
        ]
        for force in forces: force.position.to(self._units_of_length)  # all positions must have the same units
        forces = {f.name: f for f in forces}
        return forces

    def __get_external_moments(self, *loadings) -> dict[str, Moment]:
        moments = [
            m.to(self._units_of_moment)  # all moments must have the same units
            for m in loadings
            if isinstance(m, Moment)
        ]
        for moment in moments: moment.position.to(self._units_of_length)  # all positions must have the same units
        moments = {t.name: t for t in moments}
        return moments

    def __get_external_distributed_loads(self, *loadings) -> dict[str, DistributedLoad1D]:
        distributed_loads = [
            dl.to(self._units_of_load)  # all loads must have the same units
            for dl in loadings
            if isinstance(dl, DistributedLoad1D)
        ]
        for dl in distributed_loads: dl.positions(self._units_of_length)  # all positions must have the same units
        distributed_loads = {dl.name: dl for dl in distributed_loads}
        return distributed_loads

    def __contains_unknowns(self) -> bool:
        # Checks if there are any unknown forces/moments acting on the body.
        loadings = list(self.external_forces.values())
        loadings.extend(self.external_moments.values())
        for load in loadings:
            if load.is_symbolic():
                return True
        return False

    @classmethod
    def units(cls, units_of_force: str, units_of_length: str) -> None:
        """Sets the units of force and of length to be used when instantiating
        a `System` object.
        """
        cls._units_of_force = units_of_force
        cls._units_of_length = units_of_length
        cls._units_of_moment = f"{units_of_force} * {units_of_length}"
        cls._units_of_load = f"{units_of_force} / {units_of_length}"

    def solve(self) -> dict[str, Force | Moment] | None:
        """Solves the system for any unknown forces and/or moments (moments).

        Returns
        -------
        A dictionary of which the keys are the names of the forces or moments
        and the values are the corresponding `Force` or `Moment` objects.
        Returns `None` if there are no unknown forces and/or moments (moments)
        acting on the body.
        """
        if self.__contains_unknowns():
            components_dict = self.__decompose_into_components()
            equations = self.__create_equations(components_dict)
            sol_dict = self.__solve_with_sympy(equations)
            sol_dict = self.__create_vectors(sol_dict)
            self.__replace_unknowns(sol_dict)
            return sol_dict

    def __decompose_into_components(self):
        # Decomposes the forces and moments in the system into their cartesian
        # components which are held in a dictionary where unknown and known
        # components are separated.
        keys = ['F_x', 'F_y', 'F_z', 'M_x', 'M_y', 'M_z']
        values = [{'unknown': [], 'known': []} for _ in range(len(keys))]
        components_dict = dict(zip(keys, values))

        def __add_components_to_dict(keys, components):
            for key, component in zip(keys, components):
                if isinstance(component, sp.Expr):
                    components_dict[key]['unknown'].append(component)
                else:
                    components_dict[key]['known'].append(component)

        def __decompose_force(force: Force):
            keys = ['F_x', 'F_y', 'F_z', 'M_x', 'M_y', 'M_z']
            components = list(force.component_values)
            components.extend(list(force.moment().component_values))
            __add_components_to_dict(keys, components)

        def __decompose_moment(moment: Moment):
            keys = ['M_x', 'M_y', 'M_z']
            components = list(moment.component_values)
            __add_components_to_dict(keys, components)

        def __decompose_distributed_load(distr_load: DistributedLoad1D):
            keys = ['F_x', 'F_y', 'F_z', 'M_x', 'M_y', 'M_z']
            resultant = distr_load.resultant()
            components = list(resultant.component_values)
            components.extend(list(resultant.moment().component_values))
            __add_components_to_dict(keys, components)

        for force in self.external_forces.values(): __decompose_force(force)
        for moment in self.external_moments.values(): __decompose_moment(moment)
        for distr_load in self.external_distributed_loads.values(): __decompose_distributed_load(distr_load)
        return components_dict

    @staticmethod
    def __create_equations(components_dict: dict) -> list[sp.Eq]:
        # Creates the system of equations that will be solved with Sympy.
        def __create_equation(d: dict) -> sp.Eq | None:
            unknown = sum(expr for expr in d['unknown'])
            if unknown == 0:
                return None
            else:
                known = sum(value for value in d['known'])
                eq = sp.Eq(unknown + known, 0)
                return eq

        eqs = []
        for d in components_dict.values():
            eq = __create_equation(d)
            if eq is not None:
                eqs.append(eq)
        return eqs

    @staticmethod
    def __solve_with_sympy(equations: list[sp.Eq]) -> dict:
        # Solves the system of equations for the unknown components.
        sol_dict = sp.solve(equations, dict=True)
        if isinstance(sol_dict, list):
            sol_dict = sol_dict[0]
        return sol_dict

    def __create_vectors(self, sol_dict: dict) -> dict[str, Force | Moment]:
        # Create `Force` and/or `Moment` objects with the solved components
        # returned by Sympy in `sol_dict`.
        solutions = {}
        for symbol, solution in sol_dict.items():
            s = str(symbol).split('.')
            solution = float(solution)
            try:
                name, component = s[0], s[1]
            except IndexError:
                # If only the magnitude was an unknown, the symbol won't contain
                # a dot (see class `Vector`, method `__symbolic_vector`).
                name = s[0]
                solutions[name] = {'magnitude': solution}
            else:
                if name in solutions.keys():
                    # There is already at least one component added...
                    solutions[name][component] = solution
                else:
                    solutions[name] = {component: solution}

        unknowns = {f.name: f for f in self.external_forces.values() if f.is_symbolic()}
        unknowns.update({t.name: t for t in self.external_moments.values() if t.is_symbolic()})

        def __create_vector(
            name: str,
            components: dict[str, float],
            unknown: Force | Moment
        ) -> Force | Moment:
            _Vector = Force if isinstance(unknown, Force) else Moment
            mag = components.get('magnitude', None)
            if mag is not None:
                # Only the magnitude was an unknown...
                if mag < 0:
                    mag = abs(mag)
                    theta = unknown.theta + Angle(180)
                    if unknown.gamma.magnitude != 0.0:
                        gamma = unknown.gamma + Angle(180)
                    else:
                        gamma = unknown.gamma
                else:
                    theta = unknown.theta
                    gamma = unknown.gamma
                vec = _Vector(
                    magnitude=mag,
                    theta=theta,
                    gamma=gamma,
                    position=unknown.position,
                    units=unknown.units,
                    name=name
                )
                return vec
            else:
                vec_x = components.get('x', 0.0)
                vec_y = components.get('y', 0.0)
                vec_z = components.get('z', 0.0)
                vec = _Vector.create_from_components(
                    vec_x, vec_y, vec_z,
                    position=unknown.position,
                    units=unknown.units,
                    name=name
                )
                # noinspection PyTypeChecker
                return vec

        for name, components in solutions.items():
            unknown = unknowns[name]
            vector = __create_vector(name, components, unknown)
            solutions[name] = vector  # either a `Force` or a `Moment` object

        return solutions

    def __replace_unknowns(self, solutions: dict[str, Force | Moment]) -> None:
        # Replaces the unknown forces or moments with their solution.
        for name, sol in solutions.items():
            if isinstance(sol, Force):
                self.external_forces[name] = sol
            if isinstance(sol, Moment):
                self.external_moments[name] = sol
        return None


class Beam(System):
    """Defines a beam: a  structural element that is in general a slender,
    straight bar having a constant cross-sectional area.
    """
    num_decimals: int = 9

    def __init__(
        self,
        length: Quantity,
        loadings: list[Force | Moment | DistributedLoad1D],
        units: tuple[str, str] | None = None,
        num_sections: int = 50
    ) -> None:
        """Creates a `Beam` object.

        Parameters
        ----------
        length:
            Length of the beam.
        loadings:
            List of external loadings acting on the beam (`Force`-objects,
            `Moment`-objects, and/or `DistributedLoad1D`-objects).
            Unknown reaction forces or moments are solved on instantiation.
        units:
            Tuple of two strings. The first string contains the units of force
            to be used (default units of force are 'N', Newtons). The second
            string contains the units of length (default units of length are
            'm', meters). If `units` is left to `None`, the default units are
            used.
            From the given units of force and of length, the units of
            distributed loads and of moments (moments) are derived.
            Note that the units are set at class level, i.e. when the units are
            set on one `Beam` object, they will be inherited as default units
            when instantiating other subsequent `Beam` objects.
        num_sections:
            To determine the profiles of the normal force, shear force and
            bending moment along the beam, the beam is cut at a number of
            equally spaced cross-sections. The number of cuts can be specified
            via parameter `num_sections` (default number of cuts is 50).

        Notes
        -----
        The longitudinal axis of the beam is taken as the x-axis of the
        coordinate system.
        """
        super().__init__(loadings, units)
        self.length = length
        self.num_sections = num_sections
        self._length: float = length.to(self._units_of_length).m
        # Solve for any unknown external reaction forces and/or moments:
        self.solve()
        # Create profiles of the resultant internal loadings along the
        # longitudinal axis of the beam (normal force, shear force, bending
        # moment, and torsional moment):
        t = self.__profiles_of_internal_loadings()
        self._x_arr = t[0]
        self._N_arr = t[1]
        self._V_arr = t[2]
        self._M_arr = t[3]
        self._T_arr = t[4]
        self._V_interp = interp1d(self._x_arr, self._V_arr)
        self._M_interp = interp1d(self._x_arr, self._M_arr)
        self._N_interp = interp1d(self._x_arr, self._N_arr)
        self._T_interp = interp1d(self._x_arr, self._T_arr)

    def cut(
        self,
        x: Quantity,
        view: str = 'left'
    ) -> tuple[Force, Moment]:
        """Makes a cross-sectional cut through the beam at the given x-position
        and determines the resultant internal force and/or moment acting at this
        section.
        Parameter `view` indicates whether the cross-section of the left part of
        the cut beam is to be regarded (this is the default) or the
        cross-section of its right part (`view = 'right'`).

        Returns
        -------
        Tuple with the internal force (`Force`-object) and internal moment
        (`Moment`-object) acting at the viewed cross-section.
        """
        x = x.to(self._units_of_length).m
        y = 0.0
        z = 0.0
        if view == 'left':
            x_min = 0.0
            x_max = x
        else:
            x_min = x
            x_max = self._length

        # List with the external forces between `x_min` and `x_max`.
        forces = []
        for force in self.external_forces.values():
            # noinspection PyProtectedMember
            if x_min <= force.position._x <= x_max:
                forces.append(force)

        # Add resultants of distributed linear loads between x_min and x_max
        # to the list of forces.
        for distr_load in self.external_distributed_loads.values():
            x1 = distr_load.x_coords[0].m
            x2 = distr_load.x_coords[-1].m
            if x1 >= x_min and x2 <= x_max:
                force = distr_load.resultant().to(self._units_of_force)
                forces.append(force)
            elif x1 < x_max < x2:
                force = distr_load.resultant(x1, x_max).to(self._units_of_force)
                forces.append(force)
            elif x1 < x_min < x2:
                force = distr_load.resultant(x_min, x2).to(self._units_of_force)
                forces.append(force)

        # List with the external moments between x_min and x_max.
        moments = []
        for moment in self.external_moments.values():
            # noinspection PyProtectedMember
            if x_min <= moment.position._x <= x_max:
                moments.append(moment)

        # Get components of external forces.
        if forces:
            tupF_x, tupF_y, tupF_z = zip(*[f.component_values for f in forces])
            tupM_x, tupM_y, tupM_z = zip(*[f.moment().component_values for f in forces])
            F_x, F_y, F_z = sum(tupF_x), sum(tupF_y), sum(tupF_z)
            M_x, M_y, M_z = sum(tupM_x), sum(tupM_y), sum(tupM_z)
        else:
            F_x = F_y = F_z = 0.0
            M_x = M_y = M_z = 0.0

        # Get components of external moments.
        if moments:
            tupT_x, tupT_y, tupT_z = zip(*[t.component_values for t in moments])
            T_x, T_y, T_z = sum(tupT_x), sum(tupT_y), sum(tupT_z)
        else:
            T_x = T_y = T_z = 0.0
        M_x += T_x
        M_y += T_y
        M_z += T_z

        # Solve for the components of the internal force and moment at the given
        # section:
        int_F_x, int_F_y, int_F_z = sp.symbols(['F_x', 'F_y', 'F_z'])
        int_M_x, int_M_y, int_M_z = sp.symbols(['M_x', 'M_y', 'M_z'])

        int_M_x += -z * int_F_y + y * int_F_z
        int_M_y += z * int_F_x - x * int_F_z
        int_M_z += -y * int_F_x + x * int_F_y

        eqs = [
            sp.Eq(int_F_x, -F_x),
            sp.Eq(int_F_y, -F_y),
            sp.Eq(int_F_z, -F_z),
            sp.Eq(int_M_x, -M_x),
            sp.Eq(int_M_y, -M_y),
            sp.Eq(int_M_z, -M_z)
        ]
        sol_dict = sp.solve(eqs)
        sol_dict = {
            str(key): round(float(val), self.num_decimals)
            for key, val in sol_dict.items()
        }

        int_F_x = sol_dict['F_x']
        int_F_y = sol_dict['F_y']
        int_F_z = sol_dict['F_z']
        int_M_x = sol_dict['M_x']
        int_M_y = sol_dict['M_y']
        int_M_z = sol_dict['M_z']

        # Create `Force` and `Moment` objects:
        int_F = Force.create_from_components(
            int_F_x, int_F_y, int_F_z,
            position=Position(x, y, z, units=self._units_of_length),
            units=self._units_of_force,
        )
        int_M = Moment.create_from_components(
            int_M_x, int_M_y, int_M_z,
            position=Position(x, y, z, units=self._units_of_length),
            units=self._units_of_moment
        )
        return int_F, int_M

    def __profiles_of_internal_loadings(self) -> tuple[np.ndarray, ...]:
        # Calculates the internal forces and moments at multiple, equally spaced
        # sections and returns a Numpy array of the x-positions of these
        # sections, together with the arrays of the normal forces, the shear
        # forces and bending moments.
        x_arr = Q_(
            np.linspace(0.0, self._length, self.num_sections, endpoint=True),
            self._units_of_length
        )
        if self.num_sections <= 350:
            int_F_lst, int_M_lst = zip(*[self.cut(x) for x in x_arr])
        else:
            # if the number of sections is greater than 350, `cut` is run in
            # parallel to be faster.
            with ProcessPoolExecutor() as executor:
                int_F_lst, int_M_lst = zip(*list(executor.map(self.cut, x_arr)))
        N_lst, V_lst, _ = zip(*[int_F.component_values for int_F in int_F_lst])
        N_arr = np.array(N_lst)
        V_arr = np.array(V_lst)
        T_lst, _, M_lst = zip(*[int_M.component_values for int_M in int_M_lst])
        M_arr = np.array(M_lst)
        T_arr = np.array(T_lst)
        return x_arr.m, N_arr, V_arr, M_arr, T_arr

    @property
    def shear_diagram(self) -> LineChart:
        """Returns a `LineChart` object with the shear force diagram of the
        beam.
        """
        diagram = LineChart()
        diagram.add_xy_data(
            label='shear',
            x1_values=self._x_arr,
            y1_values=self._V_arr,
            style_props={'drawstyle': 'steps-post'}
        )
        diagram.x1.add_title(f"x, {self._units_of_length}")
        diagram.y1.add_title(f"shear force, {self._units_of_force}")
        return diagram

    @property
    def moment_diagram(self) -> LineChart:
        """Returns a `LineChart` object with the moment diagram of the
        beam.
        """
        diagram = LineChart()
        diagram.add_xy_data(
            label='bending',
            x1_values=self._x_arr,
            y1_values=self._M_arr,
            style_props={'drawstyle': 'steps-post'}
        )
        diagram.x1.add_title(f"x, {self._units_of_length}")
        diagram.y1.add_title(f"bending moment, {self._units_of_moment}")
        return diagram

    @property
    def normal_force_diagram(self) -> LineChart:
        """Returns a `LineChart` object with the normal force diagram of the
        beam.
        """
        diagram = LineChart()
        diagram.add_xy_data(
            label='normal force',
            x1_values=self._x_arr,
            y1_values=self._N_arr,
            style_props={'drawstyle': 'steps-post'}
        )
        diagram.x1.add_title(f"x, {self._units_of_length}")
        diagram.y1.add_title(f"normal force, {self._units_of_force}")
        return diagram

    @property
    def torque_diagram(self) -> LineChart:
        """Returns a `LineChart` object with the torque diagram of the beam."""
        diagram = LineChart()
        diagram.add_xy_data(
            label='torque',
            x1_values=self._x_arr,
            y1_values=self._T_arr,
            style_props={'drawstyle': 'steps-post'}
        )
        diagram.x1.add_title(f"x, {self._units_of_length}")
        diagram.y1.add_title(f"torque, {self._units_of_moment}")
        return diagram

    def V(self, x: Quantity) -> Quantity:
        """Returns the resultant internal shear force (`Quantity` object) at the
         section at position `x` (`Quantity` object) along the longitudinal axis
        of the beam (determined by linear interpolation).
        """
        V = self._V_interp(x.to(self._units_of_length).m)
        return Q_(V, self._units_of_force)

    def M(self, x: Quantity) -> Quantity:
        """Returns the resultant internal bending moment (`Quantity` object) at
        the section at position `x` (`Quantity` object) along the longitudinal
        axis of the beam (determined by linear interpolation).
        """
        M = self._M_interp(x.to(self._units_of_length).m)
        return Q_(M, self._units_of_moment)

    def N(self, x: Quantity) -> Quantity:
        """Returns the resultant internal normal force (`Quantity` object) at
        the section at position `x` (`Quantity` object) along the longitudinal
        axis of the beam (determined by linear interpolation).
        """
        N = self._N_interp(x.to(self._units_of_length).m)
        return Q_(N, self._units_of_force)

    def T(self, x: Quantity) -> Quantity:
        """Returns the resultant internal torque (`Quantity` object) at the
        section at position `x` (`Quantity` object) along the longitudinal
        axis of the beam (determined by linear interpolation).
        """
        T = self._T_interp(x.to(self._units_of_length).m)
        return Q_(T, self._units_of_moment)

    def V_max(self) -> tuple[Quantity, Quantity]:
        """Returns a 2-tuple: the first element is the x-position (`Quantity`
        object) of the section where the resultant internal shear force reaches
        a maximum (or minimum). The second element is this maximum shear force
        (`Quantity` object).
        """
        V_max = np.max(self._V_arr)
        i = np.argmax(self._V_arr)
        V_min = np.min(self._V_arr)
        if np.abs(V_min) > V_max:
            V_max = V_min
            i = np.argmin(self._V_arr)
        x = self._x_arr[i]
        return Q_(x, self._units_of_length), Q_(V_max, self._units_of_force)

    def M_max(self) -> tuple[Quantity, Quantity]:
        """Returns a 2-tuple: the first element is the x-position (`Quantity`
        object) of the section where the resultant internal bending moment
        reaches a maximum (or minimum). The second element is this maximum
        bending moment (`Quantity` object).
        """
        M_max = np.max(self._M_arr)
        i = np.argmax(self._M_arr)
        M_min = np.min(self._M_arr)
        if np.abs(M_min) > M_max:
            M_max = M_min
            i = np.argmin(self._M_arr)
        x = self._x_arr[i]
        return Q_(x, self._units_of_length), Q_(M_max, self._units_of_moment)

    def N_max(self) -> tuple[Quantity, Quantity]:
        """Returns a 2-tuple: the first element is the x-position (`Quantity`
        object) of the section where the resultant internal normal force reaches
        a maximum (or minimum). The second element is this maximum normal force
        (`Quantity` object).
        """
        N_max = np.max(self._N_arr)
        i = np.argmax(self._N_arr)
        N_min = np.min(self._N_arr)
        if np.abs(N_min) > N_max:
            N_max = N_min
            i = np.argmin(self._N_arr)
        x = self._x_arr[i]
        return Q_(x, self._units_of_length), Q_(N_max, self._units_of_force)

    def T_max(self) -> tuple[Quantity, Quantity]:
        """Returns a 2-tuple: the first element is the x-position (`Quantity`
        object) of the section where the resultant internal torque reaches
        a maximum (or minimum). The second element is this maximum torque
        (`Quantity` object).
        """
        T_max = np.max(self._T_arr)
        i = np.argmax(self._T_arr)
        T_min = np.min(self._T_arr)
        if np.abs(T_min) > T_max:
            T_max = T_min
            i = np.argmin(self._T_arr)
        x = self._x_arr[i]
        return Q_(x, self._units_of_length), Q_(T_max, self._units_of_moment)

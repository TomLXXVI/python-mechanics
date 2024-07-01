import sympy as sp
from mechanics import Quantity
from .force import Force, Moment, Torque
from .distributed_load import DistributedLoad1D
from .units import Units

Q_ = Quantity


class System:
    
    def __init__(
        self,
        *elements: Force | Torque | DistributedLoad1D
    ) -> None:
        self.forces = [f for f in elements if isinstance(f, Force)]
        self.torques = [t for t in elements if isinstance(t, (Torque, Moment))]
        self.forces.extend([dl.resultant() for dl in elements if isinstance(dl, DistributedLoad1D)])
        self.__set_units()
        self._undetermined_forces: dict[str, Force] = {}
        self._undetermined_torques: dict[str, Torque] = {}

    def solve(self) -> dict[str, Force | Torque]:
        dicts = self.__separate()
        eqs = self.__create_equations(dicts)
        sol_dict = sp.solve(eqs, dict=True)
        if isinstance(sol_dict, list):
            sol_dict = sol_dict[0]
        solutions = self.__get_solutions(sol_dict)
        return solutions

    def __separate(self) -> tuple[dict, ...]:
        F_x_dict = {'unknown': [], 'known': []}
        F_y_dict = {'unknown': [], 'known': []}
        F_z_dict = {'unknown': [], 'known': []}
        M_x_dict = {'unknown': [], 'known': []}
        M_y_dict = {'unknown': [], 'known': []}
        M_z_dict = {'unknown': [], 'known': []}
        for force in self.forces:
            F_x, F_y, F_z = force.components
            M_x, M_y, M_z = force.moment().components
            if isinstance(F_x, sp.Expr):
                F_x_dict['unknown'].append(F_x)
                self._undetermined_forces[force.name] = force
            else:
                F_x_dict['known'].append(F_x.to(self.__u_force).m)
            if isinstance(F_y, sp.Expr):
                F_y_dict['unknown'].append(F_y)
                self._undetermined_forces[force.name] = force
            else:
                F_y_dict['known'].append(F_y.to(self.__u_force).m)
            if isinstance(F_z, sp.Expr):
                F_z_dict['unknown'].append(F_z)
                self._undetermined_forces[force.name] = force
            else:
                F_z_dict['known'].append(F_z.to(self.__u_force).m)
            if isinstance(M_x, sp.Expr):
                M_x_dict['unknown'].append(M_x)
            else:
                M_x_dict['known'].append(M_x.to(self.__u_moment).m)
            if isinstance(M_y, sp.Expr):
                M_y_dict['unknown'].append(M_y)
            else:
                M_y_dict['known'].append(M_y.to(self.__u_moment).m)
            if isinstance(M_z, sp.Expr):
                M_z_dict['unknown'].append(M_z)
            else:
                M_z_dict['known'].append(M_z.to(self.__u_moment).m)
        for torque in self.torques:
            M_x, M_y, M_z = torque.components
            if isinstance(M_x, sp.Expr):
                M_x_dict['unknown'].append(M_x)
                self._undetermined_torques[torque.name] = torque
            else:
                M_x_dict['known'].append(M_x.to(self.__u_moment).m)
            if isinstance(M_y, sp.Expr):
                M_y_dict['unknown'].append(M_y)
                self._undetermined_torques[torque.name] = torque
            else:
                M_y_dict['known'].append(M_y.to(self.__u_moment).m)
            if isinstance(M_z, sp.Expr):
                M_z_dict['unknown'].append(M_z)
                self._undetermined_torques[torque.name] = torque
            else:
                M_z_dict['known'].append(M_z.to(self.__u_moment).m)
        return (
            F_x_dict,
            F_y_dict,
            F_z_dict,
            M_x_dict,
            M_y_dict,
            M_z_dict
        )

    def __set_units(self) -> None:
        self.__u_pos = Units.u_pos
        self.__u_force = Units.u_force
        self.__u_moment = Units.u_moment

    @staticmethod
    def __create_equation(d: dict[str, list]) -> sp.Eq | None:
        unknown = sum(expr for expr in d['unknown'])
        if unknown == 0:
            return None
        else:
            known = sum(known for known in d['known'])
            equation = sp.Eq(unknown + known, 0)
            return equation

    def __create_equations(
        self,
        dicts: tuple[dict[str, list], ...]
    ) -> list[sp.Eq]:
        eqs = []
        for d in dicts:
            eq = self.__create_equation(d)
            if eq is not None:
                eqs.append(eq)
        return eqs

    def __get_solutions(self, sol_dict) -> dict[str, Force | Torque]:
        sol_dict = {str(sym): val for sym, val in sol_dict.items()}
        solutions = {}
        for key, val in sol_dict.items():
            s = key.split('.')
            val = float(val)
            try:
                name, component = s[0], s[1]
            except IndexError:
                name = s[0]
                solutions[name] = {'magnitude': float(val)}
            else:
                if name in solutions.keys():
                    solutions[name][component] = val
                else:
                    solutions[name] = {component: val}
        for name, component_dict in solutions.items():
            try:
                unknown_force = self._undetermined_forces[name]
            except KeyError:
                unknown_torque = self._undetermined_torques[name]
                solved_torque = self.__create_solved_torque(name, component_dict, unknown_torque)
                solutions[name] = solved_torque
            else:
                solved_force = self.__create_solved_force(name, component_dict, unknown_force)
                solutions[name] = solved_force
        return solutions

    def __create_solved_force(
        self,
        name: str,
        component_dict: dict[str, float],
        unknown_force: Force
    ) -> Force:
        F_mag = component_dict.get('magnitude', None)
        if F_mag is not None:
            if F_mag < 0:
                F_mag = abs(F_mag)
                theta = unknown_force.theta + Q_(180, 'deg')
                if unknown_force.gamma.m != 0:
                    gamma = unknown_force.gamma + Q_(180, 'deg')
                else:
                    gamma = unknown_force.gamma
            else:
                theta = unknown_force.theta
                gamma = unknown_force.gamma
            solved_force = Force(
                action_point=unknown_force.action_point,
                magnitude=Q_(F_mag, self.__u_force),
                theta=theta,
                gamma=gamma,
                name=name
            )
            return solved_force
        else:
            F_x = component_dict.get('x', 0.0)
            F_y = component_dict.get('y', 0.0)
            F_z = component_dict.get('z', 0.0)
            solved_force = Force.create_from_components(
                action_point=unknown_force.action_point,
                F_x=Q_(F_x, self.__u_force),
                F_y=Q_(F_y, self.__u_force),
                F_z=Q_(F_z, self.__u_force),
                name=name
            )
            return solved_force

    def __create_solved_torque(
        self,
        name: str,
        component_dict: dict[str, float],
        unknown_torque: Torque
    ) -> Torque:
        M_mag = component_dict.get('magnitude', None)
        if M_mag is not None:
            if M_mag < 0:
                M_mag = abs(M_mag)
                theta = unknown_torque.theta + Q_(180, 'deg')
                if unknown_torque.gamma.m != 0:
                    gamma = unknown_torque.gamma + Q_(180, 'deg')
                else:
                    gamma = unknown_torque.gamma
            else:
                theta = unknown_torque.theta
                gamma = unknown_torque.gamma
            solved_torque = Torque(
                action_point=unknown_torque.action_point,
                magnitude=Q_(M_mag, self.__u_moment),
                theta=theta,
                gamma=gamma,
                name=name
            )
            return solved_torque
        else:
            M_x = component_dict.get('x', 0.0)
            M_y = component_dict.get('y', 0.0)
            M_z = component_dict.get('z', 0.0)
            solved_torque = Torque.create_from_components(
                action_point=unknown_torque.action_point,
                M_x=Q_(M_x, self.__u_moment),
                M_y=Q_(M_y, self.__u_moment),
                M_z=Q_(M_z, self.__u_moment),
                name=name
            )
            return solved_torque

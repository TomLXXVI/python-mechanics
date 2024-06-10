import sympy as sp
from mechanics import Quantity
from .force import Force, Torque
from .distributed_load import DistributedLoad2D


Q_ = Quantity


class System:
    
    def __init__(
        self,
        forces: list[Force] | None = None,
        torques: list[Torque] | None = None,
        distributed_loads: list[DistributedLoad2D] | None = None
    ) -> None:
        self.forces = forces if forces is not None else []
        self.torques = torques if torques is not None else []
        if distributed_loads is not None:
            self.add_distributed_loads(*distributed_loads)
        self._undetermined_forces: dict[str, Force] = {}
        self._undetermined_torques: dict[str, Torque] = {}
    
    def add_forces(self, *forces) -> None:
        self.forces.extend(forces)

    def add_torques(self, *torques) -> None:
        self.torques.extend(torques)

    def add_distributed_loads(self, *distributed_loads) -> None:
        self.forces.extend([dl.resultant() for dl in distributed_loads])

    def solve(self) -> tuple[dict[str, Force], dict[str, Torque]]:
        dicts = self.__separate()
        eqs = self.__create_equations(dicts)
        sol_dict = sp.solve(eqs, dict=True)
        if isinstance(sol_dict, list):
            sol_dict = sol_dict[0]
        solved_forces, solved_torques = self.__get_solutions(sol_dict)
        return solved_forces, solved_torques

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
                F_x_dict['known'].append(F_x.to('N').m)
            if isinstance(F_y, sp.Expr):
                F_y_dict['unknown'].append(F_y)
                self._undetermined_forces[force.name] = force
            else:
                F_y_dict['known'].append(F_y.to('N').m)
            if isinstance(F_z, sp.Expr):
                F_z_dict['unknown'].append(F_z)
                self._undetermined_forces[force.name] = force
            else:
                F_z_dict['known'].append(F_z.to('N').m)
            if isinstance(M_x, sp.Expr):
                M_x_dict['unknown'].append(M_x)
            else:
                M_x_dict['known'].append(M_x.to('N * m').m)
            if isinstance(M_y, sp.Expr):
                M_y_dict['unknown'].append(M_y)
            else:
                M_y_dict['known'].append(M_y.to('N * m').m)
            if isinstance(M_z, sp.Expr):
                M_z_dict['unknown'].append(M_z)
            else:
                M_z_dict['known'].append(M_z.to('N * m').m)
        for torque in self.torques:
            M_x, M_y, M_z = torque.components
            if isinstance(M_x, sp.Expr):
                M_x_dict['unknown'].append(M_x)
                self._undetermined_torques[torque.name] = torque
            else:
                M_x_dict['known'].append(M_x.to('N * m').m)
            if isinstance(M_y, sp.Expr):
                M_y_dict['unknown'].append(M_y)
                self._undetermined_torques[torque.name] = torque
            else:
                M_y_dict['known'].append(M_y.to('N * m').m)
            if isinstance(M_z, sp.Expr):
                M_z_dict['unknown'].append(M_z)
                self._undetermined_torques[torque.name] = torque
            else:
                M_z_dict['known'].append(M_z.to('N * m').m)
        return (
            F_x_dict,
            F_y_dict,
            F_z_dict,
            M_x_dict,
            M_y_dict,
            M_z_dict
        )

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

    def __get_solutions(self, sol_dict) -> tuple[dict[str, Force], dict[str, Torque]]:
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
        solved_forces, solved_torques = {}, {}
        for name, component_dict in solutions.items():
            try:
                unknown_force = self._undetermined_forces[name]
            except KeyError:
                unknown_torque = self._undetermined_torques[name]
                solved_torque = self.__create_solved_torque(name, component_dict, unknown_torque)
                solved_torques[name] = solved_torque
            else:
                solved_force = self.__create_solved_force(name, component_dict, unknown_force)
                solved_forces[name] = solved_force
        return solved_forces, solved_torques

    @staticmethod
    def __create_solved_force(
        name: str,
        component_dict: dict[str, float],
        unknown_force: Force
    ) -> Force:
        F_mag = component_dict.get('magnitude', None)
        if F_mag is not None:
            solved_force = Force(
                action_point=unknown_force.action_point,
                magnitude=Q_(F_mag, 'N'),
                theta=unknown_force.theta,
                gamma=unknown_force.gamma,
                name=name
            )
            return solved_force
        else:
            F_x = component_dict.get('x', 0.0)
            F_y = component_dict.get('y', 0.0)
            F_z = component_dict.get('z', 0.0)
            solved_force = Force.create_from_components(
                action_point=unknown_force.action_point,
                F_x=Q_(F_x, 'N'),
                F_y=Q_(F_y, 'N'),
                F_z=Q_(F_z, 'N'),
                name=name
            )
            return solved_force

    @staticmethod
    def __create_solved_torque(
        name: str,
        component_dict: dict[str, float],
        unknown_torque: Torque
    ) -> Torque:
        M_mag = component_dict.get('magnitude', None)
        if M_mag is not None:
            solved_torque = Torque(
                action_point=unknown_torque.action_point,
                magnitude=Q_(M_mag, 'N * m'),
                theta=unknown_torque.theta,
                gamma=unknown_torque.gamma,
                name=name
            )
            return solved_torque
        else:
            M_x = component_dict.get('x', 0.0)
            M_y = component_dict.get('y', 0.0)
            M_z = component_dict.get('z', 0.0)
            solved_torque = Torque.create_from_components(
                action_point=unknown_torque.action_point,
                M_x=Q_(M_x, 'N * m'),
                M_y=Q_(M_y, 'N * m'),
                M_z=Q_(M_z, 'N * m'),
                name=name
            )
            return solved_torque

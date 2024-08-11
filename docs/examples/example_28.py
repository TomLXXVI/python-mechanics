# Demo based on example 7.1 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.

# Part 2: Calculate shear stresses in T-shaped section

import numpy as np
from mechanics import Quantity
from mechanics.charts import LineChart
from mechanics.strength.stress import TransverseShear
from example_27 import TShape


Q_ = Quantity


if __name__ == '__main__':

    # Create the shape of the beam's cross-section
    t_shape = TShape(
        width=Q_(150, 'mm'),
        height=Q_(180, 'mm'),
        web_thickness=Q_(30, 'mm'),
        flange_thickness=Q_(30, 'mm')
    )
    t_shape.plot().show()

    # Apply the resultant internal shear force to the section:
    shear = TransverseShear(V=Q_(-19.5, 'kN'), shape=t_shape)

    # Get the shear stress at each vertex of the shape:
    for i, vertex in enumerate(t_shape.vertices):
        print(
            f"vertex {i} @ "
            f"(z = {vertex[0]:~P.0f}, y = {vertex[1]:~P.0f}): "
            f"tau = {shear.tau(vertex[1]).to('MPa'):~P.2f}"
        )

    # Get profile of shear stress in the section:
    y_coords = Q_(np.linspace(-120, 60, endpoint=True), 'mm')
    tau_arr = Quantity.from_list([shear.tau(y).to('MPa') for y in y_coords])
    profile = LineChart()
    profile.add_xy_data(
        label='tau',
        x1_values=y_coords.m,
        y1_values=np.abs(tau_arr.m),
        style_props={'drawstyle': 'steps-post'}
    )
    profile.x1.add_title('y-coordinate, mm')
    profile.y1.add_title('tau, MPa')
    profile.show()

    # Get shear stress at given y-coordinate:
    y = Q_(45, 'mm')
    tau = shear.tau(y).to('MPa')
    print(f"tau @ y {y:~P.0f} = {tau:~P.2f}")

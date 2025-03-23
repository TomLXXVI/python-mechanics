# Demo based on example 7.3 from Hibbeler, R. C. (2017). Mechanics of Materials
# in SI Units, 10th Edition.

import numpy as np
from mechanics import Quantity
from mechanics.charts import LineChart
from mechanics.strength import TransverseShear
from mechanics.geometry import HShape


Q_ = Quantity


# Create the shape of the beam's cross-section:
h_shape = HShape(
    height=Q_(240, 'mm'),
    width=Q_(300, 'mm'),
    web_thickness=Q_(15, 'mm'),
    flange_thickness=Q_(20, 'mm')
)
h_shape.plot().show()


# Apply the resultant internal shear force to the section:
shear = TransverseShear(V=Q_(-80, 'kN'), shape=h_shape)


# Get profile of shear stress in the section:
y_coords = Q_(np.linspace(-120, 120, endpoint=True), 'mm')
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

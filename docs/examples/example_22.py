# Demo of class `CShape` (see section.geometry.polygons.py)
from mechanics import Quantity
from mechanics.geometry import CShape

Q_ = Quantity


# Draw C-shape:

channel1 = CShape(
    height=Q_(280, 'mm'),
    width=Q_(200, 'mm'),
    web_thickness=Q_(20, 'mm'),
    flange_thickness=Q_(15, 'mm')
)

channel1.plot().show()


# Draw C-shape turned by 90° degrees clockwise:

channel2 = CShape(
    height=Q_(280, 'mm'),
    width=Q_(200, 'mm'),
    web_thickness=Q_(20, 'mm'),
    flange_thickness=Q_(15, 'mm'),
    orientation=Q_(-90, 'deg')
)

channel2.plot().show()

"""
A collection of useful functions to make live easier...
"""
import numpy as np
from mechanics import Quantity

Q_ = Quantity


def create_angle(
    v: float | Quantity,
    h: float | Quantity,
    quadrant: int = 1
) -> Quantity:
    """Creates an angle (instance of class `Quantity`) given the slope
    specified by a vertical distance `v` and horizontal distance `h`. By
    specifying the quadrant, the slope angle can be positioned in a right-handed
    (x,y) coordinate system (the x-axis runs horizontally pointing to the right
    and the y-axis runs vertically pointing upward). The first quadrant is
    situated between 0° and 90° (where both the x- and y-coordinates have
    positive values). The second quadrant is situated between 90° and 180°
    (where the x-coordinates have a negative value and the y-coordinates have a
    positive value). The third quadrant is situated between 180° and 270° (where
    both the x- and y-coordinates have negative values). The fourth quadrant is
    situated between 270° and 360° (where the x-coordinates have a positive
    value and the y-coordinates have a negative value).
    """
    if isinstance(v, Quantity) and isinstance(h, Quantity):
        v = v.m
        h = h.to(v.units).m
    a = Q_(np.arctan2(v, h), 'rad')
    match quadrant:
        case 1:
            return a
        case 2:
            return np.pi - a
        case 3:
            return a + np.pi
        case 4:
            return -a

from mechanics import Quantity

Q_ = Quantity


Position = Quantity | tuple[Quantity, ...]

ORIGIN = (Q_(0, 'm'), Q_(0, 'm'), Q_(0, 'm'))


def set_position(position: Position) -> tuple[Quantity, Quantity, Quantity]:
    """Returns a position in 2D- or 3D-space as tuple of 3 `Quantity`
    objects with the same units of length. The first element specifies the
    x-coordinate of the position, the second element the y-coordinate, and the
    third element the z-coordinate.

    Parameter `position` can be specified in multiple ways:
    1.  as a single `Quantity` object with units of length which indicates a
        position on the x-axis of the coordinate system.
    2.  as a tuple containing only 1 `Quantity` object, also indicating a
        position on the x-axis of the coordinate system.
    3.  as a `Quantity` object of the form `Q_([x, y], <length_unit>)` which
        indicates a position in the xy-plane of the coordinate system.
    4.  as a tuple containing 2 `Quantity` objects with units of length, also
        indicating a position in the xy-plane of the coordinate system.
    5.  as a `Quantity` object of the form `Q_([x, y, z], <length_unit>)` which
        indicates a position in the xyz-space of the coordinate system.
    6.  as tuple containing 3 `Quantity` objects with units of length, also
        indicating a position in the xyz-space of the coordinate system.
    """
    if isinstance(position, tuple):
        if len(position) == 1:
            pos = (position[0], Q_(0, position[0].units), Q_(0, position[0].units))
        elif len(position) == 2:
            pos = (position[0], position[1], Q_(0, position[0].units))
        elif len(position) == 3:
            pos = position
        else:
            raise ValueError("too many coordinates (maximum is 3: x, y, and z)")
    elif isinstance(position, Quantity):
        try:
            _ = position.size
        except AttributeError:
            pos = (position, Q_(0, position.units), Q_(0, position.units))
        else:
            if position.size == 1:
                try:
                    pos = (position[0], Q_(0, position[0].units), Q_(0, position[0].units))
                except IndexError:
                    pos = (position, Q_(0, position.units), Q_(0, position.units))
            elif position.size == 2:
                pos = (position[0], position[1], Q_(0, position[0].units))
            elif position.size == 3:
                pos = (position[0], position[1], position[2])
            else:
                raise ValueError("too many coordinates (maximum is 3: x, y, and z)")
    else:
        raise ValueError("not a valid position value")
    return pos

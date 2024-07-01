from typing import Type
from .shapes import Dimensions, Shape, Circle, Annulus
from .polygons import Rectangle, HShape, CShape, ZShape
from .hollow_polygons import HollowRectangle


def create_geometry(shape: Type[Shape], dim: Dimensions) -> Shape:
    if shape is Rectangle:
        shape_obj = Rectangle(
            width=dim.width,
            height=dim.height
        )
        return shape_obj
    elif shape is Circle:
        shape_obj = Circle(
            radius=dim.radius or dim.diameter / 2
        )
        return shape_obj
    elif shape is Annulus:
        shape_obj = Annulus(
            radius_in=dim.inner_radius,
            radius_out=dim.outer_radius
        )
        return shape_obj
    elif shape is HShape:
        shape_obj = HShape(
            height=dim.height,
            width=dim.width,
            web_thickness=dim.web_thickness,
            flange_thickness=dim.flange_thickness
        )
        return shape_obj
    elif shape is CShape:
        shape_obj = CShape(
            height=dim.height,
            width=dim.width,
            web_thickness=dim.web_thickness,
            flange_thickness=dim.flange_thickness,
            orientation=dim.orientation
        )
        return shape_obj
    elif shape is ZShape:
        shape_obj = ZShape(
            height=dim.height,
            width=dim.width,
            web_thickness=dim.web_thickness,
            flange_thickness=dim.flange_thickness,
            orientation=dim.orientation
        )
        return shape_obj
    elif shape is HollowRectangle:
        shape_obj = HollowRectangle(
            width=dim.width,
            height=dim.height,
            thickness=dim.thickness
        )
        return shape_obj

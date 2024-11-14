#!/usr/bin/env python

from sympy import pi, sin, cos, tan, Float

"""
@file rough_overlap_check.py
@brief Script to generate an SVG file containing multiple polygons and circles.
@detail This script calculates the vertices of polygons and parameters of circles symbolically,
        and renders them as an SVG file using sympy for symbolic computation.
@version 1.0.0
@date 2024-11-14
@author Takumi Shiota
"""

__version__ = "1.0.0"
__author__ = "Takumi Shiota"
__date__ = "2024-11-14"

# SVG templates
SVG_HEADER_TEMPLATE = """<?xml version="1.0" encoding="utf-8"?>
<svg xmlns="http://www.w3.org/2000/svg" version="1.1" viewBox="{view_box}">
<style>
    .no_fill_black_stroke {{ fill: none; stroke: #000000; stroke-width: 0.03; }}
    .circle_style {{ fill: none; stroke: #FF0000; stroke-width: 0.03; }}
</style>
"""

SVG_FOOTER = """</svg>
"""


class Polygon:
    """
    @class Polygon
    @brief Represents a polygon and calculates its vertices symbolically.
    """
    def __init__(self, n, side_length=1):
        """
        @brief Initializes a polygon with the given number of sides and side length.
        @param n Number of sides (must be >= 3).
        @param side_length Length of each side.
        @exception ValueError Raised if n < 3.
        """
        if n < 3:
            raise ValueError("A polygon must have at least 3 sides.")
        self.n = n
        self.side_length = side_length

    def calculate_polygon_vertices(self, cx, cy):
        """
        @brief Calculates the vertices of a regular polygon.
        @param cx X-coordinate of the polygon's center.
        @param cy Y-coordinate of the polygon's center.
        @return List of symbolic (x, y) tuples representing the vertices.
        """
        angle_step = 2 * pi / self.n
        radius = self.side_length / (2 * sin(pi / self.n))
        offset_angle = angle_step / 2

        vertices = [
            (cx + radius * cos(offset_angle + i * angle_step),
             cy + radius * sin(offset_angle + i * angle_step))
            for i in range(self.n)
        ]

        return vertices


class Circle:
    """
    @class Circle
    @brief Represents a circle and calculates its properties symbolically.
    """
    def __init__(self, cx, cy, radius):
        """
        @brief Initializes a circle with the given center and radius.
        @param cx X-coordinate of the circle's center.
        @param cy Y-coordinate of the circle's center.
        @param radius Radius of the circle.
        """
        self.cx = cx
        self.cy = cy
        self.radius = radius

    def get_parameters(self):
        """
        @brief Returns the parameters of the circle.
        @return Tuple containing (cx, cy, radius).
        """
        return (self.cx, self.cy, self.radius)


class SvgDrawer:
    """
    @class SvgDrawer
    @brief Creates and renders SVG files with multiple polygons and circles.
    """
    def __init__(self, filename):
        self.filename = filename

    def calculate_viewbox(self, shapes, margin=0.2, precision=100):
        """
        @brief Calculates a dynamic viewBox based on the given shapes.
        @param shapes List of shapes (polygons and circles).
        @param margin Extra space to include around the shapes.
        @param precision Number of decimal digits for evaluation.
        @return A formatted string for the viewBox.
        """
        x_coords = []
        y_coords = []

        for shape in shapes:
            if isinstance(shape, list):  # Assuming list of vertices for polygons
                for x, y in shape:
                    x_num = x.evalf(precision)
                    y_num = y.evalf(precision)
                    x_coords.append(float(x_num))
                    y_coords.append(float(y_num))
            elif isinstance(shape, Circle):
                cx_num = shape.cx.evalf(precision)
                cy_num = shape.cy.evalf(precision)
                r_num = shape.radius.evalf(precision)
                x_coords.extend([float(cx_num - r_num), float(cx_num + r_num)])
                y_coords.extend([float(cy_num - r_num), float(cy_num + r_num)])

        min_x = min(x_coords) - float(margin)
        max_x = max(x_coords) + float(margin)
        min_y = min(y_coords) - float(margin)
        max_y = max(y_coords) + float(margin)

        width = max_x - min_x
        height = max_y - min_y

        return f"{min_x} {min_y} {width} {height}"

    def draw_shapes(self, shapes, precision=100):
        """
        @brief Generates an SVG file containing the given shapes.
        @param shapes List containing polygons (as lists of vertices) and Circle objects.
        @param precision Number of decimal digits for evaluation.
        """
        view_box = self.calculate_viewbox(shapes, precision=precision)

        with open(self.filename, 'w') as file:
            file.write(SVG_HEADER_TEMPLATE.format(view_box=view_box))

            for shape in shapes:
                if isinstance(shape, list):  # Assuming list of vertices for polygons
                    evaluated_polygon = [
                        (x.evalf(precision), y.evalf(precision))
                        for x, y in shape
                    ]
                    points = " ".join(f"{x},{y}" for x, y in evaluated_polygon)
                    file.write(f"""<polygon class="no_fill_black_stroke" points="{points}" />\n""")
                elif isinstance(shape, Circle):
                    cx, cy, r = shape.get_parameters()
                    cx_num = cx.evalf(precision)
                    cy_num = cy.evalf(precision)
                    r_num = r.evalf(precision)
                    file.write(f"""<circle class="circle_style" cx="{cx_num}" cy="{cy_num}" r="{r_num}" />\n""")

            file.write(SVG_FOOTER)


def main():
    """
    @brief Main function to execute the script.
    """
    output_filename = "out.svg"
    n = int(input("Enter the number of sides for the polygon (n >= 3): "))

    if n < 3:
        print("Error: A polygon must have at least 3 sides.")
        return
    
    a = int(input("a: "))
    if a < 0:
        return
    
    k = float(input("k: "))
    if k <= 0:
        return
    
    c1_pos = float(input("c1_pos: "))
    if c1_pos <= a:
        return

    cx1 = Float(0)
    cy1 = Float(0)

    # Create a polygon and calculate its vertices
    polygon = Polygon(n)
    polygon_vertices = polygon.calculate_polygon_vertices(cx=cx1, cy=cy1)

    # Calculate the incircle of the polygon
    incircle_radius = Float(1) / (2 * tan(pi / n))  # Radius of the incircle
    circle1 = Circle(cx=cx1, cy=cy1, radius=incircle_radius)
    rectangle_a = [
        (incircle_radius, Float(-0.5)),
        (incircle_radius, Float(-0.5) + a),
        (incircle_radius + k, Float(-0.5) + a),
        (incircle_radius + k, Float(-0.5))
    ]

    cx2 = cx1 + Float(2) * incircle_radius + k
    cy2 = a - Float(1)

    circle2 = Circle(cx=cx2, cy=cy2, radius=incircle_radius)

    th = (c1_pos * Float(2) * pi) / n

    cx3 = cx2 - incircle_radius * cos(th)
    cy3 = cy2 - incircle_radius * sin(th)

    circle3 = Circle(cx=cx3, cy=cy3, radius=Float(0.01))

    cx4 = cx3 - (c1_pos - a + Float(0.5)) * sin(th)
    cy4 = cy3 + (c1_pos - a + Float(0.5)) * cos(th)

    cx5 = cx4 - k * cos(th)
    cy5 = cy4 - k * sin(th)

    cx6 = cx3 + Float(0.5) * sin(th)
    cy6 = cy3 - Float(0.5) * cos(th)

    cx7 = cx6 - k * cos(th)
    cy7 = cy6 - k * sin(th)

    rectangle_c1 = [
        (cx6, cy6),
        (cx4, cy4),
        (cx5, cy5),
        (cx7, cy7)
    ]

    # Prepare shapes to draw
    shapes = [polygon_vertices, circle1, rectangle_a, circle2, circle3, rectangle_c1]

    # Draw the shapes
    svg_drawer = SvgDrawer(output_filename)
    svg_drawer.draw_shapes(shapes)

    # Calculate the squared distance from the origin
    distance_squared1 = (cx4**2 + cy4**2).evalf(100)
    distance_squared2 = (cx5**2 + cy5**2).evalf(100)

    # Compare with the squared radius
    if distance_squared1 < incircle_radius**2 or distance_squared2 < incircle_radius**2:
        print("Included")
    else:
        print("Not included")

####################
if __name__ == "__main__":
    main()

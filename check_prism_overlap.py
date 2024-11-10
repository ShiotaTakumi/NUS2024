#! /usr/bin/env python

from sympy import pi, sin, cos

"""
@file prism_overlap.py
@brief Script to generate an SVG file containing multiple polygons and a rectangle.
@detail This script calculates the vertices of a polygon (e.g., prism base) and a rectangle symbolically,
        and renders them as an SVG file.
@version 1.0.6
@date 2024-11-10
@author Takumi Shiota
"""

__version__ = "1.0.6"
__author__ = "Takumi Shiota"
__date__ = "2024-11-10"

# SVG templates
SVG_HEADER_TEMPLATE = """<?xml version="1.0" encoding="utf-8"?>
<svg xmlns="http://www.w3.org/2000/svg" version="1.1" viewBox="{view_box}">
<style>.no_fill_black_stroke {{fill:none; stroke:#000000; stroke-width:0.01;}}</style>
"""

SVG_FOOTER = """</svg>
"""


class Polygon:
    """
    @class Polygon
    @brief Represents a polygon and provides methods to calculate its vertices and related shapes.
    """
    def __init__(self, n, height, a_height, side_length=1):
        """
        @brief Initializes a polygon and related rectangle parameters.
        @param n Number of sides of the polygon (must be >= 3).
        @param height Height of the prism (used for rectangle width).
        @param a_height Height of Rectangle A.
        @param side_length Length of each side of the polygon (default: 1).
        @exception ValueError Raised if n < 3.
        """
        if n < 3:
            raise ValueError("A polygon must have at least 3 sides.")
        self.n = n
        self.side_length = side_length
        self.height = height
        self.a_height = a_height

    def calculate_b1_vertices(self, cx, cy):
        """
        @brief Calculates the vertices of a regular polygon (base of prism).
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

    def calculate_a_vertices(self, start):
        """
        @brief Calculates the vertices of Rectangle A.
        @param start Starting vertex (x, y) of the rectangle.
        @return List of (x, y) tuples representing the vertices of the rectangle.
        """
        width = self.height  # Rectangle width (based on prism height)
        height = self.a_height  # Rectangle height (user-defined)
        vertices = [
            start,
            (start[0] + width, start[1]),
            (start[0] + width, start[1] + height),
            (start[0], start[1] + height)
        ]
        return vertices


class SvgDrawer:
    """
    @class SvgDrawer
    @brief Creates and renders SVG files containing multiple polygons.
    """
    def __init__(self, filename):
        """
        @brief Initializes the SVG drawer with the output filename.
        @param filename Name of the output SVG file.
        """
        self.filename = filename

    def calculate_viewbox(self, polygons, margin=0.2, precision=100):
        """
        @brief Calculates a dynamic viewBox based on the given polygons.
        @param polygons List of lists of symbolic (x, y) tuples defining the polygons.
        @param margin Extra space to include around the polygons.
        @param precision Number of decimal digits for evaluation.
        @return A formatted string for the viewBox.
        """
        evaluated_vertices = []
        for polygon in polygons:
            for x, y in polygon:
                x_num = x.evalf(precision)
                y_num = y.evalf(precision)
                evaluated_vertices.append((float(x_num), float(y_num)))

        x_coords = [x for x, _ in evaluated_vertices]
        y_coords = [y for _, y in evaluated_vertices]

        min_x = min(x_coords) - float(margin)
        max_x = max(x_coords) + float(margin)
        min_y = min(y_coords) - float(margin)
        max_y = max(y_coords) + float(margin)

        width = max_x - min_x
        height = max_y - min_y

        return f"{min_x} {min_y} {width} {height}"

    def draw_unfolding(self, polygons, precision=100):
        """
        @brief Generates an SVG file containing the given polygons.
        @param polygons List of lists of symbolic (x, y) tuples defining the polygons.
        @param precision Number of decimal digits for evaluation.
        """
        view_box = self.calculate_viewbox(polygons, precision=precision)

        with open(self.filename, 'w') as file:
            file.write(SVG_HEADER_TEMPLATE.format(view_box=view_box))

            for polygon in polygons:
                evaluated_polygon = [
                    (x.evalf(precision), y.evalf(precision))
                    for x, y in polygon
                ]
                points = " ".join(f"{x},{y}" for x, y in evaluated_polygon)
                file.write(f"""<polygon class="no_fill_black_stroke" points="{points}" />\n""")

            file.write(SVG_FOOTER)


def main():
    """
    @brief Main function to execute the script.
    """
    output_filename = "out.svg"

    # Input number of sides for the polygon
    n = int(input(f"Enter the number of sides for the polygon (n >= 3): "))
    if n < 3:
        print(f"Error: A polygon must have at least 3 sides.")
        return

    # Input the height of the prism
    height = float(input(f"Enter the height of the prism (e.g., 1.5): "))
    if height <= 0:
        print(f"Error: The height of the prism must be greater than 0.")
        return

    # Input the height of Rectangle A
    a_height = int(input(f"Enter the height of Rectangle A (n >= a_len): "))
    if a_height > n:
        print(f"Error: Height of Rectangle A must not exceed {n}.")
        return

    # Create polygon and calculate vertices
    polygon = Polygon(n, height, a_height)
    b1_vertices = polygon.calculate_b1_vertices(cx=0, cy=0)  # Vertices of the polygon (base of prism)
    a_vertices = polygon.calculate_a_vertices(start=b1_vertices[n-1])  # Rectangle A vertices

    # Render SVG with the calculated shapes
    svg_drawer = SvgDrawer(output_filename)
    svg_drawer.draw_unfolding([b1_vertices, a_vertices])


if __name__ == "__main__":
    main()

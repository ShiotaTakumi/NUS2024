#! /usr/bin/env python

from sympy import Float, pi, sin, cos

"""
@file prism_overlap.py
@brief Script to generate an SVG file containing multiple polygons.
@detail This script calculates the vertices of polygons and renders them in an SVG file.
@version 1.0.4
@date 2024-11-10
@author Takumi Shiota
"""

__version__ = "1.0.4"
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
    @brief Represents a polygon and calculates its vertices.
    """
    def __init__(self, n, side_length=1, precision=100):
        """
        @brief Initializes a polygon with the given number of sides and side length.
        @param n Number of sides (must be >= 3).
        @param side_length Length of each side.
        @param precision Number of decimal digits for calculations (default: 100).
        @exception ValueError Raised if n < 3.
        """
        if n < 3:
            raise ValueError("A polygon must have at least 3 sides.")
        self.n = n
        self.side_length = Float(side_length, precision)
        self.precision = precision

    def calculate_polygon_vertices(self, cx, cy):
        """
        @brief Calculates the vertices of a regular polygon centered at (cx, cy).
        @param cx X-coordinate of the center.
        @param cy Y-coordinate of the center.
        @return List of (x, y) tuples representing the vertices of the polygon.
        """
        cx = Float(cx, self.precision)
        cy = Float(cy, self.precision)

        # Angle between consecutive vertices (in radians)
        angle_step = Float(2, self.precision) * pi / Float(self.n, self.precision)

        # Radius of the circumcircle of the polygon
        radius = self.side_length / (Float(2, self.precision) * sin(pi / Float(self.n, self.precision)).evalf(self.precision))

        # Offset to rotate the polygon so the first vertex aligns properly
        offset_angle = angle_step / Float(2, self.precision)

        # Calculate vertices in a counterclockwise order
        vertices = [
            (cx + radius * cos(offset_angle + i * angle_step).evalf(self.precision),
             cy + radius * sin(offset_angle + i * angle_step).evalf(self.precision))
            for i in range(self.n)
        ]
        return vertices


class SvgDrawer:
    """
    @class SvgDrawer
    @brief Creates and renders SVG files with multiple polygons.
    """
    def __init__(self, filename):
        """
        @brief Initializes the SVG drawer with the output filename.
        @param filename Name of the output SVG file.
        """
        self.filename = filename

    def calculate_viewbox(self, polygons, margin=Float(0.2, 100)):
        """
        @brief Calculates a dynamic viewBox based on the given polygons.
        @param polygons List of lists of (x, y) tuples defining the polygons.
        @param margin Extra space to include around the polygons.
        @return A formatted string for the viewBox.
        """
        all_vertices = [vertex for polygon in polygons for vertex in polygon]
        min_x = min(x for x, _ in all_vertices)
        max_x = max(x for x, _ in all_vertices)
        min_y = min(y for _, y in all_vertices)
        max_y = max(y for _, y in all_vertices)

        # Add margin to the boundaries
        min_x -= margin
        max_x += margin
        min_y -= margin
        max_y += margin

        width = max_x - min_x
        height = max_y - min_y

        return f"{min_x.evalf(100)} {min_y.evalf(100)} {width.evalf(100)} {height.evalf(100)}"

    def draw_unfolding(self, polygons):
        """
        @brief Generates an SVG file containing the given polygons.
        @param polygons List of lists of (x, y) tuples defining the polygons.
        """
        view_box = self.calculate_viewbox(polygons)

        with open(self.filename, 'w') as file:
            # Write the SVG header
            file.write(SVG_HEADER_TEMPLATE.format(view_box=view_box))

            # Write each polygon as a separate <polygon> element
            for polygon in polygons:
                points = " ".join(f"{x},{y.evalf(100)}" for x, y in polygon)
                file.write(f"""<polygon class="no_fill_black_stroke" points="{points}" />\n""")

            # Write the SVG footer
            file.write(SVG_FOOTER)


def main():
    """
    @brief Main function to execute the script.
    """
    # Output SVG file
    output_filename = "out.svg"

    # Get the number of sides from the user
    n = int(input("Enter the number of sides for the polygon (n >= 3): "))

    # Validate the input
    if n < 3:
        print("Error: A polygon must have at least 3 sides.")
        return

    # Define polygons
    polygons = [
        Polygon(n).calculate_polygon_vertices(cx=0, cy=0)
    ]

    # Draw polygons in an SVG file
    svg_drawer = SvgDrawer(output_filename)
    svg_drawer.draw_unfolding(polygons)


####################
if __name__ == "__main__":
    main()

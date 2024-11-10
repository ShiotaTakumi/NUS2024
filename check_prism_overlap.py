#! /usr/bin/env python

"""
@file prism_overlap.py
@brief Script to draw an SVG with multiple polygons
@detail This script generates an SVG file containing multiple polygons defined by given vertices.
@version 1.0.2
@date 2024-11-10
@author Takumi Shiota
"""

__version__ = "1.0.2"
__author__ = "Takumi Shiota"
__date__ = "2024-11-10"

# SVG header template with a placeholder for the dynamic viewBox
SVG_HEADER_TEMPLATE = """<?xml version="1.0" encoding="utf-8"?>
<svg xmlns="http://www.w3.org/2000/svg" version="1.1" viewBox="{view_box}">
<style>.no_fill_black_stroke {{fill:none; stroke:#000000; stroke-width:0.01;}}</style>
"""

# SVG footer
SVG_FOOTER = """</svg>
"""


class SvgDrawer:
    """
    @class SvgDrawer
    @brief A class to handle the creation and rendering of SVG files.
    """
    def __init__(self, filename):
        """
        @brief Constructor for SvgDrawer class
        @param filename Name of the output SVG file
        """
        self.filename = filename

    def calculate_viewbox(self, polygons, margin=0.2):
        """
        @brief Calculate a dynamic viewBox based on the given polygons
        @param polygons List of lists of (x, y) tuples defining the vertices of the polygons
        @param margin Additional space to include around the polygons in the viewBox
        @return A formatted string for the viewBox
        """
        # Flatten the list of polygons and find the min/max coordinates
        all_vertices = [vertex for polygon in polygons for vertex in polygon]
        min_x = min(x for x, _ in all_vertices)
        max_x = max(x for x, _ in all_vertices)
        min_y = min(y for _, y in all_vertices)
        max_y = max(y for _, y in all_vertices)

        # Add margin to each side
        min_x -= margin
        max_x += margin
        min_y -= margin
        max_y += margin

        # Calculate width and height of the viewBox
        width = max_x - min_x
        height = max_y - min_y

        return f"{min_x} {min_y} {width} {height}"

    def draw_unfolding(self, polygons):
        """
        @brief Create and save the SVG file containing multiple polygons
        @param polygons List of lists of (x, y) tuples defining the polygons
        """
        # Calculate the dynamic viewBox based on all polygons
        view_box = self.calculate_viewbox(polygons)

        # Write the SVG content to the file
        with open(self.filename, 'w') as file:
            # Write the header with the calculated viewBox
            file.write(SVG_HEADER_TEMPLATE.format(view_box=view_box))

            # Write each polygon as a separate <polygon> element
            for polygon in polygons:
                points = " ".join(f"{x},{y}" for x, y in polygon)
                file.write(f"""<polygon class="no_fill_black_stroke" points="{points}" />\n""")

            # Write the SVG footer
            file.write(SVG_FOOTER)


def main():
    """
    @brief Main function to execute the script
    """
    # Name of the output SVG file
    output_filename = "out.svg"

    # Define multiple polygons
    polygons = [
        [(0, 0), (0, 1), (1, 1), (1, 0)],   # First polygon
        [(1, 0), (1, 1), (1.8660, 0.5)],    # Second polygon
        [(0, 1), (1, 1), (1, 2), (0, 2)]    # Third polygon
    ]

    # Create an instance of SvgDrawer and render the SVG file
    svg_drawer = SvgDrawer(output_filename)
    svg_drawer.draw_unfolding(polygons)


####################
if __name__ == "__main__":
    main()

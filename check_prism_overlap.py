#! /usr/bin/env python

from sympy import pi, sin, cos, Float

"""
@file prism_overlap.py
@brief Script to generate an SVG file containing multiple polygons and rectangles.
@detail This script calculates the vertices of polygons (e.g., prism bases) and rectangles symbolically,
        and renders them as an SVG file. It also includes functionality to check for edge intersections
        between polygons.
@version 1.1.2
@date 2024-11-12
@author Takumi Shiota
"""

__version__ = "1.1.2"
__author__ = "Takumi Shiota"
__date__ = "2024-11-12"

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
    
    def calculate_b2_vertices(self, cx, cy):
        """
        @brief Calculates the vertices of a regular polygon (top base of the prism).
        @param cx X-coordinate of the center of the polygon.
        @param cy Y-coordinate of the center of the polygon.
        @return List of symbolic (x, y) tuples representing the vertices of the polygon.
        """
        angle_step = 2 * pi / self.n
        radius = self.side_length / (2 * sin(pi / self.n))

        cx += radius * cos(angle_step / 2)
        cy -= 0.5

        offset_angle = pi + angle_step / 2

        vertices = [
            (cx + radius * cos(offset_angle + i * angle_step),
             cy + radius * sin(offset_angle + i * angle_step))
            for i in range(self.n)
        ]
        return vertices
    
    def calculate_c1_vertices_v1(self, polygon_b2, c1_pos):
        """
        @brief Calculates the vertices of Rectangle C1 (version 1).
        @param polygon_b2 List of vertices of Polygon B2 (top base of the prism).
        @param c1_pos Position index for Rectangle C1.
        @return List of (x, y) tuples representing the vertices of the rectangle.
        """
        width = c1_pos + 1 - self.a_height
        height = self.height
        angle = ((2 * pi / self.n) * c1_pos) - pi

        cx = polygon_b2[c1_pos][0]
        cy = polygon_b2[c1_pos][1]

        vertices = [
            (cx, cy),
            (cx + width * sin(angle), cy - width * cos(angle)),
            (cx + width * sin(angle) + height * cos(angle), cy - width * cos(angle) + height * sin(angle)),
            (cx + height * cos(angle), cy + height * sin(angle))
        ]
        return vertices

    def calculate_c1_vertices_v2(self, polygon_b2, c1_pos):
        """
        @brief Calculates the vertices of Rectangle C1 (version 2).
        @param polygon_b2 List of vertices of Polygon B2 (top base of the prism).
        @param c1_pos Position index for Rectangle C1.
        @return List of (x, y) tuples representing the vertices of the rectangle.
        """
        width = self.n - c1_pos
        height = self.height
        angle = ((2 * pi / self.n) * (self.n - c1_pos)) - pi / 2

        cx = polygon_b2[c1_pos - 1][0]
        cy = polygon_b2[c1_pos - 1][1]

        vertices = [
            (cx, cy),
            (cx - width * cos(angle), cy + width * sin(angle)),
            (cx - width * cos(angle) + height * sin(angle), cy + width * sin(angle) + height * cos(angle)),
            (cx + height * sin(angle), cy + height * cos(angle))
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


class PolygonIntersection:
    """
    @class PolygonIntersection
    @brief Class for managing edge combinations and intersection detection between polygons.
    """

    def __init__(self, tolerance=1e-10):
        """
        @brief Initializes the PolygonIntersection class.
        @param tolerance Precision threshold for floating-point comparisons.
        """
        self.tolerance = tolerance

    def generate_edge_combinations(self, vertices1, vertices2):
        """
        @brief Generates all edge combinations between two polygons.
        @param vertices1 List of vertices (x, y) for the first polygon.
        @param vertices2 List of vertices (x, y) for the second polygon.
        @return List of tuples, where each tuple contains an edge from vertices1 and an edge from vertices2.
        """
        def edges_from_vertices(vertices):
            """
            @brief Generates edges (pairs of vertices) from a list of vertices.
            @param vertices List of vertices (x, y).
            @return List of edges, where each edge is a tuple of two vertices.
            """
            return [(vertices[i], vertices[(i + 1) % len(vertices)]) for i in range(len(vertices))]

        edges1 = edges_from_vertices(vertices1)
        edges2 = edges_from_vertices(vertices2)

        return [(edge1, edge2) for edge1 in edges1 for edge2 in edges2]
        
    def points_equal(self, point1, point2):
        """
        @brief Checks if two points are equal within a given tolerance.
        @param point1 Tuple (x1, y1) representing the first point.
        @param point2 Tuple (x2, y2) representing the second point.
        @return True if the points are equal, False otherwise.
        """
        x1, y1 = point1
        x2, y2 = point2
        return abs(x1 - x2) < self.tolerance and abs(y1 - y2) < self.tolerance
    
    def point_on_edge(self, point, edge):
        """
        @brief Checks if a point is on a given edge.
        @param point Tuple (px, py) representing the point.
        @param edge Tuple of two points ((x1, y1), (x2, y2)) representing the edge.
        @return True if the point lies on the edge, False otherwise.
        """
        px, py = point
        (x1, y1), (x2, y2) = edge

        # Check if the cross product is near zero (point lies on the line)
        cross_product = (px - x1) * (y2 - y1) - (py - y1) * (x2 - x1)
        if abs(cross_product) > self.tolerance:
            return False

        # Check if the point is within the edge bounds
        dot_product = (px - x1) * (x2 - x1) + (py - y1) * (y2 - y1)
        if dot_product < 0:
            return False

        squared_length = (x2 - x1) ** 2 + (y2 - y1) ** 2
        if dot_product > squared_length:
            return False

        return True
    

    def edge_intersection(self, edge1, edge2, evalf_precision=100):
        """
        @brief Checks if two edges intersect using cross products and checks for overlaps.
        @param edge1 A tuple of two points ((x1, y1), (x2, y2)) representing the first edge.
        @param edge2 A tuple of two points ((x3, y3), (x4, y4)) representing the second edge.
        @param evalf_precision Number of decimal digits for numerical evaluation.
        @return True if the edges intersect or overlap, False otherwise.
        """
        (x1, y1), (x2, y2) = edge1
        (x3, y3), (x4, y4) = edge2

        # Calculate determinants for the four cases
        d1 = ((x3 - x1) * (y2 - y1) - (y3 - y1) * (x2 - x1)).evalf(evalf_precision)
        d2 = ((x4 - x1) * (y2 - y1) - (y4 - y1) * (x2 - x1)).evalf(evalf_precision)
        d3 = ((x1 - x3) * (y4 - y3) - (y1 - y3) * (x4 - x3)).evalf(evalf_precision)
        d4 = ((x2 - x3) * (y4 - y3) - (y2 - y3) * (x4 - x3)).evalf(evalf_precision)

        # If the edges are collinear (all cross products are zero),
        # check if their bounding boxes overlap to determine if the edges overlap.
        if d1 == 0 and d2 == 0 and d3 == 0 and d4 == 0:
            # Calculate the bounds for both edges
            min_x1, max_x1 = sorted([x1, x2])
            min_y1, max_y1 = sorted([y1, y2])
            min_x2, max_x2 = sorted([x3, x4])
            min_y2, max_y2 = sorted([y3, y4])
            return not (max_x1 < min_x2 or max_x2 < min_x1 or max_y1 < min_y2 or max_y2 < min_y1)
        
        # Check for vertex-to-vertex contact
        vv_in_touch = (
            self.points_equal((x1, y1), (x3, y3)) or
            self.points_equal((x1, y1), (x4, y4)) or
            self.points_equal((x2, y2), (x3, y3)) or
            self.points_equal((x2, y2), (x4, y4))
        )
        
        # Check for vertex-to-edge contact
        ve_in_touch = (
            self.point_on_edge((x1, y1), edge2) or
            self.point_on_edge((x2, y2), edge2) or
            self.point_on_edge((x3, y3), edge1) or
            self.point_on_edge((x4, y4), edge1)
        )

        # Check for proper intersection
        proper_intersection = (d1 * d2 < 0 and d3 * d4 < 0)

        # Return True if edges intersect or touch
        return proper_intersection or vv_in_touch or ve_in_touch
    

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
    a_height = int(input(f"Enter the height of Rectangle A (n >= a_height): "))
    if a_height > n:
        print(f"Error: Height of Rectangle A must not exceed {n}.")
        return

    # Input the position of Rectangle C1
    c1_pos = int(input(f"Enter the position of Rectangle C1 (c1_pos > a_height): "))
    if c1_pos <= a_height:
        print(f"Error: Position of Rectangle C1 must be greater than {a_height}.")
        return

    # Create polygon and calculate vertices
    polygon = Polygon(n, height, a_height)
    b1_vertices = polygon.calculate_b1_vertices(cx=0, cy=0)  # Vertices of Polygon B1 (base of the prism)
    a_vertices = polygon.calculate_a_vertices(start=b1_vertices[n-1])  # Vertices of Rectangle A
    b2_vertices = polygon.calculate_b2_vertices(cx=a_vertices[2][0], cy=a_vertices[2][1])  # Vertices of Polygon B2 (top base of the prism)
    c1_vertices = polygon.calculate_c1_vertices_v1(polygon_b2=b2_vertices, c1_pos=c1_pos)  # Vertices of Rectangle C1 (I-1)
    # c1_vertices = polygon.calculate_c1_vertices_v2(polygon_b2=b2_vertices, c1_pos=c1_pos)  # Vertices of Rectangle C1 (I-2)

    # Test polygons for intersection testing (uncomment to use)
    # test_poly1 = [(Float(0), Float(0)), (Float(2), Float(0)), (Float(2), Float(2)), (Float(0), Float(2))]
    # test_poly2 = [(Float(1), Float(1)), (Float(3), Float(1)), (Float(3), Float(3)), (Float(1), Float(3))]
    # test_poly2 = [(Float(2), Float(0)), (Float(4), Float(0)), (Float(4), Float(2)), (Float(2), Float(2))]
    # test_poly2 = [(Float(2), Float(2)), (Float(4), Float(2)), (Float(4), Float(4)), (Float(2), Float(4))]
    # test_poly2 = [(Float(2), Float(1)), (Float(3), Float(0)), (Float(4), Float(1)), (Float(3), Float(2))]
    # test_poly2 = [(Float(3), Float(1)), (Float(5), Float(1)), (Float(5), Float(3)), (Float(3), Float(3))]

    # Render SVG with the calculated shapes
    svg_drawer = SvgDrawer(output_filename)
    svg_drawer.draw_unfolding([b1_vertices, a_vertices, b2_vertices, c1_vertices])
    # svg_drawer.draw_unfolding([test_poly1, test_poly2])  # Test edge combination generation

    # Initialize the PolygonIntersection class for edge combination generation
    intersec_checker = PolygonIntersection()

    # Generate all edge combinations between Polygon B1 and Rectangle C1
    edge_comb = intersec_checker.generate_edge_combinations(b1_vertices, c1_vertices)
    # edge_comb = intersec_checker.generate_edge_combinations(test_poly1, test_poly2)

    # Check if any edge combination intersects
    any_intersects = False

    for edge1, edge2 in edge_comb:
        intersects = intersec_checker.edge_intersection(edge1=edge1, edge2=edge2)
        if intersects:
            any_intersects = True
            break

    # Output the result with appropriate messages
    if any_intersects:
        print("Overlapping partial edge unfolding")
    else:
        print("Non-overlapping partial edge unfolding")


####################
if __name__ == "__main__":
    main()

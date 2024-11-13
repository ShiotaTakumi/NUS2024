#! /usr/bin/env python

from sympy import pi, sin, cos, Float

"""
@file prism_overlap.py
@brief Script to generate an SVG file containing multiple polygons and rectangles.
@detail This script calculates the vertices of polygons (e.g., prism bases) and rectangles symbolically,
        and renders them as an SVG file. It also includes functionality to check for edge intersections
        between polygons.
@version 1.1.3
@date 2024-11-13
@author Takumi Shiota
"""

__version__ = "1.1.3"
__author__ = "Takumi Shiota"
__date__ = "2024-11-13"

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

    def __init__(self, n, h, a):
        """
        @brief Initializes prism's edge unfolding parameter.
        @param n The number of sides for the polygon.
        @param h The height of the prism.
        @param a The length of one side of Rectangle A.
        """
        self.n = n
        self.h = h
        self.a = a

    def calculate_b1_vertices(self):
        """
        @brief Calculates the vertices of polygon B1.
        @return List of symbolic (x, y) tuples representing vertices of B1.
        """
        angle_step = 2 * pi / self.n
        radius = 1 / (2 * sin(pi / self.n))
        offset_angle = angle_step / 2

        b1_0_x = Float(0)
        b1_0_y = Float(0)

        vertices = [
            (b1_0_x + radius * cos(offset_angle + i * angle_step),
             b1_0_y + radius * sin(offset_angle + i * angle_step))
            for i in range(self.n)
        ]
        return vertices

    def calculate_a_vertices(self, b1_vertices):
        """
        @brief Calculates the vertices of rectangle A.
        @param b1_vertices List of vertices of polygon B1.
        @return List of symbolic (x, y) tuples representing vertices of rectangle A.
        """
        a0_x = b1_vertices[self.n-1][0]
        a0_y = b1_vertices[self.n-1][1]

        vertices = [
            (a0_x, a0_y),  # a0
            (a0_x + self.h, a0_y),  # a1
            (a0_x + self.h, a0_y + self.a),  # a2
            (a0_x, a0_y + self.a)  # a3
        ]
        return vertices
    
    def calculate_b2_vertices(self, a_vertices):
        """
        @brief Calculates the vertices of polygon B2.
        @param a_vertices List of vertices of rectangle A.
        @return List of symbolic (x, y) tuples representing vertices of polygon B2.
        """
        angle_step = 2 * pi / self.n
        radius = 1 / (2 * sin(pi / self.n))

        b2_0_x = a_vertices[2][0] + radius * cos(angle_step / 2)
        b2_0_y = a_vertices[2][1] - radius * sin(angle_step / 2)

        offset_angle = pi + angle_step / 2

        vertices = [
            (b2_0_x + radius * cos(offset_angle + i * angle_step),
             b2_0_y + radius * sin(offset_angle + i * angle_step))
            for i in range(self.n)
        ]
        return vertices
    
    def calculate_c2_vertices_I_1(self, polygon_b2, c2):
        """
        @brief Calculates the vertices of polygon C2.
        @param polygon_b2 List of vertices of Polygon B2.
        @param c2 Position of Rectangle C2.
        @return List of symbolic (x, y) tuples representing vertices of rectangle C2.
        """
        strip = (c2 - self.a) + 1
        th = ((2 * pi / self.n) * c2)

        c2_0_x = polygon_b2[c2][0]
        c2_0_y = polygon_b2[c2][1]

        vertices = [
            (c2_0_x, c2_0_y),  # c2_0
            (c2_0_x - strip * cos(th), c2_0_y + strip * sin(th)),  # c2_1
            (c2_0_x - strip * cos(th) - self.h * sin(th), c2_0_y + strip * sin(th) - self.h * cos(th)),  # c2_2
            (c2_0_x - self.h * sin(th), c2_0_y - self.h * cos(th))  # c2_3
        ]
        return vertices
    
    def calculate_c2_vertices_I_2(self, polygon_b2, c2):
        """
        @brief Calculates the vertices of polygon C2.
        @param polygon_b2 List of vertices of Polygon B2.
        @param c2 Position of Rectangle C2.
        @return List of symbolic (x, y) tuples representing vertices of rectangle C2.
        """
        strip = ((self.n - 1) - c2) + 1
        th = (2 * pi / self.n) * (self.n - c2)

        c2_0_x = polygon_b2[c2 - 1][0]
        c2_0_y = polygon_b2[c2 - 1][1]

        vertices = [
            (c2_0_x, c2_0_y),  # c2_0
            (c2_0_x - strip * sin(th), c2_0_y - strip * cos(th)),  # c2_1
            (c2_0_x - strip * sin(th) - self.h * cos(th), c2_0_y - strip * cos(th) + self.h * sin(th)),  # c2_2
            (c2_0_x - self.h * cos(th), c2_0_y + self.h * sin(th))  # c2_3
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

    # Ask the user for the type of overlap to check
    type = input("Enter the type of overlap to check (I-1, I-2): ").strip()

    # Input number of sides for the polygon
    n = int(input(f"Enter number of sides for the polygon n (n >= 3): "))
    if n < 3:
        print(f"Error: n < 3.")
        return

    # Input the height of the prism
    h = float(input(f"Enter the height of the prism h (h > 0): "))
    if h <= 0:
        print(f"Error: h <= 0.")
        return

    # Input the length of one side of Rectangle A
    a = int(input(f"Enter the length of one side of Rectangle A ({n} > a): "))
    if a > n:
        print(f"Error: a <= {n}.")
        return

    # Input the position of Rectangle C2
    c2 = int(input(f"Enter the position of Rectangle C2 ({a} < c2 < {n}): "))
    if c2 <= a or c2 >= n:
        print(f"Error: c2 <= {a} or c2 >= {n}.")
        return

    # Create each polygon's vertices
    polygon = Polygon(n, h, a)
    b1_vertices = polygon.calculate_b1_vertices()  # Vertices of polygon B1
    a_vertices = polygon.calculate_a_vertices(b1_vertices)  # Vertices of rectangle A

    if type == 'I-1':
        b2_vertices = polygon.calculate_b2_vertices(a_vertices)  # Vertices of polygon B2
        c2_vertices = polygon.calculate_c2_vertices_I_1(b2_vertices, c2)  # Vertices of rectangle C2
    elif type == 'I-2':
        b2_vertices = polygon.calculate_b2_vertices(a_vertices)  # Vertices of polygon B2
        c2_vertices = polygon.calculate_c2_vertices_I_2(b2_vertices, c2)  # Vertices of rectangle C2
    else:
        # Polygons for intersection testing (uncomment to use)
        test_poly1 = [(Float(0), Float(0)), (Float(2), Float(0)), (Float(2), Float(2)), (Float(0), Float(2))]
        test_poly2 = [(Float(1), Float(1)), (Float(3), Float(1)), (Float(3), Float(3)), (Float(1), Float(3))]
        # test_poly2 = [(Float(2), Float(0)), (Float(4), Float(0)), (Float(4), Float(2)), (Float(2), Float(2))]
        # test_poly2 = [(Float(2), Float(2)), (Float(4), Float(2)), (Float(4), Float(4)), (Float(2), Float(4))]
        # test_poly2 = [(Float(2), Float(1)), (Float(3), Float(0)), (Float(4), Float(1)), (Float(3), Float(2))]
        # test_poly2 = [(Float(3), Float(1)), (Float(5), Float(1)), (Float(5), Float(3)), (Float(3), Float(3))]

    # Render SVG with the calculated shapes
    svg_drawer = SvgDrawer(output_filename)
    if type == 'I-1' or type == 'I-2':
        svg_drawer.draw_unfolding([b1_vertices, a_vertices, b2_vertices, c2_vertices])
    else:
        svg_drawer.draw_unfolding([test_poly1, test_poly2])  # Test edge combination generation

    # Initialize the PolygonIntersection class for edge combination generation
    intersec_checker = PolygonIntersection()

    # Generate all edge combinations between polygon X and polygon Y
    if type == 'I-1' or type == 'I-2':
        edge_comb = intersec_checker.generate_edge_combinations(b1_vertices, c2_vertices)
    else: 
        edge_comb = intersec_checker.generate_edge_combinations(test_poly1, test_poly2)

    # Check if any edge combination intersects
    any_intersects = False

    for edge1, edge2 in edge_comb:
        intersects = intersec_checker.edge_intersection(edge1=edge1, edge2=edge2)
        if intersects:
            any_intersects = True
            break

    # Output the result with appropriate messages
    if any_intersects:
        print("Overlapping")
    else:
        print("Not overlapping")


####################
if __name__ == "__main__":
    main()

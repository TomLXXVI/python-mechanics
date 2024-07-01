from mechanics import Quantity
from mechanics.section.geometry import Shape, Polygon, HollowPolygon
from mechanics.utils import LineSegment


Q_ = Quantity


class TransverseShear:

    def __init__(
        self,
        V: Quantity,
        shape: Shape
    ) -> None:
        self.V = V
        self.shape = shape
        self.tau_max = self.tau(y=Q_(0.0, 'mm'))  # tau at centroid of cross-section

    def shear_flow(self, y: Quantity) -> Quantity:
        """Returns the shear force per unit length of beam at height `y`
        relative to the centroid of the cross-section.
        """
        if isinstance(self.shape, Polygon):
            pg = TransverseShear._cut_polygon(y, self.shape)
            A_acc = pg.area
            _, y_acc = pg.centroid
            # returns the coordinates of the centroid of the cut polygon with
            # respect to the centroid of the original shape.
            Q = A_acc * abs(y_acc)
            I_xx = self.shape.moment_of_inertia_xx
            q = self.V * Q / I_xx
            return q

    def tau(self, y: Quantity) -> Quantity:
        """Returns the shear stress at height `y` relative to the centroid of
        the cross-section.
        """
        if isinstance(self.shape, (Polygon, HollowPolygon)):
            pg = TransverseShear._cut_polygon(y, self.shape)
            A_acc = pg.area
            _, y_acc = pg.centroid
            Q = A_acc * abs(y_acc)
            I_xx = self.shape.moment_of_inertia_xx
            q = self.V * Q / I_xx
            t = TransverseShear.__determine_t(y, pg)
            if t is None:
                y = y - Q_(0.01, y.units)
                pg = TransverseShear._cut_polygon(y, self.shape)
                t = TransverseShear.__determine_t(y, pg)
            tau = q / t
            return tau

    @staticmethod
    def __create_polygon_line_segments(polygon: Polygon) -> list[LineSegment]:
        i_max = len(polygon.vertices) - 1
        num = len(polygon.vertices)
        segments = []
        for i in range(num):
            c = i
            n = i + 1 if i + 1 <= i_max else 0
            vertex_1 = (
                round(float(polygon.vertices[c][0].m), 12),
                round(float(polygon.vertices[c][1].m), 12)
            )
            vertex_2 = (
                round(float(polygon.vertices[n][0].m), 12),
                round(float(polygon.vertices[n][1].m), 12)
            )
            segments.append(LineSegment(vertex_1, vertex_2))
        return segments

    @staticmethod
    def __create_line_segments(
        polygon: Polygon | HollowPolygon
    ) -> list[LineSegment]:
        if isinstance(polygon, Polygon):
            segments = TransverseShear.__create_polygon_line_segments(polygon)
            return segments
        elif isinstance(polygon, HollowPolygon):
            outer_segments = TransverseShear.__create_polygon_line_segments(polygon.outer_polygon)
            inner_segments = TransverseShear.__create_polygon_line_segments(polygon.inner_polygon)
            inner_segments = [s.reverse() for s in reversed(inner_segments)]
            segments = outer_segments + inner_segments
            return segments

    @staticmethod
    def __get_intersections(
        y: float,
        segments: list[LineSegment]
    ) -> list[tuple[float, float]]:
        x_coords = [
            x
            for segment in segments
            if (x := segment.x(y)) is not None
        ]
        intersections = [(x, y) for x in x_coords]
        return intersections

    @staticmethod
    def __order_segments(y: float, segments: list[LineSegment]):
        # Remove line segments that lay above (or below) the cut line:
        valid_segments = []
        for segment in segments:
            if segment.is_horizontal():
                if y >= 0 and segment.y1 < y:
                    continue
                elif y < 0 and segment.y1 > y:
                    continue
                else:
                    valid_segments.append(segment)
            else:
                valid_segments.append(segment)
        # Order line segments in consecutive order:
        current_segment = valid_segments.pop(0)
        ordered_segments = [current_segment]
        while valid_segments:
            flags = [False, False]
            # Look for a line segment of which the start point coincides with
            # the end point of the current segment:
            for segment in valid_segments:
                if segment.p1 == current_segment.p2:
                    i = valid_segments.index(segment)
                    ordered_segments.append(valid_segments.pop(i))
                    flags[1] = True
                    break
            # Look for a line segment of which the end point coincides with the
            # start point of the current segment:
            for segment in valid_segments:
                if segment.p2 == current_segment.p1:
                    i = valid_segments.index(segment)
                    j = ordered_segments.index(current_segment)
                    ordered_segments.insert(j, valid_segments.pop(i))
                    flags[0] = True
                    break

            if flags[1] is True:
                current_segment = ordered_segments[-1]
            elif flags[0] is True:
                j = ordered_segments.index(current_segment)
                current_segment = ordered_segments[j-1]
            else:
                current_segment = valid_segments.pop(0)
                ordered_segments.append(None)
                ordered_segments.append(current_segment)
        ordered_segments = [s for s in ordered_segments if s is not None]
        return ordered_segments

    @staticmethod
    def __get_vertices(
        y: float,
        segments: list[LineSegment]
    ) -> list[tuple[float, float]]:
        vertices = []
        if y >= 0:
            for segment in segments:
                if segment.y1 >= y:
                    vertices.append(segment.p1)
                if segment.y2 >= y:
                    vertices.append(segment.p2)
        else:
            for segment in segments:
                if segment.y1 <= y:
                    vertices.append(segment.p1)
                if segment.y2 <= y:
                    vertices.append(segment.p2)
        return vertices

    @staticmethod
    def __remove_duplicates(
        vertices: list[tuple[float, float]]
    ) -> list[tuple[float, float]]:
        unique_vertices = [vertices[0]]
        for vertex in vertices:
            flags = []
            for unique_vertex in unique_vertices:
                cond1 = vertex[0] == unique_vertex[0]
                cond2 = vertex[1] == unique_vertex[1]
                if not (cond1 and cond2):
                    flags.append(True)
                else:
                    flags.append(False)
            if all(flags):
                unique_vertices.append(vertex)
        return unique_vertices

    @staticmethod
    def __convert_to_quantity(
        vertices: list[tuple[float, float]],
        units: str
    ) -> list[Quantity]:
        vertices = [Q_(vertex, units) for vertex in vertices]
        return vertices

    @staticmethod
    def _cut_polygon(y: Quantity, shape: Polygon | HollowPolygon) -> Polygon:
        y_mag = float(y.magnitude)
        segments = TransverseShear.__create_line_segments(shape)
        intersections = TransverseShear.__get_intersections(y_mag, segments)
        for i in range(len(segments)):
            segment = segments[i]
            for intersection in intersections:
                if segment.contains((intersection[0], intersection[1])):
                    if segment.slope < 0:
                        if y_mag >= 0:
                            new_segment = LineSegment(
                                p1=segment.p1,
                                p2=(intersection[0], intersection[1])
                            )
                        else:
                            new_segment = LineSegment(
                                p1=(intersection[0], intersection[1]),
                                p2=segment.p2
                            )
                        segments[i] = new_segment
                    if segment.slope > 0:
                        if y_mag >= 0:
                            new_segment = LineSegment(
                                p1=(intersection[0], intersection[1]),
                                p2=segment.p2
                            )
                        else:
                            new_segment = LineSegment(
                                p1=segment.p1,
                                p2=(intersection[0], intersection[1])
                            )
                        segments[i] = new_segment
                    if segment.slope == 0:
                        pass
        segments = TransverseShear.__order_segments(y_mag, segments)
        vertices = TransverseShear.__get_vertices(y_mag, segments)
        vertices = TransverseShear.__remove_duplicates(vertices)
        vertices = TransverseShear.__convert_to_quantity(vertices, y.units)
        polygon = Polygon(vertices)
        return polygon

    @staticmethod
    def __determine_t(y: Quantity, cut_polygon: Polygon) -> Quantity | None:
        # Get horizontal line segments from cut polygon:
        horizontal_segments = [
            segment
            for segment in TransverseShear.__create_line_segments(cut_polygon)
            if segment.is_horizontal()
        ]
        # Determine t:
        if y.m >= 0:
            # cut above or through centroid of original shape
            min_y1 = min(s.y1 for s in horizontal_segments)
            segments_min_y1 = [s for s in horizontal_segments if s.y1 == min_y1]
            sum_of_lengths, total_length = TransverseShear.__get_sum_of_lengths(segments_min_y1)
            if len(segments_min_y1) == 1:
                t = segments_min_y1[0].length
                return Q_(t, y.units)
            elif (len(segments_min_y1) > 1) and (sum_of_lengths == total_length):
                # multiple connected horizontal line segments at min_y1
                return None
            else:
                # multiple disconnected horizontal line segments at min_y1
                t = sum_of_lengths
                return Q_(t, y.units)
        if y.m < 0:
            # cut below centroid of original shape
            max_y1 = max(s.y1 for s in horizontal_segments)
            segments_max_y1 = [s for s in horizontal_segments if s.y1 == max_y1]
            if len(segments_max_y1) == 1:
                t = segments_max_y1[0].length
                return Q_(t, y.units)
            elif len(segments_max_y1) > 1:
                if TransverseShear.__segments_do_overlap(segments_max_y1):
                    # multiple overlapping horizontal line segments at max_y1
                    # E.g. in case of a C-shape: the cut will result in two
                    # horizontal segments that overlap.
                    segments_max_y1.sort(key=lambda s: s.length)
                    l_max = segments_max_y1[-1].length
                    l_delta = sum(s.length for s in segments_max_y1[:-1])
                    t = l_max - l_delta
                else:
                    # multiple disconnected horizontal line segments at max_y1
                    sum_of_lengths, _ = TransverseShear.__get_sum_of_lengths(segments_max_y1)
                    t = sum_of_lengths
                return Q_(t, y.units)
            else:
                raise ValueError("`t` in shear formula cannot be determined")

    @staticmethod
    def __get_sum_of_lengths(
        segments: list[LineSegment]
    ) -> tuple[float, float]:
        segments.sort(key=lambda s: s.x1)
        first_x = segments[0].x1
        last_x = segments[-1].x2
        total_length = last_x - first_x
        sum_of_lengths = sum(s.length for s in segments)
        return sum_of_lengths, total_length

    @staticmethod
    def __segments_do_overlap(segments: list[LineSegment]) -> bool:
        i_max = len(segments) - 1
        for i in range(len(segments)):
            curr_segment = segments[i]
            other_segments = []
            if 0 < i < i_max:
                other_segments = segments[:i] + segments[i+1:]
            elif i == 0:
                other_segments = segments[1:]
            elif i == i_max:
                other_segments = segments[:-1]
            for segment in other_segments:
                if curr_segment.coincide(segment):
                    return True
        return False

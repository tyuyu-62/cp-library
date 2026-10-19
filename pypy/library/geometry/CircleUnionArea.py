# Library: https://github.com/tyuyu-62/cp-library
from math import acos, atan2, cos, pi, sin


def circle_union_area(circles) -> float:
    """Return the area covered by the union of circles given as (x, y, radius) triples.
    Nonpositive-radius circles contribute no area, and an empty input returns 0.0.

    Constraints:
        Circle centers and radii are finite numeric values.

    Time:
        O(N^2 log N), where N is the number of circles.

    Space:
        O(N) auxiliary space.
    """
    circles = [(x, y, radius) for x, y, radius in circles if radius > 0]
    if not circles:
        return 0.0
    parent = list(range(len(circles)))
    sizes = [1] * len(circles)
    arcs = []

    def leader(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    result = 0.0
    full_turn = 2.0 * pi
    for index, (x, y, radius) in enumerate(circles):
        if radius <= 0:
            continue
        events = []
        initial = 0
        covered = False
        for other_index, (other_x, other_y, other_radius) in enumerate(circles):
            if index == other_index or other_radius <= 0:
                continue
            dx = other_x - x
            dy = other_y - y
            distance2 = dx * dx + dy * dy
            if distance2 < (radius + other_radius) ** 2:
                first = leader(index)
                second = leader(other_index)
                if first != second:
                    if sizes[first] < sizes[second]:
                        first, second = second, first
                    parent[second] = first
                    sizes[first] += sizes[second]
            if distance2 == 0:
                if other_radius > radius or (
                    other_radius == radius and other_index < index
                ):
                    covered = True
                    break
                continue
            distance = distance2 ** 0.5
            if distance + radius <= other_radius:
                covered = True
                break
            if distance + other_radius <= radius or distance >= radius + other_radius:
                continue
            angle = atan2(dy, dx)
            cosine = (
                (radius * radius + distance2 - other_radius * other_radius)
                / (2.0 * radius * distance)
            )
            delta = acos(max(-1.0, min(1.0, cosine)))
            start = (angle - delta) % full_turn
            end = (angle + delta) % full_turn
            if start <= end:
                events.append((start, 1))
                events.append((end, -1))
            else:
                initial += 1
                events.append((start, 1))
                events.append((end, -1))
        if covered:
            continue
        if not events:
            result += pi * radius * radius
            continue
        events.sort()
        arc_width = 0.0
        arc_sine = 0.0
        arc_cosine = 0.0
        previous = 0.0
        active = initial
        event_index = 0
        while event_index < len(events):
            angle = events[event_index][0]
            if active == 0:
                width = angle - previous
                arc_width += width
                arc_sine += sin(angle) - sin(previous)
                arc_cosine += cos(previous) - cos(angle)
            while event_index < len(events) and events[event_index][0] == angle:
                active += events[event_index][1]
                event_index += 1
            previous = angle
        if active == 0:
            width = full_turn - previous
            arc_width += width
            arc_sine -= sin(previous)
            arc_cosine += cos(previous) - 1.0
        arcs.append((index, arc_width, arc_sine, arc_cosine))
    for index, width, sine, cosine in arcs:
        x, y, radius = circles[index]
        origin_x, origin_y = circles[leader(index)][:2]
        result += 0.5 * (
            radius * radius * width
            + radius * (x - origin_x) * sine
            + radius * (y - origin_y) * cosine
        )
    return result

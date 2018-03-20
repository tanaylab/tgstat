/*
 * Point.h
 *
 *  Created on: Jan 3, 2012
 *      Author: hoichman
 */

#ifndef POINT_H_
#define POINT_H_

#include <stdio.h>
#include <sys/types.h>

#include <algorithm>

#include "Rectangle.h"

using namespace std;

#pragma pack(push)
#pragma pack(8)

struct Point {
	int64_t x;
	int64_t y;

	Point() {}

	Point(const Point &p) : x(p.x), y(p.y) {}

	Point(int64_t _x, int64_t _y) : x(_x), y(_y) {}

	int64_t get_x1() const { return x; }
	int64_t get_x2() const { return x + 1; }
	int64_t get_y1() const { return y; }
	int64_t get_y2() const { return y + 1; }

	bool do_intersect(const Rectangle &rect) const {
		return x >= rect.x1 && x < rect.x2 && y >= rect.y1 && y < rect.y2;
	}

	bool is_inside(const Rectangle &rect) const { return do_intersect(rect); }

	const Rectangle intersect(const Rectangle &rect) const {
		return Rectangle(max(x, rect.x1), max(y, rect.y1), min(x + 1, rect.x2), min(y + 1, rect.y2));
	}

	Rectangle intersect(const Rectangle &r1, const Rectangle &r2) const {
		return Rectangle(max(x, max(r1.x1, r2.x1)), max(y, max(r1.y1, r2.y1)), min(x + 1, min(r1.x2, r2.x2)), min(y + 1, min(r1.y2, r2.y2)));
	}

	int64_t intersected_area(const Rectangle &rect) const {
		return do_intersect(rect) ? 1L : 0L;
	}

	int64_t intersected_area(const Rectangle &r1, const Rectangle &r2) const {
		return x >= max(r1.x1, r2.x1) && x < min(r1.x2, r2.x2) && y >= max(r1.y1, r2.y1) && y < min(r1.y2, r2.y2) ? 1L : 0L;
	}

	int64_t xdist(const Rectangle &rect, bool touch_is_at_dist_one = false) const {
		if (x >= rect.x2)
			return x - rect.x2 + touch_is_at_dist_one;
		if (x < rect.x1)
			return rect.x1 - x - touch_is_at_dist_one;
		return 0;
	}

	int64_t ydist(const Rectangle &rect, bool touch_is_at_dist_one = false) const {
		if (y >= rect.y2)
			return y - rect.y2 + touch_is_at_dist_one;
		if (y < rect.y1)
			return rect.y1 - y - touch_is_at_dist_one;
		return 0;
	}

	int64_t manhattan_dist(const Rectangle &rect, bool touch_is_at_dist_one = false) const { return xdist(rect, touch_is_at_dist_one) + ydist(rect, touch_is_at_dist_one); }

	char *debug_str() const {
		static char str[200];
		sprintf(str, "(%lld, %lld)", x, y);
		return str;
	}
};

#pragma pack(pop)

#endif /* POINT_H_ */

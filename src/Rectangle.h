/*
 * Rectangle.h
 *
 *  Created on: Jan 3, 2012
 *      Author: hoichman
 */

#ifndef RECTANGLE_H_
#define RECTANGLE_H_

#include <stdio.h>
#include <sys/types.h>

#include <algorithm>
#include <vector>

using namespace std;

#pragma pack(push)
#pragma pack(8)

struct Rectangle {
	int64_t x1;
	int64_t y1;
	int64_t x2;
	int64_t y2;

	Rectangle() {}

	Rectangle(const Rectangle &r) : x1(r.x1), y1(r.y1), x2(r.x2), y2(r.y2) {}

	Rectangle(int64_t _x1, int64_t _y1, int64_t _x2, int64_t _y2) : x1(_x1), y1(_y1), x2(_x2), y2(_y2) {}

	int64_t get_x1() const { return x1; }
	int64_t get_x2() const { return x2; }
	int64_t get_y1() const { return y1; }
	int64_t get_y2() const { return y2; }

	bool operator==(const Rectangle &r) const { return x1 == r.x1 && x2 == r.x2 && y1 == r.y1 && y2 == r.y2; }

	bool is_point() const { return x1 == x2 - 1 && y1 == y2 - 1; }

	bool do_intersect(const Rectangle &rect) const {
		return max(x1, rect.x1) < min(x2, rect.x2) && max(y1, rect.y1) < min(y2, rect.y2);
	}

	bool is_inside(const Rectangle &rect) const {
		return x1 >= rect.x1 && y1 >= rect.y1 && x2 <= rect.x2 && y2 <= rect.y2;
	}

	Rectangle intersect(const Rectangle &rect) const {
		return Rectangle(max(x1, rect.x1), max(y1, rect.y1), min(x2, rect.x2), min(y2, rect.y2));
	}

	Rectangle intersect(const Rectangle &r1, const Rectangle &r2) const {
		return Rectangle(max(x1, max(r1.x1, r2.x1)), max(y1, max(r1.y1, r2.y1)), min(x2, min(r1.x2, r2.x2)), min(y2, min(r1.y2, r2.y2)));
	}

	int64_t intersected_area(const Rectangle &rect) const {
		int64_t width = min(x2, rect.x2) - max(x1, rect.x1);
		int64_t height = min(y2, rect.y2) - max(y1, rect.y1);
		return width > 0 && height > 0 ? width * height : 0;
	}

	int64_t intersected_area(const Rectangle &r1, const Rectangle &r2) const {
		int64_t width = min(x2, min(r1.x2, r2.x2)) - max(x1, max(r1.x1, r2.x1));
		int64_t height = min(y2, min(r1.y2, r2.y2)) - max(y1, max(r1.y1, r2.y1));
		return width > 0 && height > 0 ? width * height : 0;
	}

	int64_t xdist(const Rectangle &rect, bool touch_is_at_dist_one = false) const {
		if (x1 >= rect.x2)
			return x1 - rect.x2 + touch_is_at_dist_one;
		if (x2 <= rect.x1)
			return rect.x1 - x2 + touch_is_at_dist_one;
		return 0;
	}

	int64_t ydist(const Rectangle &rect, bool touch_is_at_dist_one = false) const {
		if (y1 >= rect.y2)
			return y1 - rect.y2 + touch_is_at_dist_one;
		if (y2 <= rect.y1)
			return rect.y1 - y2 + touch_is_at_dist_one;
		return 0;
	}

	int64_t manhattan_dist(const Rectangle &rect, bool touch_is_at_dist_one = false) const { return xdist(rect, touch_is_at_dist_one) + ydist(rect, touch_is_at_dist_one); }

	int64_t width() const { return x2 - x1; }
	int64_t height() const { return y2 - y1; }
	int64_t area() const { return width() * height(); }
	bool    is_non_empty_area() const { return x2 > x1 && y2 > y1; }

	char *debug_str() const {
		static char str[200];
		sprintf(str, "(%ld - %ld) (%ld - %ld)", x1, x2, y1, y2);
		return str;
	}
};

#pragma pack(pop)

typedef vector<Rectangle> Rectangles;

#endif /* RECTANGLE_H_ */

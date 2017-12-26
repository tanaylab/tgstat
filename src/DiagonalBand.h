#ifndef DIAGONALBAND_H_
#define DIAGONALBAND_H_

#include "Rectangle.h"

// Defines a band bounded by two diagonals at 45 degrees angle each: [d1, d2)

struct DiagonalBand {
	int64_t d1;
	int64_t d2;

	DiagonalBand() : d1(0), d2(0) {}

	// Creates a 45 degrees band by giving two distances from a 45 degree diagonal.
	// Diagonal is defined by the delta of its x and y coordinates (x-y).
	// d1 must be lesser than d2.
	DiagonalBand(int64_t _d1, int64_t _d2) : d1(_d1), d2(_d2) {}

	bool do_intersect(const Rectangle &rect) const;

	bool do_contain(const Rectangle &rect) const;

	// shrinks the rectangle to the minimal area such that it still contains the area of its original intersection with the band
	void shrink2intersected(Rectangle &rect) const;

	int64_t intersected_area(const Rectangle &shrinked_rect) const;

	bool is_non_empty_area() const { return d1 < d2; }

	int64_t x1(int64_t y1) const { return y1 + d1; }
	int64_t y1(int64_t x1) const { return x1 - d1; }
	int64_t x2(int64_t y2) const { return y2 + d2; }
	int64_t y2(int64_t x2) const { return x2 - d2; }
};


//--------------------------------------- IMPLEMENTATION ------------------------------------------------

inline bool DiagonalBand::do_intersect(const Rectangle &rect) const
{
	return rect.x2 - rect.y1 > d1 && rect.x1 - rect.y2 + 1 < d2;
}

inline bool DiagonalBand::do_contain(const Rectangle &rect) const
{
	return rect.x1 - rect.y2 + 1 >= d1 && rect.x2 - rect.y1 <= d2;
}

inline void DiagonalBand::shrink2intersected(Rectangle &rect) const
{
	if (rect.x1 - rect.y1 < d1)
		rect.x1 = x1(rect.y1);
	else if (rect.x1 - rect.y1 > d2)
		rect.y1 = y2(rect.x1);

	if (rect.x2 - rect.y2 < d1)
		rect.y2 = y1(rect.x2);
	else if (rect.x2 - rect.y2 > d2)
		rect.x2 = x2(rect.y2);
}

inline int64_t DiagonalBand::intersected_area(const Rectangle &shrinked_rect) const
{
	int64_t area = shrinked_rect.area();

	// substract the area of the triangle above d1
	if (shrinked_rect.x1 - shrinked_rect.y2 + 1 < d1) {
		int64_t n = x1(shrinked_rect.y2) - shrinked_rect.x1;
		area -= (n * n - n) >> 1;
	}

	// substract the area of the triangle below d2
	if (shrinked_rect.x2 - shrinked_rect.y1 > d2) {
		int64_t n = shrinked_rect.x2 - x2(shrinked_rect.y1);
		area -= (n * n + n) >> 1;
	}

	return area;
}

#endif /* DIAGONALBAND_H_ */

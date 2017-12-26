#ifndef SEGMENT_H_INCLUDED
#define SEGMENT_H_INCLUDED

#include <math.h>
#include <stdlib.h>

#pragma pack(push)
#pragma pack(8)

struct Segment {
	int64_t start;
	int64_t end;

	Segment() {}

	Segment(const Segment &s) : start(s.start), end(s.end) {}

	Segment(int64_t _start, int64_t _end) : start(_start), end(_end) {}

	bool operator==(const Segment &s) const { return start == s.start && end == s.end; }

	int64_t range() const { return end - start; }

	// returns distance from point coordinate to this segment. If touch_is_at_dist_one is true the distance between coordinate == end is 1 otherwise 0.
	double dist2coord(int64_t coord, bool touch_is_at_dist_one = false) const { return coord >= start && coord < end ? 0 : min(coord - start, end - coord + touch_is_at_dist_one); }

	// returns distance from this segment to s. If touch_is_at_dist_one is true the distance between touching segments is 1 otherwise 0.
	int64_t dist2segment(const Segment &s, bool touch_is_at_dist_one = false) const {
		return max(start, s.start) < min(end, s.end) ? 0 : min(llabs(s.start - end + touch_is_at_dist_one), llabs(s.end - start + touch_is_at_dist_one));
	}

	// returns distance from point coordinate to the center of segment; returns NaN if the coordinate is outside of segment
	double dist2center(int64_t coord) const { coord < start || coord >= end ? numeric_limits<double>::quiet_NaN() : fabs(coord - (start + end) / 2); }

	// returns true if the two segments touch or overlap each other
	bool do_touch(const Segment &s) const { return max(start, s.start) <= min(end, s.end); }

	// returns true if the two segments overlap each other
	bool do_overlap(const Segment &s) const { return max(start, s.start) < min(end, s.end); }

	// returns true if the segment fully overlaps the given segment
	bool do_contain(const Segment &s) const { return start <= s.start && end >= s.end; }

	char *debug_str() const {
		static char str[200];
		sprintf(str, "(%ld - %ld)", start, end);
		return str;
	}
};

#pragma pack(pop)

typedef vector<Segment> Segments;

#endif


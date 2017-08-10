#ifndef GINTERVAL_H_
#define GINTERVAL_H_

#include <math.h>
#include <sys/types.h>
#include <iostream>
#include <limits>
#include <vector>

#include "GenomeChromKey.h"
#include "Segment.h"
#include "TGLException.h"

//------------------------------------- GInterval -----------------------------------------------
// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

struct GInterval : public Segment {
	int     chromid;
	char    strand;
	void   *udata;

	enum Errors { BAD_INTERVAL, BAD_STRAND };
	enum { CHROM, START, END, NUM_COLS };

	static const char *COL_NAMES[];

	GInterval() : Segment(-1, -1), chromid(-1), strand(0), udata(NULL) {}
	GInterval(int _chromid, int64_t _start, int64_t _end, char _strand, void *_udata = NULL) : Segment(_start, _end), chromid(_chromid), strand(_strand), udata(_udata) {}

	// compares two intervals by chrom1 and chrom2 and then by start coordinate
	bool operator<(const GInterval &interv) const;

	template <class T> static void *cast2udata(T v);

	int64_t range() const { return end - start; }

	// verifies basic interval correctness
	void verify(const GenomeChromKey &chromkey, bool check_chrom_boundary = true) const;

	// returns distance from point coordinate to this interval; centerbins is the value of margin
	double dist2coord(int64_t coord, double margin) const;

	// returns distance from this interval to interv. If touch_is_at_dist_one is true the distance between touching intervals is 1 otherwise 0.
	int64_t dist2interv(const GInterval &interv, bool touch_is_at_dist_one = false) const;

	// returns distance from point coordinate to the center of interval; returns NaN if the coordinate is outside of interval
	double dist2center(int64_t coord) const;

	// returns true if the two intervals touch or overlap each other
	bool do_touch(const GInterval &interv) const { return chromid == interv.chromid && Segment::do_touch(interv); }

	// returns true if the two intervals overlap each other
	bool do_overlap(const GInterval &interv) const { return chromid == interv.chromid && Segment::do_overlap(interv); }

	// returns true if the two intervals overlap each other
	bool do_overlap(const Segment &segment) const { return Segment::do_overlap(segment); }

	// returns true if the interval fully overlaps the given interval
	bool do_contain(const GInterval &interv) const { return chromid == interv.chromid && Segment::do_contain(interv); }

	// returns intersection(this, interv) / union(this, interv)
	double coverage_ratio(const GInterval &interv) const;

	static char char2strand(char c);

	char *debug_str(const GenomeChromKey &chromkey) const {
		static char str[200];
		sprintf(str, "(%s, %ld, %ld)", chromkey.id2chrom(chromid).c_str(), start, end);
		return str;
	}

	char *debug_str() const {
		static char str[200];
		sprintf(str, "(%d, %ld, %ld)", chromid, start, end);
		return str;
	}
};


//----------------------------------- IMPLEMENTATION --------------------------------------------

inline bool GInterval::operator<(const GInterval &interv) const
{
	return chromid < interv.chromid || chromid == interv.chromid && start < interv.start;
}

template <class T> void *GInterval::cast2udata(T v)
{
	void *p = NULL;
	memcpy(&p, &v, sizeof(v));
	return p;
}

inline double GInterval::dist2center(int64_t coord) const
{
	if (coord < start || coord >= end)
		return numeric_limits<double>::quiet_NaN();

	if (strand == 1)
		return (double)(coord - (start + end) / 2);
	if (strand == -1)
		return (double)((start + end) / 2 - coord);
	return fabs(coord - (start + end) / 2);
}

inline void GInterval::verify(const GenomeChromKey &chromkey, bool check_chrom_boundary) const
{
	if (start < 0)
		TGLError<GInterval>(BAD_INTERVAL, "Interval (%s, %ld, %ld): start coordinate must be greater or equal than zero", chromkey.id2chrom(chromid).c_str(), start, end);
	if (start >= end)
		TGLError<GInterval>(BAD_INTERVAL, "Interval (%s, %ld, %ld): start coordinate must be lesser than end coordinate", chromkey.id2chrom(chromid).c_str(), start, end);
	if (check_chrom_boundary && (uint64_t)end > chromkey.get_chrom_size(chromid))
		TGLError<GInterval>(BAD_INTERVAL, "Interval (%s, %ld, %ld): end coordinate exceeds chromosome boundaries", chromkey.id2chrom(chromid).c_str(), start, end);
}

#endif /* GINTERVAL_H_ */

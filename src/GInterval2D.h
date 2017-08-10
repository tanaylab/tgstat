/*
 * GInterval2D.h
 *
 *  Created on: Jan 16, 2012
 *      Author: hoichman
 */

#ifndef GINTERVAL2D_H_
#define GINTERVAL2D_H_

#include <vector>

#include "GenomeChromKey.h"
#include "Rectangle.h"
#include "TGLException.h"

//------------------------------------- GInterval2D ----------------------------------------------
// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

class GInterval2D : public Rectangle {
public:
	enum Errors { BAD_INTERVAL };
	enum { CHROM1, START1, END1, CHROM2, START2, END2, NUM_COLS };

	static const char *COL_NAMES[];

	const int     &chromid1() const { return m_chromid1; }
	const int     &chromid2() const { return m_chromid2; }
	const int64_t &start1() const { return x1; }
	const int64_t &end1() const { return x2; }
	const int64_t &start2() const { return y1; }
	const int64_t &end2() const { return y2; }
	void *udata() const { return m_udata; }

	int     &chromid1() { return m_chromid1; }
	int     &chromid2() { return m_chromid2; }
	int64_t &start1() { return x1; }
	int64_t &end1() { return x2; }
	int64_t &start2() { return y1; }
	int64_t &end2() { return y2; }
	void   *&udata() { return m_udata; }

	GInterval2D() : Rectangle(-1, -1, -1, -1), m_chromid1(-1), m_chromid2(-1), m_udata(NULL) {}

	GInterval2D(int _chromid1, int _chromid2, const Rectangle &rect, void *_udata = NULL) :
		Rectangle(rect), m_chromid1(_chromid1), m_chromid2(_chromid2), m_udata(_udata) {}

	GInterval2D(int _chromid1, int64_t _start1, int64_t _end1, int _chromid2, int64_t _start2, int64_t _end2, void *_udata = NULL) :
		Rectangle(_start1, _start2, _end1, _end2), m_chromid1(_chromid1), m_chromid2(_chromid2), m_udata(_udata) {}

	void set(int _chromid1, int64_t _start1, int64_t _end1, int _chromid2, int64_t _start2, int64_t _end2, void *_udata = NULL);

	// verifies basic interval correctness
	void verify(const GenomeChromKey &chromkey) const;

	// compares two intervals by chrom1 and chrom2
	bool operator<(const GInterval2D &interv) const;

	bool           intersects_diagonal() const;

	double surface() const { return (end1() - start1()) * (end2() - start2()); }

	bool is_same_chrom(const GInterval2D &interv) const { return chromid1() == interv.chromid1() && chromid2() == interv.chromid2(); }

	char *debug_str() const {
		static char str[200];
		sprintf(str, "(%d, %ld, %ld) (%d, %ld, %ld)", chromid1(), start1(), end1(), chromid2(), start2(), end2());
		return str;
	}

	char *debug_str(const GenomeChromKey &chromkey) const {
		static char str[200];
		sprintf(str, "(%s, %ld, %ld) (%s, %ld, %ld)", chromkey.id2chrom(chromid1()).c_str(), start1(), end1(), chromkey.id2chrom(chromid2()).c_str(), start2(), end2());
		return str;
	}

private:
	int     m_chromid1;
	int     m_chromid2;
	void   *m_udata;
};


//------------------------------------- ChromPair ------------------------------------------------

struct ChromPair {
	ChromPair(int _chromid1, int _chromid2) : chromid1(_chromid1), chromid2(_chromid2) {}
	ChromPair(const ChromPair &obj) : chromid1(obj.chromid1), chromid2(obj.chromid2) {}

	bool operator==(const ChromPair &obj) const { return chromid1 == obj.chromid1 && chromid2 == obj.chromid2; }
	bool operator<(const ChromPair &obj) const { return chromid1 < obj.chromid1 || chromid1 == obj.chromid1 && chromid2 < obj.chromid2; }

	int chromid1;
	int chromid2;
};


//----------------------------------- IMPLEMENTATION --------------------------------------------

inline bool GInterval2D::operator<(const GInterval2D &interv) const
{
	return chromid1() < interv.chromid1() || chromid1() == interv.chromid1() && chromid2() < interv.chromid2();
}

inline void GInterval2D::set(int _chromid1, int64_t _start1, int64_t _end1, int _chromid2, int64_t _start2, int64_t _end2, void *_udata)
{
	x1 = _start1;
	y1 = _start2;
	x2 = _end1;
	y2 = _end2;
	m_chromid1 = _chromid1;
	m_chromid2 = _chromid2;
	m_udata = _udata;
}

inline bool GInterval2D::intersects_diagonal() const
{
	// Diagonal is intersected if and only if there is a point x,y such that x==y.
	// Therefore interval (x1,x2) must intersect (y1,y2) in 1D.
	return m_chromid1 == m_chromid2 && max(x1, y1) < min(x2, y2);
}

inline void GInterval2D::verify(const GenomeChromKey &chromkey) const
{
	if (start1() < 0 || start2() < 0)
		TGLError<GInterval2D>(BAD_INTERVAL, "Interval (%s, %ld, %ld, %s, %ld, %ld): coordinate must be greater or equal than zero",
				chromkey.id2chrom(chromid1()).c_str(), start1(), end1(), chromkey.id2chrom(chromid2()).c_str(), start2(), end2());

	if (start1() >= end1() || start2() >= end2())
		TGLError<GInterval2D>(BAD_INTERVAL, "Interval (%s, %ld, %ld, %s, %ld, %ld): start coordinate must be lesser than end coordinate",
				chromkey.id2chrom(chromid1()).c_str(), start1(), end1(), chromkey.id2chrom(chromid2()).c_str(), start2(), end2());

	if ((uint64_t)end1() > chromkey.get_chrom_size(chromid1()) || (uint64_t)end2() > chromkey.get_chrom_size(chromid2()))
		TGLError<GInterval2D>(BAD_INTERVAL, "Interval (%s, %ld, %ld, %s, %ld, %ld): coordinate exceeds the chromosome boundaries",
				chromkey.id2chrom(chromid1()).c_str(), start1(), end1(), chromkey.id2chrom(chromid2()).c_str(), start2(), end2());
}

#endif /* GINTERVAL2D_H_ */

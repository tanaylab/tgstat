#ifndef _GINTERVALS2D_H_INCLUDED_
#define _GINTERVALS2D_H_INCLUDED_

#include "GIntervalsFetcher2D.h"

//------------------------------------- GIntervals2D ---------------------------------------------
// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

class GIntervals2D : public std::vector<GInterval2D>, public GIntervalsFetcher2D {
public:
	GIntervals2D() : std::vector<GInterval2D>(), GIntervalsFetcher2D(INTERVALS2D), m_cur_chromid1(-1), m_cur_chromid2(-1), m_num_chroms(0) {}
	GIntervals2D(size_type n) : std::vector<GInterval2D>(n), GIntervalsFetcher2D(INTERVALS2D), m_cur_chromid1(-1), m_cur_chromid2(-1), m_num_chroms(0) {}
	GIntervals2D(const std::vector<GInterval2D> &v) : std::vector<GInterval2D>(v), GIntervalsFetcher2D(INTERVALS2D), m_cur_chromid1(-1), m_cur_chromid2(-1), m_num_chroms(0) {}

	void clear();

	//-------------------------------- GIntervalsFetcher2D interface -----------------------------------

	virtual GIntervalsFetcher2D *create_masked_copy(const set<ChromPair> &chrompairs_mask) const;

	virtual void seal();

	virtual size_t size() const { return std::vector<GInterval2D>::size(); }

	virtual size_t size(int chromid1, int chromid2) const;

	virtual int num_chrom_pairs() const;

	virtual double surface() const;                           // complexity: O(n)
	virtual double surface(int chromid1, int chromid2) const; // complexity: O(n)

	virtual void begin_iter();

	// intervals are expected to be sorted by chromosome!
	virtual void begin_chrom_iter(int chromid1, int chromid2);

	virtual bool next();
	virtual bool next_in_chrom();

	virtual bool isend() const { return m_iinterval >= end(); }
	virtual bool isend_chrom() const;

	virtual const_iterator get_chrom_begin() const;
	virtual const_iterator get_chrom_end() const;

	virtual bool get_next_chroms(int *chromid1, int *chromid2);

	virtual size_t iter_index() const { return m_iinterval - begin(); }

	virtual size_t iter_chrom_index() const { return m_iter_chrom_index; }

	virtual const GInterval2D &cur_interval() const { return *m_iinterval; }

	virtual void sort(Compare_t = compare_for_sort);

	virtual void verify_no_overlaps(const GenomeChromKey &chromkey, const char *error_prefix = "") const;

private:
	mutable const_iterator m_iinterval;
	size_t                 m_iter_chrom_index;
	int                    m_cur_chromid1;
	int                    m_cur_chromid2;
	mutable int            m_num_chroms;

	// Holds the start location of each chromosome. Chromosomes that do not appear among the intervals
	// point to the start location of the next closest existing chromosome.
	// This order ensures iter_index() function works when begin_chrom() is called.
	mutable vector<const_iterator> m_chrom2itr;

	int chroms2idx(int chromid1, int chromid2) const { return chromid1 * m_num_chroms + chromid2; }
	void build_chrom_map() const;
};


//------------------------------------- IMPLEMENTATION ------------------------------------------

inline void GIntervals2D::clear()
{
	std::vector<GInterval2D>::clear();
	seal();
}

inline GIntervalsFetcher2D *GIntervals2D::create_masked_copy(const set<ChromPair> &chrompairs_mask) const
{
	GIntervals2D *obj = new GIntervals2D();
	for (const_iterator iinterv = begin(); iinterv < end(); ++iinterv) {
		if (chrompairs_mask.find(ChromPair(iinterv->chromid1(), iinterv->chromid2())) != chrompairs_mask.end())
			obj->push_back(*iinterv);
	}
	obj->seal();
	return obj;
}

inline void GIntervals2D::seal()
{
	m_num_chroms = 0;
	m_chrom2itr.clear();
	m_iinterval = begin();
}

inline size_t GIntervals2D::size(int chromid1, int chromid2) const
{
	build_chrom_map();
	if (chromid1 >= m_num_chroms || chromid2 >= m_num_chroms) 
		return 0;

	int idx = chroms2idx(chromid1, chromid2);
	if ((size_t)idx == m_chrom2itr.size() - 1) 
		return end() - m_chrom2itr[idx];
	return m_chrom2itr[idx + 1] - m_chrom2itr[idx];
}

inline int GIntervals2D::num_chrom_pairs() const
{
	int num_chrom_pairs = 0;
	build_chrom_map();
	for (int chromid1 = 0; chromid1 < m_num_chroms; ++chromid1) {
		for (int chromid2 = 0; chromid2 < m_num_chroms; ++chromid2) {
			if (size(chromid1, chromid2)) 
				++num_chrom_pairs;
		}
	}
	return num_chrom_pairs;
}

inline void GIntervals2D::begin_iter()
{
	m_cur_chromid1 = -1;
	m_cur_chromid2 = -1;
	m_iinterval = begin();
	m_iter_chrom_index = 0;
}

inline bool GIntervals2D::next()
{
	++m_iinterval;
	bool retv = isend();
	if (!retv && m_iinterval->is_same_chrom(*(m_iinterval - 1)))
		++m_iter_chrom_index;
	else
		m_iter_chrom_index = 0;
	return !retv;
}

inline bool GIntervals2D::next_in_chrom()
{
	if (!isend_chrom())
		++m_iinterval;
	return !isend_chrom();
}

inline bool GIntervals2D::isend_chrom() const
{
	return m_cur_chromid1 < 0 || m_cur_chromid2 < 0 || m_iinterval >= end() ||
		m_iinterval->chromid1() != m_cur_chromid1 || m_iinterval->chromid2() != m_cur_chromid2;
}

inline bool GIntervals2D::get_next_chroms(int *chromid1, int *chromid2)
{
	build_chrom_map();
	if (*chromid2 < m_num_chroms - 1)
		++*chromid2;
	else {
		++*chromid1;
		*chromid2 = 0;
	}
	return *chromid1 < m_num_chroms && *chromid2 < m_num_chroms;
}

inline void GIntervals2D::build_chrom_map() const
{
	if (m_chrom2itr.empty() && size()) {
		m_num_chroms = 0;
		for (const_iterator iinterv = begin(); iinterv < end(); ++iinterv) {
			m_num_chroms = max(m_num_chroms, iinterv->chromid1() + 1);
			m_num_chroms = max(m_num_chroms, iinterv->chromid2() + 1);
		}

		m_chrom2itr.resize(m_num_chroms * m_num_chroms, end());

		for (const_iterator iinterv = begin(); iinterv != end(); ++iinterv) {
			int idx = chroms2idx(iinterv->chromid1(), iinterv->chromid2());
			if (m_chrom2itr[idx] == end())
				m_chrom2itr[idx] = iinterv;
		}

		// make non-existing chromosomes point to the start location of the next closest existing ones
		for (vector<const_iterator>::reverse_iterator ichrom2itr = m_chrom2itr.rbegin() + 1; ichrom2itr < m_chrom2itr.rend(); ++ichrom2itr) {
			if (*ichrom2itr == end())
				*ichrom2itr = *(ichrom2itr - 1);
			else if (*ichrom2itr > *(ichrom2itr - 1))
				TGLError<GIntervals2D>(UNSORTED_INTERVALS, "Intervals are not sorted");
		}
	}
}

#endif


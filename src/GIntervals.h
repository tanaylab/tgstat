#ifndef _GINTERVALS_H_INCLUDED_
#define _GINTERVALS_H_INCLUDED_

#include "GIntervalsFetcher1D.h"

//------------------------------------- GIntervals ----------------------------------------------
// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

class GIntervals : public std::vector<GInterval>, public GIntervalsFetcher1D {
public:
	GIntervals() : std::vector<GInterval>(), GIntervalsFetcher1D(INTERVALS1D), m_cur_chromid(-1) {}
	GIntervals(size_type n) : std::vector<GInterval>(n), GIntervalsFetcher1D(INTERVALS1D), m_cur_chromid(-1) {}
	GIntervals(const std::vector<GInterval> &v) : std::vector<GInterval>(v), GIntervalsFetcher1D(INTERVALS1D), m_cur_chromid(-1) {}

	bool empty() const { return std::vector<GInterval>::empty(); }

	void clear();

	// returns an interval that fully overlaps the given interval or NULL if such does not exist; intervals are expected to be sorted
	const GInterval *containing_interval(const GInterval &interv);

	// calculates union of two sets of intervals
	static void unify(const GIntervals &intervs1, const GIntervals &intervs2, GIntervals &res_intervs);

	// calculates intersection of two sets of intervals
	static void intersect(const GIntervals &intervs1, const GIntervals &intervs2, GIntervals &res_intervs);

	// calculates diff of two sets of intervals (set1-set2)
	static void diff(const GIntervals &intervs1, const GIntervals &intervs2, GIntervals &res_intervs);

	void read(const GenomeChromKey &chromkey, std::istream &tab, int nostrand=0); // throws TGLException<GInterval>
    void read_bed(const GenomeChromKey &chromkey, std::istream &bed); // throws TGLException<GInterval>
	void write(const GenomeChromKey &chromkey, std::ostream &tab);

	//-------------------------------- GIntervalsFetcher1D interface -----------------------------------

	virtual GIntervalsFetcher1D *create_masked_copy(const set<int> &chromids_mask) const;

	virtual void seal();

	virtual size_t size() const { return std::vector<GInterval>::size(); }

	virtual size_t size(int chromid) const;

	virtual int num_chroms() const;

	virtual int64_t range() const;              // complexity: O(n)
	virtual int64_t range(int chromid) const;   // complexity: O(n)

	virtual void begin_iter();

	// intervals are expected to be sorted by chromosome!
	virtual void begin_chrom_iter(int chromid);

	virtual bool next();
	virtual bool next_in_chrom();

	virtual bool isend() const { return m_iinterval >= end(); }
	virtual bool isend_chrom() const { return m_iinterval >= end() || m_iinterval->chromid != m_cur_chromid; }

	virtual const_iterator get_chrom_begin() const;
	virtual const_iterator get_chrom_end() const;

	virtual size_t iter_index() const { return m_iinterval - begin(); }

	virtual size_t iter_chrom_index() const { return m_iter_chrom_index; }

	virtual const GInterval &cur_interval() const { return *m_iinterval; }

	virtual void sort(Compare_t = compare_by_start_coord);

	virtual void unify_overlaps(bool unify_touching_intervals = true);

	virtual void verify_no_overlaps(const GenomeChromKey &chromkey, const char *error_prefix = "") const;

protected:
	mutable const_iterator m_iinterval;
	int                    m_cur_chromid;
	size_t                 m_iter_chrom_index;

	// Holds the start location of each chromosome. Chromosomes that do not appear among the intervals
	// point to the start location of the next closest existing chromosome.
	// This order ensures iter_index() function works when begin_chrom() is called.
	mutable vector<const_iterator> m_chrom2itr;

	void build_chrom_map() const;
};

//------------------------------------- IMPLEMENTATION ------------------------------------------

inline void GIntervals::clear()
{
	std::vector<GInterval>::clear();
	seal();
}

inline GIntervalsFetcher1D *GIntervals::create_masked_copy(const set<int> &chromids_mask) const
{
	GIntervals *obj = new GIntervals();
	for (const_iterator iinterv = begin(); iinterv < end(); ++iinterv) {
		if (chromids_mask.find(iinterv->chromid) != chromids_mask.end())
			obj->push_back(*iinterv);
	}
	obj->seal();
	return obj;
}

inline void GIntervals::seal()
{
	m_chrom2itr.clear();
	m_iinterval = begin();
}

inline size_t GIntervals::size(int chromid) const
{
	build_chrom_map();
	if ((size_t)chromid >= m_chrom2itr.size()) 
		return 0;
	if ((size_t)chromid == m_chrom2itr.size() - 1) 
		return end() - m_chrom2itr[chromid];
	return m_chrom2itr[chromid + 1] - m_chrom2itr[chromid];
}

inline int GIntervals::num_chroms() const
{
	int num_chroms = 0;
	build_chrom_map();
	for (size_t chromid = 0; chromid < m_chrom2itr.size(); ++chromid) {
		if (size(chromid)) 
			++num_chroms;
	}
	return num_chroms;
}

inline void GIntervals::begin_iter()
{
	m_cur_chromid = -1;
	m_iinterval = begin();
	m_iter_chrom_index = 0;
}

inline bool GIntervals::next()
{
	++m_iinterval;
	bool retv = isend();
	if (!retv && m_iinterval->chromid == (m_iinterval - 1)->chromid) 
		++m_iter_chrom_index;
	else
		m_iter_chrom_index = 0;
	return !retv;
}

inline bool GIntervals::next_in_chrom()
{
	if (!isend_chrom()) 
		++m_iinterval;
	return !isend_chrom();
}

inline void GIntervals::build_chrom_map() const
{
	if (m_chrom2itr.empty() && size()) {
		for (const_iterator iinterv = begin(); iinterv < end(); ++iinterv) {
			m_chrom2itr.resize(max((size_t)(iinterv->chromid + 1), m_chrom2itr.size()), end());
			if (m_chrom2itr[iinterv->chromid] == end())
				m_chrom2itr[iinterv->chromid] = iinterv;
		}

		// make non-existing chromosomes point to the start location of the next closest existing ones
		for (vector<const_iterator>::reverse_iterator ichrom2itr = m_chrom2itr.rbegin() + 1; ichrom2itr < m_chrom2itr.rend(); ++ichrom2itr) {
			if (*ichrom2itr == end())
				*ichrom2itr = *(ichrom2itr - 1);
			else if (*ichrom2itr > *(ichrom2itr - 1))
				TGLError<GIntervals>(UNSORTED_INTERVALS, "Intervals are not sorted");
		}
	}
}

#endif


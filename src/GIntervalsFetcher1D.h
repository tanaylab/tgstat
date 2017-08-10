#ifndef _GINTERVALSFETCHER1D_H_INCLUDED_
#define _GINTERVALSFETCHER1D_H_INCLUDED_

#include <set>

#include "GInterval.h"
#include "GIntervalsFetcher.h"

class GIntervalsFetcher1D : public GIntervalsFetcher<GInterval> {
public:
	enum Errors { OVERLAPPING_INTERVAL, UNSORTED_INTERVALS };

	typedef bool (*Compare_t)(const GInterval &, const GInterval &);

	GIntervalsFetcher1D(Type type) : GIntervalsFetcher<GInterval>(type) {}
	virtual ~GIntervalsFetcher1D() {}

	// creates a copy of the object but only for the specified chromosome
	GIntervalsFetcher1D *create_masked_copy(int chromid) const;

	// creates a copy of the object but only for the selected chromosomes
	virtual GIntervalsFetcher1D *create_masked_copy(const set<int> &chromids_mask) const = 0;

	// returns number of intervals in the intervals set
	virtual size_t size() const = 0;

	// returns number of intervals for the given chromosome (intervals must be sorted)
	virtual size_t size(int chromid) const = 0;

	virtual void begin_chrom_iter(int chromid) = 0;

	// Get iterators for the current chromosome (begin_iter or begin_chrom_iter must be called before)
	// in some implementations these two functions might throw "unsupported" exception
	virtual vector<GInterval>::const_iterator get_chrom_begin() const = 0;
	virtual vector<GInterval>::const_iterator get_chrom_end() const = 0;

	// returns number of unique chromosomes appearing in the intervals (intervals must be sorted)
	virtual int num_chroms() const = 0;

	// Returns sum of ranges of all intervals
	virtual int64_t range() const = 0;

	// Returns sum of ranges of all intervals in the chromosome
	virtual int64_t range(int chromid) const = 0;

	// Sorts the intervals
	virtual void sort(Compare_t = compare_by_start_coord) = 0;

	// unite overlapping intervals (intervs are expected to be already sorted)
	virtual void unify_overlaps(bool unify_touching_intervals = true) = 0;

	// Verifies that there are no overlaps between the intervals; if intervals overlap an exception is thrown.
	// Intervals are expected to be already sorted.
	virtual void verify_no_overlaps(const GenomeChromKey &chromkey, const char *error_prefix = "") const = 0;

	// compares two intervals by chrom and then by start coord, returns true if interv1<interv2
	static bool compare_by_start_coord(const GInterval &interv1, const GInterval &interv2);

	// compares two intervals by chrom and then by end coord, returns true if interv1<interv2
	static bool compare_by_end_coord(const GInterval &interv1, const GInterval &interv2);
};


//------------------------------------- IMPLEMENTATION ------------------------------------------

inline GIntervalsFetcher1D *GIntervalsFetcher1D::create_masked_copy(int chromid) const
{
	set<int> mask;
	mask.insert(chromid);
	return create_masked_copy(mask);
}

#endif


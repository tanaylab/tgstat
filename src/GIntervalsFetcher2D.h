#ifndef _GINTERVALSFETCHER2D_H_INCLUDED_
#define _GINTERVALSFETCHER2D_H_INCLUDED_

#include <set>

#include "GInterval2D.h"
#include "GIntervalsFetcher.h"

class GIntervalsFetcher2D : public GIntervalsFetcher<GInterval2D> {
public:
	enum Errors { OVERLAPPING_INTERVAL, UNSORTED_INTERVALS };

	typedef bool (*Compare_t)(const GInterval2D &, const GInterval2D &);

	GIntervalsFetcher2D(Type type) : GIntervalsFetcher<GInterval2D>(type) {}
	virtual ~GIntervalsFetcher2D() {}

	// creates a copy of the object but only for the specified chromosome pair
	virtual GIntervalsFetcher2D *create_masked_copy(int chromid1, int chromid2) const;

	// creates a copy of the object but only for the selected chromosome pairs
	virtual GIntervalsFetcher2D *create_masked_copy(const set<ChromPair> &chrompairs_mask) const = 0;

	// returns number of intervals in the intervals set
	virtual size_t size() const = 0;

	// returns number of intervals for the given chromosome pair (intervals must be sorted)
	virtual size_t size(int chromid1, int chromid2) const = 0;

	// returns number of unique chromosomes pairs appearing in the intervals
	virtual int num_chrom_pairs() const = 0;

	// Returns sum of intervals surfaces
	virtual double surface() const = 0;

	// Returns sum of intervals surfaces for given chromosomes
	virtual double surface(int chromid1, int chromid2) const = 0;

	virtual void begin_chrom_iter(int chromid1, int chromid2) = 0;

	// Get iterators for the current chromosome (begin_iter or begin_chrom_iter must be called before)
	// in some implementations these two functions might throw "unsupported" exception
	virtual vector<GInterval2D>::const_iterator get_chrom_begin() const = 0;
	virtual vector<GInterval2D>::const_iterator get_chrom_end() const = 0;

	// returns false if end reached and no more chrom pairs exist
	virtual bool get_next_chroms(int *chromid1, int *chromid2) = 0;

	// Sorts the intervals
	virtual void sort(Compare_t = compare_for_sort) = 0;

	// Verifies that there are no overlaps between the intervals; if intervals overlap an exception is thrown.
	// Intervals are expected to be already sorted.
	virtual void verify_no_overlaps(const GenomeChromKey &chromkey, const char *error_prefix = "") const = 0;

	// compares two intervals by chrom1 and chrom2
	static bool compare_for_sort(const GInterval2D &interv1, const GInterval2D &interv2);
};


//------------------------------------- IMPLEMENTATION ------------------------------------------

inline GIntervalsFetcher2D *GIntervalsFetcher2D::create_masked_copy(int chromid1, int chromid2) const
{
	set<ChromPair> mask;
	mask.insert(ChromPair(chromid1, chromid2));
	return create_masked_copy(mask);
}

#endif


/*
 * BinFinder.h
 *
 *  Created on: Dec 29, 2010
 *      Author: hoichman
 */

#ifndef BINFINDER_H_
#define BINFINDER_H_

#include <cmath>
#include <vector>

#include "TGLException.h"

// -------------------- BinFinder  -----------------------
// BinFinder is used to convert a value to bin given the division of break points.
// If the bin size is identical for all bins the complexity is O(1), otherwise it is O(logN).
//
// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

class BinFinder {
public:
	enum Errors { BAD_NUM_BREAKS, NOT_UNIQUE_BREAKS, UNSORTED_BREAKS };

	BinFinder() {}
	BinFinder(const std::vector<double> &breaks, bool include_lowest = false) { init(breaks, include_lowest); }
	BinFinder(const double *breaks, unsigned num_breaks, bool include_lowest = false) { init(breaks, num_breaks, include_lowest); }

	// Given x1, x2, x3, x4 as breaks the range is split to (x1,x2], (x2,x3], (x3,x4]
	// This is complient with R's "cut" function.
	// The breaks must be sorted and unique. There must be 2 or more breaks.
	void init(const std::vector<double> &breaks, bool include_lowest = false) { init(&*breaks.begin(), (unsigned)breaks.size(), include_lowest); }
	void init(const double *breaks, unsigned num_breaks, bool include_lowest = false);

	// returns the bin or -1 if the value is out of range
	int val2bin(double val) const;

	const std::vector<double> &get_breaks() const { return m_breaks; }
	unsigned                   get_numbins() const { return m_breaks.size() - 1; }

private:
	std::vector<double> m_breaks;
	double              m_binsize;  // 0 if binsize is not the same for all the bins
	bool                m_include_lowest;
};


// -------------------- implementation  -----------------------

inline int BinFinder::val2bin(double val) const
{
	if (m_include_lowest && val == m_breaks.front())
		return 0;

	if (std::isnan(val) || val <= m_breaks.front() || val > m_breaks.back())
		return -1;

	if (m_binsize) // are we using the same bin size for all bins?
		return std::min((int)ceil((val - m_breaks.front()) / (double)m_binsize) - 1, (int)get_numbins() - 1);

	// perform binary search
	unsigned start_bin = 0;
	unsigned end_bin = get_numbins();

	while (end_bin - start_bin > 1) {
		unsigned middle_bin = (start_bin + end_bin) / 2;

		if (val <= m_breaks[middle_bin])
			end_bin = middle_bin;
		else
			start_bin = middle_bin;
	}

	return start_bin;
}

#endif /* BINFINDER_H_ */

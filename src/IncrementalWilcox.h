#ifndef INCREMENTALWILCOX_H_
#define INCREMENTALWILCOX_H_

#include <map>

// IncrementalWilcox performs Wilcoxon test for two windows.
//
// Rather than performing the test for two windows given by containers, IncrementalWilcox
// works incrementally, i.e. values are added or deleted from the window one by one using update() function
// and the recomputed P-val is returned. The complexity of IncrementalWilcox::update() is O(N) where
// N = n1 + n2, and n1, n2 are the sizes of the windows. The complexity of computing Wilcoxon test from
// scratch is O(NlogN) at its best, hence IncrementalWilcox provides better complexity when used for sliding windows.

class IncrementalWilcox {
	std::map<double, unsigned> m_v2cnt[2];  // sorted unique values, each holding a counter of num appearances
	double                u[2];        // U as defined in Wilcox algorithm
	double                n[2];        // n as defined in Wilcox algorithm
	double                m_z;
	bool                  m_one_tailed;
	mutable double        m_pval;
	double                m_std_dev_u; // in fact this is 1 / [standard deviation of U]
	double                m_mean_u;

public:
	static const unsigned MIN_RELIABLE_WINSIZE;

	IncrementalWilcox(bool one_tailed);

	// use reset() if you want to reset the algorithm and clear the windows that were accumulated so far
	void reset();

	// old_v1 - value that should be deleted from window 1, use numeric_limits<double>::quiet_NaN() if you don't want to delete a value
	// new_v1 - value that should be added to window 1, use numeric_limits<double>::quiet_NaN() if you don't want to add a value
	void update(double old_v1, double new_v1, double old_v2, double new_v2);

	// P-val of the test, can be -1 if the test is unreliable
	double pval() const;

	// P-val of the test for finding peaks (v2 > v1), can be -1 if the test is unreliable
	double pval_highs() const { return m_pval != -1 && u[0] >= m_mean_u ? 0.5 : pval(); }

	// P-val of the test for finding lows (v2 < v1), can be -1 if the test is unreliable
	double pval_lows() const { return m_pval != -1 && u[1] >= m_mean_u ? 0.5 : pval(); }

	// z-score; can be 1 if the test is unreliable
	double z() const { return m_z; }

	// z-score for finding highs (v2 > v1); can be 1 if the test is unreliable
	double z_highs() const { return m_z != 1 && u[0] >= m_mean_u ? 0 : m_z; }

	// z-score for finding lows (v2 < v1); can be 1 if the test is unreliable
	double z_lows() const { return m_z != 1 && u[1] >= m_mean_u ? 0 : m_z; }

	unsigned n1() const { return (unsigned)n[0]; }
	unsigned n2() const { return (unsigned)n[1]; }

	double   u1() const { return u[0]; }
	double   u2() const { return u[1]; }
};

#endif /* INCREMENTALWILCOX_H_ */

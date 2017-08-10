/*
 * GenomeTrackSparse.h
 *
 *  Created on: May 24, 2011
 *      Author: hoichman
 */

#ifndef GENOMETRACKSPARSE_H_
#define GENOMETRACKSPARSE_H_

#include <math.h>
#include <vector>

#include "GenomeTrack1D.h"
#include "GIntervals.h"

#ifdef __GXX_EXPERIMENTAL_CXX0X__
#ifndef isnan
#define isnan ::isnan
#endif
#endif

// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

class GenomeTrackSparse : public GenomeTrack1D {
public:
	GenomeTrackSparse();

	virtual void read_interval(const GInterval &interval);

	void init_read(const char *filename, int chromid);
	void init_write(const char *filename, int chromid);

	void write_next_interval(const GInterval &interval, float val);

	const GIntervals &get_intervals();
	const vector<float> &get_vals();

protected:
	static const int RECORD_SIZE;

	GIntervals    m_intervals;
	vector<float> m_vals;
	bool          m_loaded;
	int64_t       m_num_records;
	GIntervals::const_iterator m_icur_interval;

	void read_file_into_mem();
	void calc_vals(const GInterval &interval);
	bool check_first_overlap(const GIntervals::const_iterator &iinterval1, const GInterval &interval2);
};


//------------------------------------ IMPLEMENTATION --------------------------------

inline bool GenomeTrackSparse::check_first_overlap(const GIntervals::const_iterator &iinterval1, const GInterval &interval2)
{
	return iinterval1->do_overlap(interval2) && (iinterval1 == m_intervals.begin() || !(iinterval1 - 1)->do_overlap(interval2));
}

inline void GenomeTrackSparse::calc_vals(const GInterval &interval)
{
	float num_vs = 0;
	double mean_square_sum = 0;
	float v;

	m_last_sum = 0;
	m_last_min = numeric_limits<float>::max();
	m_last_max = -numeric_limits<float>::max();

	for (GIntervals::const_iterator iinterv = m_icur_interval; iinterv != m_intervals.end(); ++iinterv) {
		if (!iinterv->do_overlap(interval))
			break;

		v = m_vals[iinterv - m_intervals.begin()];
		if (!isnan(v)) {
			m_last_sum += v;
			m_last_min = min(m_last_min, v);
			m_last_max = max(m_last_max, v);

			if (m_functions[STDDEV])
				mean_square_sum += v * v;

			if (m_use_quantile)
				m_sp.add(v);

			++num_vs;
		}
	}

	if (num_vs > 0)
		m_last_avg = m_last_nearest = m_last_sum / num_vs;
	else
		m_last_avg = m_last_nearest = m_last_min = m_last_max = m_last_sum = numeric_limits<float>::quiet_NaN();

	// we are calaculating unbiased standard deviation:
	// sqrt(sum((x-mean)^2) / (N-1)) = sqrt(sum(x^2)/(N-1) - N*(mean^2)/(N-1))
	if (m_functions[STDDEV])
		m_last_stddev = num_vs > 1 ? sqrt(mean_square_sum / (num_vs - 1) - (m_last_avg * (double)m_last_avg) * (num_vs / (num_vs - 1))) : numeric_limits<float>::quiet_NaN();
}

#endif /* GENOMETRACKSPARSE_H_ */

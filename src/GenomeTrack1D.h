/*
 * GenomeTrack1D.h
 *
 *  Created on: Jan 25, 2012
 *      Author: hoichman
 */

#ifndef GENOMETRACK1D_H_
#define GENOMETRACK1D_H_

#include "GenomeTrack.h"
#include "GInterval.h"
#include "StreamPercentiler.h"

// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

class GenomeTrack1D : public GenomeTrack {
public:
	enum Functions { AVG, MIN, MAX, NEAREST, STDDEV, SUM, NUM_FUNCS };

	virtual ~GenomeTrack1D() {}

	int get_chrom_id() const { return m_chromid; }

	virtual void read_interval(const GInterval &interval) = 0;

	void register_function(Functions func) { m_functions[func] = true; }
	void register_quantile(uint64_t rnd_sampling_buf_size, uint64_t lowest_vals_buf_size, uint64_t highest_vals_buf_size);

	float last_avg() const { return m_last_avg; }
	float last_min() const { return m_last_min; }
	float last_max() const { return m_last_max; }
	float last_nearest() const { return m_last_nearest; }
	float last_stddev() const { return m_last_stddev; }
	float last_sum() const { return m_last_sum; }
	float last_quantile(double percentile);

	const string &file_name() const { return m_bfile.file_name(); }

protected:
	vector<bool> m_functions;
	bool         m_use_quantile;
	int          m_chromid;

	float        m_last_avg;
	float        m_last_min;
	float        m_last_max;
	float        m_last_nearest;
	float        m_last_stddev;
	float        m_last_sum;
	StreamPercentiler<float> m_sp;

	GenomeTrack1D(Type type) : GenomeTrack(type), m_use_quantile(false) { m_functions.resize(NUM_FUNCS, false); }
};


//--------------------------------------------- IMPLEMENTATION -----------------------------------------------------------

inline void GenomeTrack1D::register_quantile(uint64_t rnd_sampling_buf_size, uint64_t lowest_vals_buf_size, uint64_t highest_vals_buf_size)
{
	m_sp.init(rnd_sampling_buf_size, lowest_vals_buf_size, highest_vals_buf_size);
	m_use_quantile = true;
}

inline float GenomeTrack1D::last_quantile(double percentile)
{
	if (m_sp.stream_size()) {
		bool is_estimated;
		return m_sp.get_percentile(percentile, is_estimated);
	}
	return numeric_limits<float>::quiet_NaN();
}

#endif /* GENOMETRACK1D_H_ */

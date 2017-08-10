#include <errno.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "TGLException.h"
#include "GenomeTrackFixedBin.h"

#ifdef __GXX_EXPERIMENTAL_CXX0X__
#ifndef isnan
#define isnan ::isnan
#endif
#endif

void GenomeTrackFixedBin::read_interval(const GInterval &interval)
{
	if (m_use_quantile)
		m_sp.reset();

	// optimization of the most common case when the expression iterator starts at 0 and steps by bin_size
	if (interval.start == m_cur_coord && interval.end == m_cur_coord + m_bin_size) {
		if (read_next_bin(m_last_avg)) {
			m_last_min = m_last_max = m_last_nearest = m_last_sum = m_last_avg;
			m_last_stddev = numeric_limits<float>::quiet_NaN();
			if (m_use_quantile && !isnan(m_last_avg))
				m_sp.add(m_last_avg);
		} else
			m_last_min = m_last_max = m_last_nearest = m_last_avg = m_last_stddev = m_last_sum = numeric_limits<float>::quiet_NaN();
		return;
	}

	int64_t sbin = (int64_t)(interval.start / m_bin_size);
	int64_t ebin = (int64_t)ceil(interval.end / (double)m_bin_size);

	if (ebin == sbin + 1) {
		goto_bin(sbin);
		if (read_next_bin(m_last_avg)) {
			m_last_min = m_last_max = m_last_nearest = m_last_sum = m_last_avg;
			m_last_stddev = numeric_limits<float>::quiet_NaN();
			if (m_use_quantile && !isnan(m_last_avg))
				m_sp.add(m_last_avg);
		} else
			m_last_min = m_last_max = m_last_nearest = m_last_avg = m_last_stddev = m_last_sum = numeric_limits<float>::quiet_NaN();
	} else {
		float num_vs = 0;
		double mean_square_sum = 0;
		float v;

		m_last_sum = 0;
		m_last_min = numeric_limits<float>::max();
		m_last_max = -numeric_limits<float>::max();

		goto_bin(sbin);
		for (int64_t bin = sbin; bin < ebin; ++bin) {
			if (read_next_bin(v) && !isnan(v)) {
				m_last_sum += v;
				m_last_min = min(m_last_min, v);
				m_last_max = max(m_last_max, v);

				if (m_functions[STDDEV])
					mean_square_sum += v * v;

				if (m_use_quantile && !isnan(v))
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
}

void GenomeTrackFixedBin::init_read(const char *filename, const char *mode, int chromid)
{
	m_cur_coord = 0;

	if (m_bfile.open(filename, mode))
		TGLError<GenomeTrackFixedBin>("%s", strerror(errno));

	if (m_bfile.read(&m_bin_size, sizeof(m_bin_size)) != sizeof(m_bin_size)) {
		if (m_bfile.error())
			TGLError<GenomeTrackFixedBin>("Failed to read a dense track file %s: %s", filename, strerror(errno));
		TGLError<GenomeTrackFixedBin>("Invalid format of a dense track file %s", filename);
	}

	// determine the number of samples in the file
	double num_samples = (m_bfile.file_size() - m_bfile.tell()) / (double)sizeof(float);

	if (m_bin_size <= 0 || num_samples != (int64_t)num_samples)
		TGLError<GenomeTrackFixedBin>("Invalid format of a dense track file %s", filename);

	m_num_samples = (int64_t)num_samples;
	m_chromid = chromid;
}

void GenomeTrackFixedBin::init_write(const char *filename, unsigned bin_size, int chromid)
{
	m_num_samples = 0;
	m_cur_coord = 0;

	umask(07);

	if (m_bfile.open(filename, "wb"))
		TGLError<GenomeTrackFixedBin>("Opening a dense track file %s: %s", filename, strerror(errno));

	m_bin_size = bin_size;
	if (m_bfile.write(&m_bin_size, sizeof(m_bin_size)) != sizeof(m_bin_size)) {
		if (m_bfile.error())
			TGLError<GenomeTrackFixedBin>("Failed to write a dense track file %s: %s", filename, strerror(errno));
		TGLError<GenomeTrackFixedBin>("Failed to write a dense track file %s", filename);
	}

	m_chromid = chromid;
}

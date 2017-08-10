/*
 * GenomeTrackFixedBin.h
 *
 *  Created on: May 15, 2011
 *      Author: hoichman
 */

#ifndef GENOMETRACKFIXEDBIN_H_
#define GENOMETRACKFIXEDBIN_H_

#include <math.h>

#include "GenomeTrack1D.h"

// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

class GenomeTrackFixedBin : public GenomeTrack1D {
public:
	GenomeTrackFixedBin() : GenomeTrack1D(FIXED_BIN), m_bin_size(0), m_num_samples(0), m_cur_coord(0) {}

	virtual void read_interval(const GInterval &interval);

	void init_read(const char *filename, int chromid) { init_read(filename, "rb", chromid); }
	void init_write(const char *filename, unsigned bin_size, int chromid);

	void init_update(const char *filename, int chromid) { init_read(filename, "rb+", chromid); }

	unsigned get_bin_size() const { return m_bin_size; }
	int64_t  get_num_samples() const { return m_num_samples; }

	void goto_bin(size_t bin);

	bool read_next_bin(float &val);
	void write_next_bin(float val);
	void write_next_bins(float *vals, size_t num_vals);

protected:
	unsigned  m_bin_size;
	int64_t   m_num_samples;
	int64_t   m_cur_coord;

	void init_read(const char *filename, const char *mode, int chromid);
};


//------------------------------ IMPLEMENTATION ------------------------------------

inline void GenomeTrackFixedBin::goto_bin(size_t bin)
{
	if (m_bfile.seek((long)(bin * sizeof(float) + sizeof(m_bin_size)), SEEK_SET))
		TGLError<GenomeTrackFixedBin>("Failed to seek a dense track file %s: %s", m_bfile.file_name().c_str(), strerror(errno));
	m_cur_coord = bin * m_bin_size;
}


inline bool GenomeTrackFixedBin::read_next_bin(float &val)
{
	if (m_bfile.read(&val, sizeof(val)) != sizeof(val)) {
		if (m_bfile.error())
			TGLError<GenomeTrackFixedBin>("Failed to read a dense track file %s: %s", m_bfile.file_name().c_str(), strerror(errno));
		return false;
	}

	if (isinf(val))
		val = numeric_limits<float>::quiet_NaN();

	m_cur_coord += m_bin_size;
	return true;
}

inline void GenomeTrackFixedBin::write_next_bin(float val)
{
	if (m_bfile.write(&val, sizeof(val)) != sizeof(val)) {
		if (m_bfile.error())
			TGLError<GenomeTrackFixedBin>("Failed to write a dense track file %s: %s", m_bfile.file_name().c_str(), strerror(errno));
		TGLError<GenomeTrackFixedBin>("Failed to write a dense track file %s", m_bfile.file_name().c_str());
	}
	m_num_samples++;
	m_cur_coord += m_bin_size;
}

inline void GenomeTrackFixedBin::write_next_bins(float *vals, size_t num_vals)
{
	size_t size = sizeof(vals[0]) * num_vals;
	if (m_bfile.write(vals, size) != size) {
		if (m_bfile.error())
			TGLError<GenomeTrackFixedBin>("Failed to write a dense track file %s: %s", m_bfile.file_name().c_str(), strerror(errno));
		TGLError<GenomeTrackFixedBin>("Failed to write a dense track file %s", m_bfile.file_name().c_str());
	}
	m_num_samples += num_vals;
	m_cur_coord += m_bin_size * num_vals;
}

#endif /* GENOMETRACKFIXEDBIN_H_ */

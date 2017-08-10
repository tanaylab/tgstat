#include <sys/stat.h>
#include <sys/types.h>

#include "GenomeTrackSparse.h"

const int GenomeTrackSparse::RECORD_SIZE = 2 * sizeof(int64_t) + sizeof(float);

GenomeTrackSparse::GenomeTrackSparse() :
	GenomeTrack1D(SPARSE),
	m_loaded(false),
	m_num_records(0)
{}

void GenomeTrackSparse::init_read(const char *filename, int chromid)
{
	m_bfile.close();
	m_loaded = false;

	read_type(filename);

	// determine the number of records in the file
	double num_records = (m_bfile.file_size() - m_bfile.tell()) / (double)RECORD_SIZE;

	if (num_records != (int64_t)num_records)
		TGLError<GenomeTrackSparse>("Invalid format of a sparse track file %s", filename);

	m_num_records = (int64_t)num_records;

	m_chromid = chromid;
}

void GenomeTrackSparse::init_write(const char *filename, int chromid)
{
	m_bfile.close();
	m_loaded = false;
	write_type(filename);
	m_chromid = chromid;
}

void GenomeTrackSparse::read_file_into_mem()
{
	if (m_loaded)
		return;

	m_intervals.resize(m_num_records);
	m_vals.resize(m_num_records);

	for (int64_t i = 0; i < m_num_records; ++i) {
		GInterval &interval = m_intervals[i];

		if (m_bfile.read(&interval.start, sizeof(int64_t)) != sizeof(int64_t) ||
				m_bfile.read(&interval.end, sizeof(int64_t)) != sizeof(int64_t) ||
				m_bfile.read(&m_vals[i], sizeof(float)) != sizeof(float))
		{
			if (m_bfile.error())
				TGLError<GenomeTrackSparse>("Failed to read a sparse track file %s: %s", m_bfile.file_name().c_str(), strerror(errno));
			TGLError<GenomeTrackSparse>("Invalid format of a sparse track file %s", m_bfile.file_name().c_str());
		}

		if (isinf(m_vals[i]))
			m_vals[i] = numeric_limits<float>::quiet_NaN();

		interval.chromid = m_chromid;

		if (interval.start < 0 || interval.start >= interval.end || i && interval.start < m_intervals[i - 1].end)
			TGLError<GenomeTrackSparse>("Invalid format of a sparse track file %s", m_bfile.file_name().c_str());
	}

	m_icur_interval = m_intervals.begin();
	m_loaded = true;
}

void GenomeTrackSparse::read_interval(const GInterval &interval)
{
	m_last_avg = m_last_nearest = m_last_min = m_last_max = m_last_stddev = m_last_sum = numeric_limits<float>::quiet_NaN();

	if (m_use_quantile)
		m_sp.reset();

	read_file_into_mem();

	if (m_intervals.empty())
		return;

	if (m_intervals.front().start >= interval.end) {
		m_last_nearest = m_vals.front();
		return;
	}

	if (m_intervals.back().end <= interval.start) {
		m_last_nearest = m_vals.back();
		return;
	}

	if (check_first_overlap(m_icur_interval, interval)) {
		calc_vals(interval);
	} else if (m_icur_interval + 1 < m_intervals.end() && check_first_overlap(m_icur_interval + 1, interval)) {
		++m_icur_interval;
		calc_vals(interval);
	} else {
		// run the binary search
		GIntervals::const_iterator istart_interval = m_intervals.begin();
		GIntervals::const_iterator iend_interval = m_intervals.end();

		while (iend_interval - istart_interval > 1) {
			GIntervals::const_iterator imid_interval = istart_interval + (iend_interval - istart_interval) / 2;

			if (check_first_overlap(imid_interval, interval)) {
				m_icur_interval = imid_interval;
				calc_vals(interval);
				break;
			}

			// is mid_interval < interval?
			if (GIntervals::compare_by_start_coord(*imid_interval, interval))
				istart_interval = imid_interval;
			else
				iend_interval = imid_interval;
		}

		if (iend_interval - istart_interval == 1 && check_first_overlap(istart_interval, interval)) {
			m_icur_interval = istart_interval;
			calc_vals(interval);
		}

		if (iend_interval - istart_interval == 1)
			m_last_nearest = iend_interval == m_intervals.end() || interval.dist2interv(*istart_interval) <= interval.dist2interv(*iend_interval) ?
					m_vals[istart_interval - m_intervals.begin()] : m_vals[iend_interval - m_intervals.begin()];
	}
}

void GenomeTrackSparse::write_next_interval(const GInterval &interval, float val)
{
	size_t size = 0;
	size += m_bfile.write(&interval.start, sizeof(interval.start));
	size += m_bfile.write(&interval.end, sizeof(interval.end));
	size += m_bfile.write(&val, sizeof(val));

	if ((int)size != RECORD_SIZE) {
		if (m_bfile.error())
			TGLError<GenomeTrackSparse>("Failed to write a sparse track file %s: %s", m_bfile.file_name().c_str(), strerror(errno));
		TGLError<GenomeTrackSparse>("Failed to write a sparse track file %s", m_bfile.file_name().c_str());
	}
}

const GIntervals &GenomeTrackSparse::get_intervals()
{
	read_file_into_mem();
	return m_intervals;
}

const vector<float> &GenomeTrackSparse::get_vals()
{
	read_file_into_mem();
	return m_vals;
}

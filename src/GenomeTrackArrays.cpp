#include <sys/stat.h>
#include <sys/types.h>

#include "GenomeTrackArrays.h"

#ifdef __GXX_EXPERIMENTAL_CXX0X__
#ifndef isnan
#define isnan ::isnan
#endif
#endif

const char *GenomeTrackArrays::SLICE_FUNCTION_NAMES[GenomeTrackArrays::NUM_S_FUNCS] = { "avg", "min", "max", "stddev", "sum", "quantile" };
const int GenomeTrackArrays::RECORD_SIZE = 2 * sizeof(int64_t) + sizeof(long);

GenomeTrackArrays::GenomeTrackArrays() :
	GenomeTrack1D(ARRAYS),
	m_master_obj(NULL),
	m_last_array_vals_idx((size_t)-1),
	m_loaded(false),
	m_is_writing(false),
	m_slice_function(S_AVG),
	m_slice_percentile(0.5)
{
	m_dependent_objs.push_back(this);
}

void GenomeTrackArrays::set_master_obj(GenomeTrackArrays *master_obj)
{
	m_master_obj = master_obj;
	m_master_obj->m_dependent_objs.push_back(this);
}

void GenomeTrackArrays::init_read(const char *filename, int chromid)
{
	finish_writing();
	m_bfile.close();
	m_loaded = false;
	m_is_writing = false;
	m_intervals.clear();
	m_vals_pos.clear();
	if (!m_master_obj) 
		read_type(filename);
	m_chromid = chromid;
}

void GenomeTrackArrays::init_write(const char *filename, int chromid)
{
	finish_writing();
	m_bfile.close();
	m_loaded = false;
	m_is_writing = true;
	m_intervals.clear();
	m_vals_pos.clear();
	write_type(filename);

	m_intervals_pos = m_bfile.tell();
	if (m_bfile.write(&m_intervals_pos, sizeof(m_intervals_pos)) != sizeof(m_intervals_pos)) {
		if (m_bfile.error())
			TGLError<GenomeTrackArrays>("Failed to write %s track file %s: %s", TYPE_NAMES[ARRAYS], m_bfile.file_name().c_str(), strerror(errno));
		TGLError<GenomeTrackArrays>("Failed to write %s track file %s", TYPE_NAMES[ARRAYS], m_bfile.file_name().c_str());
	}

	m_chromid = chromid;
}

void GenomeTrackArrays::finish_writing()
{
	if (!m_is_writing) 
		return;

	m_is_writing = false;

	// write the position of intervals in the file
	m_bfile.seek(m_intervals_pos, SEEK_SET);
	m_intervals_pos = m_bfile.file_size();
	m_bfile.write(&m_intervals_pos, sizeof(m_intervals_pos));

	// write the number of intervals
	m_bfile.seek(m_intervals_pos, SEEK_SET);
	size_t num_intervals = m_intervals.size();
	m_bfile.write(&num_intervals, sizeof(num_intervals));

	// write the intervals
	for (GIntervals::const_iterator iinterv = m_intervals.begin(); iinterv != m_intervals.end(); ++iinterv) {
		size_t size = 0;
		size += m_bfile.write(&iinterv->start, sizeof(iinterv->start));
		size += m_bfile.write(&iinterv->end, sizeof(iinterv->end));
		size += m_bfile.write(&m_vals_pos[iinterv - m_intervals.begin()], sizeof(m_vals_pos.front()));

		if ((int)size != RECORD_SIZE) {
			if (m_bfile.error())
				TGLError<GenomeTrackArrays>("Failed to write %s track file %s: %s", TYPE_NAMES[ARRAYS], m_bfile.file_name().c_str(), strerror(errno));
			TGLError<GenomeTrackArrays>("Failed to write %s track file %s", TYPE_NAMES[ARRAYS], m_bfile.file_name().c_str());
		}
	}
}

void GenomeTrackArrays::read_intervals_map()
{
	if (m_loaded)
		return;

	if (m_master_obj) {
		m_intervals = m_master_obj->m_intervals;
		m_vals_pos = m_master_obj->m_vals_pos;
		m_intervals_pos = m_master_obj->m_intervals_pos;
	} else {
		// read intervals position in the file (intervals appear after the values)
		if (m_bfile.read(&m_intervals_pos, sizeof(m_intervals_pos)) != sizeof(m_intervals_pos)) {
			if (m_bfile.error())
				TGLError<GenomeTrackArrays>("Failed to read %s track file %s: %s", TYPE_NAMES[ARRAYS], m_bfile.file_name().c_str(), strerror(errno));
			TGLError<GenomeTrackArrays>("Invalid format of %s track file %s", TYPE_NAMES[ARRAYS], m_bfile.file_name().c_str());
		}

		// read number of intervals
		if (m_bfile.seek(m_intervals_pos, SEEK_SET))
			TGLError<GenomeTrackArrays>("Failed to read %s track file %s: %s", TYPE_NAMES[ARRAYS], m_bfile.file_name().c_str(), strerror(errno));

		size_t num_intervals;
		if (m_bfile.read(&num_intervals, sizeof(num_intervals)) != sizeof(num_intervals)) {
			if (m_bfile.error())
				TGLError<GenomeTrackArrays>("Failed to read %s track file %s: %s", TYPE_NAMES[ARRAYS], m_bfile.file_name().c_str(), strerror(errno));
			TGLError<GenomeTrackArrays>("Invalid format of %s track file %s", TYPE_NAMES[ARRAYS], m_bfile.file_name().c_str());
		}

		// read the intervals
		m_intervals.resize(num_intervals);
		m_vals_pos.resize(num_intervals);

		for (int64_t i = 0; i < num_intervals; ++i) {
			GInterval &interval = m_intervals[i];

			if (m_bfile.read(&interval.start, sizeof(int64_t)) != sizeof(int64_t) ||
				m_bfile.read(&interval.end, sizeof(int64_t)) != sizeof(int64_t) ||
				m_bfile.read(&m_vals_pos[i], sizeof(long)) != sizeof(long))
			{
				if (m_bfile.error())
					TGLError<GenomeTrackArrays>("Failed to read %s track file %s: %s", TYPE_NAMES[ARRAYS], m_bfile.file_name().c_str(), strerror(errno));
				TGLError<GenomeTrackArrays>("Invalid format of %s track file %s", TYPE_NAMES[ARRAYS], m_bfile.file_name().c_str());
			}

			interval.chromid = m_chromid;

			if (interval.start < 0 || interval.start >= interval.end || i && interval.start < m_intervals[i - 1].end ||
				m_vals_pos[i] < 0 || m_vals_pos[i] >= m_bfile.file_size() || i && m_vals_pos[i - 1] >= m_vals_pos[i])
				TGLError<GenomeTrackArrays>("Invalid format of %s track file %s", TYPE_NAMES[ARRAYS], m_bfile.file_name().c_str());
		}
	}

	m_icur_interval = m_intervals.begin();
	m_loaded = true;
}

void GenomeTrackArrays::read_interval(const GInterval &interval)
{
	if (m_master_obj) 
		return;

	for (vector<GenomeTrackArrays *>::iterator itrack = m_dependent_objs.begin(); itrack != m_dependent_objs.end(); ++itrack) {
		(*itrack)->m_last_avg = (*itrack)->m_last_nearest = (*itrack)->m_last_min = (*itrack)->m_last_max = (*itrack)->m_last_stddev = (*itrack)->m_last_sum = numeric_limits<float>::quiet_NaN();

		if ((*itrack)->m_use_quantile)
			(*itrack)->m_sp.reset();
	}

	if (!m_loaded)
		read_intervals_map();

	if (m_intervals.empty())
		return;

	if (m_intervals.front().start >= interval.end) {
		for (vector<GenomeTrackArrays *>::iterator itrack = m_dependent_objs.begin(); itrack != m_dependent_objs.end(); ++itrack) {
			if ((*itrack)->m_functions[NEAREST])
				(*itrack)->m_last_nearest = (*itrack)->get_sliced_val(0);
		}
		return;
	}

	if (m_intervals.back().end <= interval.start) {
		for (vector<GenomeTrackArrays *>::iterator itrack = m_dependent_objs.begin(); itrack != m_dependent_objs.end(); ++itrack) {
			if ((*itrack)->m_functions[NEAREST])
				(*itrack)->m_last_nearest = (*itrack)->get_sliced_val(m_vals_pos.size() - 1);
		}
		return;
	}

	if (check_first_overlap(m_icur_interval, interval))
		calc_vals(interval);
	else if (m_icur_interval + 1 < m_intervals.end() && check_first_overlap(m_icur_interval + 1, interval)) {
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

		if (iend_interval - istart_interval == 1) {
			for (vector<GenomeTrackArrays *>::iterator itrack = m_dependent_objs.begin(); itrack != m_dependent_objs.end(); ++itrack) {
				if ((*itrack)->m_functions[NEAREST])
					(*itrack)->m_last_nearest = iend_interval == m_intervals.end() || interval.dist2interv(*istart_interval) <= interval.dist2interv(*iend_interval) ?
						(*itrack)->get_sliced_val(istart_interval - m_intervals.begin()) : (*itrack)->get_sliced_val(iend_interval - m_intervals.begin());
			}
		}
	}
}

void GenomeTrackArrays::write_next_interval(const GInterval &interval, const ArrayVals::const_iterator &iarray_vals_begin, const ArrayVals::const_iterator &iarray_vals_end)
{
	if (iarray_vals_begin == iarray_vals_end) 
		return;

	// add the intervals to the array: intervals are written later in finish_writing()
	m_intervals.push_back(interval);
	m_vals_pos.push_back(m_bfile.tell());

	// write the values
	unsigned num_non_nan_vals = 0;
	for (ArrayVals::const_iterator iarray_val = iarray_vals_begin; iarray_val != iarray_vals_end; ++iarray_val) {
		if (!isnan(iarray_val->val))
			++num_non_nan_vals;
	}

	m_bfile.write(&num_non_nan_vals, sizeof(num_non_nan_vals));

	for (ArrayVals::const_iterator iarray_val = iarray_vals_begin; iarray_val != iarray_vals_end; ++iarray_val) {
		if (!isnan(iarray_val->val)) {
			m_bfile.write(&iarray_val->val, sizeof(iarray_val->val));
			m_bfile.write(&iarray_val->idx, sizeof(iarray_val->idx));
		}
	}

	if (m_bfile.error())
		TGLError<GenomeTrackArrays>("Failed to write %s track file %s: %s", TYPE_NAMES[ARRAYS], m_bfile.file_name().c_str(), strerror(errno));
}

const GIntervals &GenomeTrackArrays::get_intervals()
{
	read_intervals_map();
	return m_intervals;
}

void GenomeTrackArrays::read_array_vals(size_t idx)
{
	if (m_last_array_vals_idx != idx) {
		m_last_array_vals_idx = idx;
		m_bfile.seek(m_vals_pos[idx], SEEK_SET);

		unsigned num_vals = 0;
		m_bfile.read(&num_vals, sizeof(num_vals));

		m_array_vals.resize(num_vals);

		for (unsigned i = 0; i < num_vals; ++i) {
			ArrayVal &array_val = m_array_vals[i];

			m_bfile.read(&array_val.val, sizeof(array_val.val));
			if (m_bfile.read(&array_val.idx, sizeof(array_val.idx)) != sizeof(array_val.idx)) {
				if (m_bfile.error())
					TGLError<GenomeTrackArrays>("Failed to read %s track file %s: %s", TYPE_NAMES[ARRAYS], m_bfile.file_name().c_str(), strerror(errno));
				TGLError<GenomeTrackArrays>("Invalid format of %s track file %s", TYPE_NAMES[ARRAYS], m_bfile.file_name().c_str());
			}
		}

		for (vector<GenomeTrackArrays *>::iterator itrack = m_dependent_objs.begin() + 1; itrack < m_dependent_objs.end(); ++itrack) 
			(*itrack)->m_array_vals = m_array_vals;
	}
}

void GenomeTrackArrays::get_sliced_vals(GIntervals::const_iterator iinterval, vector<float> &vals, unsigned numcols)
{
	if (m_master_obj)
		m_master_obj->read_array_vals(iinterval - m_intervals.begin());
	else 
		read_array_vals(iinterval - m_intervals.begin());

	vals.resize(0);

	if (m_slice.empty()) {
		vals.resize(numcols, numeric_limits<float>::quiet_NaN());
		for (ArrayVals::const_iterator iarray_val = m_array_vals.begin(); iarray_val != m_array_vals.end(); ++iarray_val) {
			if (iarray_val->idx >= numcols) 
				TGLError<GenomeTrackArrays>("Track file %s: value index %d exceeds total number of columns %d",
											m_bfile.file_name().c_str(), iarray_val->idx, numcols);
			vals[iarray_val->idx] = iarray_val->val;
		}
	} else {
		for (size_t islice = 0; islice < m_slice.size(); ++islice)
			vals.push_back(get_array_val(islice));
	}
}

float GenomeTrackArrays::get_sliced_val(size_t idx)
{
	if (m_master_obj)
		m_master_obj->read_array_vals(idx);
	else 
		read_array_vals(idx);

	// calculate sliced value for all available array values
	if (m_slice.empty()) {
		switch (m_slice_function) {
		case S_AVG:
			{
				double sum = 0;
				for (ArrayVals::const_iterator iarray_val = m_array_vals.begin(); iarray_val != m_array_vals.end(); ++iarray_val)
					sum += iarray_val->val;
				return sum / m_array_vals.size();
			}
		case S_MIN:
			{
				float s_min = numeric_limits<float>::max();
				for (ArrayVals::const_iterator iarray_val = m_array_vals.begin(); iarray_val != m_array_vals.end(); ++iarray_val)
					s_min = min(iarray_val->val, s_min);
				return s_min;
			}
		case S_MAX:
			{
				float s_max = -numeric_limits<float>::max();
				for (ArrayVals::const_iterator iarray_val = m_array_vals.begin(); iarray_val != m_array_vals.end(); ++iarray_val)
					s_max = max(iarray_val->val, s_max);
				return s_max;
			}
		case S_STDDEV:
			{ 
				// we are calaculating unbiased standard deviation:
				// sqrt(sum((x-mean)^2) / (N-1)) = sqrt(sum(x^2)/(N-1) - N*(mean^2)/(N-1))
				if (m_array_vals.size() <= 1) 
					return numeric_limits<float>::quiet_NaN();

				long double mean_square_sum = 0;
				double sum = 0;

				for (ArrayVals::const_iterator iarray_val = m_array_vals.begin(); iarray_val != m_array_vals.end(); ++iarray_val) {
					sum += iarray_val->val;
					mean_square_sum += iarray_val->val * iarray_val->val;
				}

				double N = m_array_vals.size();
				double avg = sum / N;
				return sqrt(mean_square_sum / (N - 1) - (avg * avg) * (N / (N - 1)));
			}
		case S_SUM:
			{
				double sum = 0;
				for (ArrayVals::const_iterator iarray_val = m_array_vals.begin(); iarray_val != m_array_vals.end(); ++iarray_val)
					sum += iarray_val->val;
				return sum;
			}
		case S_QUANTILE:
			{
				m_slice_sp.reset();
				for (ArrayVals::const_iterator iarray_val = m_array_vals.begin(); iarray_val != m_array_vals.end(); ++iarray_val)
					m_slice_sp.add(iarray_val->val);
				bool is_estimated;
				return m_slice_sp.get_percentile(m_slice_percentile, is_estimated);
			}
		default:
			TGLError<GenomeTrackArrays>("Unrecognized slice function");
		}
	}

	// calculate sliced value for specified array values
	switch (m_slice_function) {
	case S_AVG:
		{
			double sum = 0;
			double N = 0;
			for (size_t islice = 0; islice < m_slice.size(); ++islice) {
				float v = get_array_val(islice);
				if (!isnan(v)) {
					sum += v;
					++N;
				}
				if (m_array_hints[islice] >= m_array_vals.size())
					break;
			}
			return N ? sum / N : numeric_limits<float>::quiet_NaN();
		}
	case S_MIN:
		{
			float s_min = numeric_limits<float>::max();
			for (size_t islice = 0; islice < m_slice.size(); ++islice) {
				float v = get_array_val(islice);
				if (!isnan(v))
					s_min = min(v, s_min);
				if (m_array_hints[islice] >= m_array_vals.size())
					break;
			}
			return s_min == numeric_limits<float>::max() ? numeric_limits<float>::quiet_NaN() : s_min;
		}
	case S_MAX:
		{
			float s_max = -numeric_limits<float>::max();
			for (size_t islice = 0; islice < m_slice.size(); ++islice) {
				float v = get_array_val(islice);
				if (!isnan(v))
					s_max = max(v, s_max);
				if (m_array_hints[islice] >= m_array_vals.size())
					break;
			}
			return s_max == -numeric_limits<float>::max() ? numeric_limits<float>::quiet_NaN() : s_max;
		}
	case S_STDDEV:
		{ 
			// we are calaculating unbiased standard deviation:
			// sqrt(sum((x-mean)^2) / (N-1)) = sqrt(sum(x^2)/(N-1) - N*(mean^2)/(N-1))
			double mean_square_sum = 0;
			double sum = 0;
			double N = 0;

			for (size_t islice = 0; islice < m_slice.size(); ++islice) {
				float v = get_array_val(islice);
				if (!isnan(v)) {
					++N;
					sum += v;
					mean_square_sum += v * v;
				}
				if (m_array_hints[islice] >= m_array_vals.size())
					break;
			}

			if (N <= 1) 
				numeric_limits<float>::quiet_NaN();

			double avg = sum / N;
			return sqrt(mean_square_sum / (N - 1) - (avg * avg) * (N / (N - 1)));
		}
	case S_SUM:
		{
			double sum = 0;
			double N = 0;
			for (size_t islice = 0; islice < m_slice.size(); ++islice) {
				float v = get_array_val(islice);
				if (!isnan(v)) {
					sum += v;
					++N;
				}
				if (m_array_hints[islice] >= m_array_vals.size())
					break;
			}
			return N ? sum : numeric_limits<float>::quiet_NaN();
		}
	case S_QUANTILE:
		{
			m_slice_sp.reset();
			for (size_t islice = 0; islice < m_slice.size(); ++islice) {
				float v = get_array_val(islice);
				if (!isnan(v))
					m_slice_sp.add(v);
				if (m_array_hints[islice] >= m_array_vals.size())
					break;
			}

			if (m_slice_sp.stream_size()) {
				bool is_estimated;
				return m_slice_sp.get_percentile(m_slice_percentile, is_estimated);
			}

			return numeric_limits<float>::quiet_NaN();
		}
	default:
		TGLError<GenomeTrackArrays>("Unrecognized slice function");
	}
}


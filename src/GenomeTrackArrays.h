#ifndef GENOMETRACKARRAYS_H_
#define GENOMETRACKARRAYS_H_

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

class GenomeTrackArrays : public GenomeTrack1D {
public:
	enum SliceFunctions { S_AVG, S_MIN, S_MAX, S_STDDEV, S_SUM, S_QUANTILE, NUM_S_FUNCS };

	static const char *SLICE_FUNCTION_NAMES[NUM_S_FUNCS];

	struct ArrayVal {
		float    val;
		unsigned idx;

		ArrayVal() {}
		ArrayVal(float _val, unsigned _idx) : val(_val), idx(_idx) {}

		bool operator<(const ArrayVal &obj) const { return idx < obj.idx; }
	};

	struct LessIdx : public binary_function<ArrayVal, unsigned, bool> {
		bool operator()(const ArrayVal &obj, unsigned idx) const { return obj.idx < idx; }
	};

	typedef vector<ArrayVal> ArrayVals;

	GenomeTrackArrays();
	virtual ~GenomeTrackArrays() { finish_writing(); }

	// Used for optimization in TrackExpressionVars to prevent multiple virtual tracks read the same file over and over.
	// It optimizes read_interval: the file is read only once and then all objects that are depending on the master calculate
	// their track value.
	void set_master_obj(GenomeTrackArrays *master_obj);

	void set_slice(const vector<unsigned> &slice);
	void set_slice_function(SliceFunctions func, const vector<unsigned> &slice);
	void set_slice_quantile(double percentile, uint64_t rnd_sampling_buf_size, uint64_t lowest_vals_buf_size,
							uint64_t highest_vals_buf_size, const vector<unsigned> &slice);

	virtual void read_interval(const GInterval &interval);

	void init_read(const char *filename, int chromid);
	void init_write(const char *filename, int chromid);

	void write_next_interval(const GInterval &interval, const ArrayVals::const_iterator &iarray_vals_begin, const ArrayVals::const_iterator &iarray_vals_end);

	const GIntervals &get_intervals();
	void get_sliced_vals(GIntervals::const_iterator iinterval, vector<float> &vals, unsigned numcols);
	
protected:
	static const int RECORD_SIZE;

	GenomeTrackArrays          *m_master_obj;
	vector<GenomeTrackArrays *> m_dependent_objs;
	size_t                      m_last_array_vals_idx;

	GIntervals                  m_intervals;
	vector<long>                m_vals_pos;
	bool                        m_loaded;
	bool                        m_is_writing;
	long                        m_intervals_pos;
	GIntervals::const_iterator  m_icur_interval;

	SliceFunctions              m_slice_function;
	double                      m_slice_percentile;
	vector<unsigned>            m_slice;
	vector<unsigned>            m_array_hints;
	StreamPercentiler<float>    m_slice_sp;
	ArrayVals                   m_array_vals;

	float                       m_num_vs;
	double                      m_mean_square_sum;

	void read_intervals_map();
	void finish_writing();
	void read_array_vals(size_t idx);
	float get_array_val(size_t islice);
	float get_sliced_val(size_t idx);
	void calc_vals(const GInterval &interval);
	bool check_first_overlap(const GIntervals::const_iterator &iinterval1, const GInterval &interval2);
};


//------------------------------------ IMPLEMENTATION --------------------------------

inline bool GenomeTrackArrays::check_first_overlap(const GIntervals::const_iterator &iinterval1, const GInterval &interval2)
{
	return iinterval1->do_overlap(interval2) && (iinterval1 == m_intervals.begin() || !(iinterval1 - 1)->do_overlap(interval2));
}

inline void GenomeTrackArrays::set_slice(const vector<unsigned> &slice)
{
	m_slice = slice;
	m_array_hints.resize(m_slice.size(), 0);
}

inline void GenomeTrackArrays::set_slice_function(SliceFunctions func, const vector<unsigned> &slice)
{
	m_slice_function = func;
	set_slice(slice);
}

inline void GenomeTrackArrays::set_slice_quantile(double percentile, uint64_t rnd_sampling_buf_size, uint64_t lowest_vals_buf_size,
												  uint64_t highest_vals_buf_size, const vector<unsigned> &slice)
{
	m_slice_function = S_QUANTILE;
	m_slice_percentile = percentile;
	m_slice_sp.init(rnd_sampling_buf_size, lowest_vals_buf_size, highest_vals_buf_size);
	set_slice(slice);
}

inline void GenomeTrackArrays::calc_vals(const GInterval &interval)
{
	float v;

	for (vector<GenomeTrackArrays *>::iterator itrack = m_dependent_objs.begin(); itrack != m_dependent_objs.end(); ++itrack) {
		(*itrack)->m_num_vs = 0;
		(*itrack)->m_mean_square_sum = 0;
		(*itrack)->m_last_sum = 0;
		(*itrack)->m_last_min = numeric_limits<float>::max();
		(*itrack)->m_last_max = -numeric_limits<float>::max();
	}

	for (GIntervals::const_iterator iinterv = m_icur_interval; iinterv != m_intervals.end(); ++iinterv) {
		if (!iinterv->do_overlap(interval))
			break;

		for (vector<GenomeTrackArrays *>::iterator itrack = m_dependent_objs.begin(); itrack != m_dependent_objs.end(); ++itrack) {
			v = (*itrack)->get_sliced_val(iinterv - m_intervals.begin());
			if (!isnan(v)) {
				(*itrack)->m_last_sum += v;
				(*itrack)->m_last_min = min((*itrack)->m_last_min, v);
				(*itrack)->m_last_max = max((*itrack)->m_last_max, v);

				if ((*itrack)->m_functions[STDDEV])
					(*itrack)->m_mean_square_sum += v * v;

				if ((*itrack)->m_use_quantile)
					(*itrack)->m_sp.add(v);

				++(*itrack)->m_num_vs;
			}
		}
	}

	for (vector<GenomeTrackArrays *>::iterator itrack = m_dependent_objs.begin(); itrack != m_dependent_objs.end(); ++itrack) {
		if ((*itrack)->m_num_vs > 0)
			(*itrack)->m_last_avg = (*itrack)->m_last_nearest = (*itrack)->m_last_sum / (*itrack)->m_num_vs;
		else
			(*itrack)->m_last_avg = (*itrack)->m_last_nearest = (*itrack)->m_last_min = (*itrack)->m_last_max = (*itrack)->m_last_sum = numeric_limits<float>::quiet_NaN();

		// we are calaculating unbiased standard deviation:
		// sqrt(sum((x-mean)^2) / (N-1)) = sqrt(sum(x^2)/(N-1) - N*(mean^2)/(N-1))
		if ((*itrack)->m_functions[STDDEV])
			(*itrack)->m_last_stddev = (*itrack)->m_num_vs > 1 ?
				sqrt((*itrack)->m_mean_square_sum / ((*itrack)->m_num_vs - 1) - ((*itrack)->m_last_avg * (double)(*itrack)->m_last_avg) * ((*itrack)->m_num_vs / ((*itrack)->m_num_vs - 1))) :
				numeric_limits<float>::quiet_NaN();
	}
}

inline float GenomeTrackArrays::get_array_val(size_t islice)
{
	unsigned &hint = m_array_hints[islice];
	unsigned slice = m_slice[islice];

	// check the hint
	if (hint < m_array_vals.size() && m_array_vals[hint].idx == slice)
		return m_array_vals[hint].val;

	// check the previous hint + 1
	unsigned prev_hint;
	if (islice) {
		prev_hint = m_array_hints[islice - 1];
		hint = prev_hint + 1;
		if (hint < m_array_vals.size() && m_array_vals[hint].idx == slice)
			return m_array_vals[hint].val;
	} else
		prev_hint = 0;

	// run binary search
	ArrayVals::const_iterator iarray_val = lower_bound(m_array_vals.begin() + prev_hint, m_array_vals.end(), slice, LessIdx());
	hint = iarray_val - m_array_vals.begin();
	return iarray_val < m_array_vals.end() && iarray_val->idx == slice ? iarray_val->val : numeric_limits<float>::quiet_NaN();
}

#endif


/*
 * StreamPercentiler.h
 *
 *  Created on: Dec 8, 2011
 *      Author: hoichman
 */

#ifndef STREAMPERCENTILER_H_
#define STREAMPERCENTILER_H_

// StreamPercentiler returns the percentile of stream values.
// There is no real need to use StreamPercentiler if the number of values in the stream is small enough to fit the memory.
// Storing all the values in a vector, sorting them and then looking at index (=vector.size()*percentile) will have exactly the same effect.
//
// However things become trickier when the stream cannot fit into the physical memory.
//
// StreamPercentiler may calculate a percentile over the data that does not fit into the physical memory.
// This is done via random sampling of the data. StreamPercentiler also supports greater precision at the edges (i.e. close to 0 or to 1)
// where usually greater resolution is needed.
// The greater resolution at the edges is achieved by monitoring the N1 lowest and N2 highest values of the stream.
// During the construction of the object the user should supply the size of the buffer for random sampling as well as N1 and N2.

#include <math.h>
#include "StreamSampler.h"

template <class T>
class StreamPercentiler {
private:
	enum { LOWEST, HIGHEST, NUM_EXTREMES };

public:
	StreamPercentiler() { init(0, 0, 0, false); }
	StreamPercentiler(size_t rnd_sampling_buf_size, size_t lowest_vals_buf_size = 0, size_t highest_vals_buf_size = 0, bool do_reserve = false);

	void init(size_t rnd_sampling_buf_size, size_t lowest_vals_buf_size = 0, size_t highest_vals_buf_size = 0, bool do_reserve = false);
	void init_with_swap(size_t stream_size, vector<T> &samples, vector<T> &lowest_vals, vector<T> &highest_vals);

	void reset();

	size_t max_rnd_sampling_buf_size() const { return m_stream_sampler.max_reservoir_size(); }
	size_t cur_rnd_sampling_buf_size() const { return m_stream_sampler.cur_reservoir_size(); }
	size_t lowest_vals_buf_size() const { return m_extreme_vals_buf_size[LOWEST]; }
	size_t highest_vals_buf_size() const { return m_extreme_vals_buf_size[HIGHEST]; }

	const vector<T> &samples() const { return m_stream_sampler.samples(); }
	const vector<T> &lowest_vals() const { return m_extreme_vals[LOWEST]; }
	const vector<T> &highest_vals() const { return m_extreme_vals[HIGHEST]; }

	size_t add(const T &sample); // adds a sample, returns the number of samples inserted so far

	size_t stream_size() const { return m_stream_sampler.stream_size(); }
	const T get_percentile(double percentile, bool &is_estimated);

private:
	typedef bool (*Compare_t)(const T &, const T &);

	StreamSampler<T> m_stream_sampler;
	size_t           m_extreme_vals_buf_size[NUM_EXTREMES];
	Compare_t        m_compare_f[NUM_EXTREMES];
	vector<T>        m_extreme_vals[NUM_EXTREMES];
	bool             m_stream_sealed;
	bool             m_heaps_activated;

	static bool myless(const T &v1, const T &v2) { return v1 <= v2; }
	static bool mygreater(const T &v1, const T &v2) { return v1 >= v2; }
};


//------------------------------ IMPLEMENTATION ----------------------------------------

template <class T>
StreamPercentiler<T>::StreamPercentiler(size_t rnd_sampling_buf_size, size_t lowest_vals_buf_size, size_t highest_vals_buf_size, bool do_reserve)
{
	init(rnd_sampling_buf_size, lowest_vals_buf_size, highest_vals_buf_size, do_reserve);
}

template <class T>
void StreamPercentiler<T>::init(size_t rnd_sampling_buf_size, size_t lowest_vals_buf_size, size_t highest_vals_buf_size, bool do_reserve)
{
	m_stream_sampler.init(rnd_sampling_buf_size, do_reserve);
	m_extreme_vals_buf_size[LOWEST] = lowest_vals_buf_size;
	m_extreme_vals_buf_size[HIGHEST] = highest_vals_buf_size;
	m_compare_f[LOWEST] = myless;
	m_compare_f[HIGHEST] = mygreater;

	if (do_reserve) {
		for (int i = 0; i < NUM_EXTREMES; ++i) {
			if (m_extreme_vals_buf_size[i])
				m_extreme_vals[i].reserve(m_extreme_vals_buf_size[i] + 1);
		}
	}

	reset();
}

template <class T>
void StreamPercentiler<T>::init_with_swap(size_t stream_size, vector<T> &samples, vector<T> &lowest_vals, vector<T> &highest_vals)
{
	sort(samples.begin(), samples.end());
	sort(lowest_vals.begin(), lowest_vals.end());
	sort(highest_vals.begin(), highest_vals.end());

	m_stream_sampler.init_with_swap(stream_size, samples);
	m_extreme_vals_buf_size[LOWEST] = lowest_vals.size();
	m_extreme_vals_buf_size[HIGHEST] = highest_vals.size();
	m_heaps_activated = !lowest_vals.empty() || !highest_vals.empty();
	m_extreme_vals[LOWEST].swap(lowest_vals);
	m_extreme_vals[HIGHEST].swap(highest_vals);
	m_compare_f[LOWEST] = myless;
	m_compare_f[HIGHEST] = mygreater;
	m_stream_sealed = true;
}

template <class T>
void StreamPercentiler<T>::reset()
{
	m_stream_sampler.reset();
	m_extreme_vals[LOWEST].clear();
	m_extreme_vals[HIGHEST].clear();

	m_stream_sealed = false;
	m_heaps_activated = false;
}

template <class T>
size_t StreamPercentiler<T>::add(const T &sample)
{
	m_stream_sealed = false;

	// the stream reached its limit and we need to add another element =>
	// it's time to create the heaps
	if (m_stream_sampler.stream_size() == m_stream_sampler.max_reservoir_size()) {
		for (int i = 0; i < NUM_EXTREMES; ++i) {
			if (!m_extreme_vals_buf_size[i])
				continue;

			if (m_extreme_vals_buf_size[i] > m_stream_sampler.stream_size())
				m_extreme_vals[i] = m_stream_sampler.samples();
			else {
				vector<T> &samples = m_stream_sampler.samples();

				m_extreme_vals[i].reserve(m_extreme_vals_buf_size[i] + 1);
				m_extreme_vals[i].resize(m_extreme_vals_buf_size[i]);

				// sort the first m_extreme_vals_buf_size[i] elements, we can do it "in place" of stream_sampler buffer
				partial_sort(samples.begin(), samples.begin() + m_extreme_vals_buf_size[i], samples.end(), m_compare_f[i]);
				copy(samples.begin(), samples.begin() + m_extreme_vals_buf_size[i], m_extreme_vals[i].begin());
				make_heap(m_extreme_vals[i].begin(), m_extreme_vals[i].end(), m_compare_f[i]);
			}
			m_heaps_activated = true;
		}
	}

	size_t num_samples = m_stream_sampler.add(sample);

	if (m_heaps_activated) {
		for (int i = 0; i < NUM_EXTREMES; ++i) {
			if (m_extreme_vals[i].size() < m_extreme_vals_buf_size[i] || m_compare_f[i](sample, m_extreme_vals[i].front())) {
				m_extreme_vals[i].push_back(sample);

				if (m_extreme_vals[i].size() == m_extreme_vals_buf_size[i])
					make_heap(m_extreme_vals[i].begin(), m_extreme_vals[i].end(), m_compare_f[i]);
				else if (m_extreme_vals[i].size() == m_extreme_vals_buf_size[i] + 1) {
					push_heap(m_extreme_vals[i].begin(), m_extreme_vals[i].end(), m_compare_f[i]);
					pop_heap(m_extreme_vals[i].begin(), m_extreme_vals[i].end(), m_compare_f[i]);
					m_extreme_vals[i].pop_back();
				}
			}
		}
	}

	return num_samples;
}

template <class T>
const T StreamPercentiler<T>::get_percentile(double percentile, bool &is_estimated)
{
	percentile = max(percentile, 0.);
	percentile = min(percentile, 1.);

	if (!m_stream_sealed) {
		vector<T> &samples = m_stream_sampler.samples();
		sort(samples.begin(), samples.end());

		// if stream size is smaller than the reservoir size, don't use the heaps
		if (m_stream_sampler.stream_size() > m_stream_sampler.max_reservoir_size()) {
			for (int i = 0; i < NUM_EXTREMES; ++i)
				sort(m_extreme_vals[i].begin(), m_extreme_vals[i].end());
		}

		m_stream_sealed = true;
	}

	if (m_stream_sampler.stream_size() <= m_stream_sampler.max_reservoir_size()) {
		double index = (m_stream_sampler.stream_size() - 1) * percentile;
		size_t index1 = (size_t)floor(index);
		size_t index2 = (size_t)ceil(index);
		double weight = index - index1;
		is_estimated = false;
		return m_stream_sampler.samples()[index1] * (1 - weight) + m_stream_sampler.samples()[index2] * weight;
	}

	if (m_heaps_activated) {
		double index = (m_stream_sampler.stream_size() - 1) * percentile;
		size_t index1 = (size_t)floor(index);
		size_t index2 = (size_t)ceil(index);
		double weight = index - index1;
		T v1, v2;

		is_estimated = false;

		if (index1 < m_extreme_vals[LOWEST].size())
			v1 = m_extreme_vals[LOWEST][index1];
		else if (index1 >= m_stream_sampler.stream_size() - m_extreme_vals[HIGHEST].size())
			v1 = m_extreme_vals[HIGHEST][index1 - m_stream_sampler.stream_size() + m_extreme_vals[HIGHEST].size()];
		else {
			is_estimated = true;
			v1 = m_stream_sampler.samples()[(size_t)floor((m_stream_sampler.max_reservoir_size() - 1) * percentile)];
		}

		if (index2 < m_extreme_vals[LOWEST].size())
			v2 = m_extreme_vals[LOWEST][index2];
		else if (index2 >= m_stream_sampler.stream_size() - m_extreme_vals[HIGHEST].size())
			v2 = m_extreme_vals[HIGHEST][index2 - m_stream_sampler.stream_size() + m_extreme_vals[HIGHEST].size()];
		else {
			is_estimated = true;
			v2 = m_stream_sampler.samples()[(size_t)ceil((m_stream_sampler.max_reservoir_size() - 1) * percentile)];
		}

		return v1 * (1 - weight) + v2 * weight;
	}

	double index = (m_stream_sampler.max_reservoir_size() - 1) * percentile;
	size_t index1 = (size_t)floor(index);
	size_t index2 = (size_t)ceil(index);
	double weight = index - index1;
	is_estimated = true;
	return m_stream_sampler.samples()[index1] * (1 - weight) + m_stream_sampler.samples()[index2] * weight;
}

#endif /* STREAMPERCENTILER_H_ */

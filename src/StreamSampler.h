/*
 * StreamSampler.h
 *
 *  Created on: Nov 28, 2011
 *      Author: hoichman
 */

#ifndef STREAMSAMPLER_H_
#define STREAMSAMPLER_H_

#include <vector>

#include "TGLException.h"

using namespace std;

// StreamSampler picks up random N samples from a stream of an arbitrary size.
// The samples are guaranteed to be selected evenly over the whole stream.

template <class T>
class StreamSampler {
public:
	StreamSampler() { init(0, false); }
	StreamSampler(size_t reservoir_size, bool do_reserve = false) { init(reservoir_size, do_reserve); }

	void init(size_t reservoir_size, bool do_reserve = false);
	void init_with_swap(size_t stream_size, vector<T> &samples);
	void reset();

	const vector<T> &samples() const { return m_samples; }
	vector<T>       &samples() { return m_samples; }
	size_t          max_reservoir_size() const { return m_reservoir_size; }
	size_t          cur_reservoir_size() const { return m_samples.size(); }
	size_t          stream_size() const { return m_stream_size; }

	size_t          add(const T &sample); // returns the number of samples inserted so far

private:
	vector<T> m_samples;
	size_t  m_reservoir_size;
	size_t  m_stream_size;
};


//------------------------------ IMPLEMENTATION ----------------------------------------

template <class T>
void StreamSampler<T>::init(size_t reservoir_size, bool do_reserve)
{
	m_reservoir_size = reservoir_size;
	if (do_reserve)
		m_samples.reserve(m_reservoir_size);
	reset();
}

template <class T>
void StreamSampler<T>::init_with_swap(size_t stream_size, vector<T> &samples)
{
	m_stream_size = stream_size;
	m_reservoir_size = samples.size();
	m_samples.swap(samples);
}

template <class T>
void StreamSampler<T>::reset()
{
	m_samples.clear();
	m_stream_size = 0;
}

template <class T>
size_t StreamSampler<T>::add(const T &sample) {
	if (m_samples.size() < m_reservoir_size)
		m_samples.push_back(sample);
	else if (drand48() * (m_stream_size + 1) < m_reservoir_size)
		m_samples[(size_t)(drand48() * m_reservoir_size)] = sample;
	return (int64_t)++m_stream_size;
}

#endif /* STREAMSAMPLER_H_ */

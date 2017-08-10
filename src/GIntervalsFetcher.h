#ifndef _GINTERVALSFETCHER_H_INCLUDED_
#define _GINTERVALSFETCHER_H_INCLUDED_

template <typename Interval>
class GIntervalsFetcher {
public:
	enum Type { INTERVALS1D, BIGSET1D, INTERVALS2D, BIGSET2D };

	GIntervalsFetcher(Type type) : m_type(type) {}
	virtual ~GIntervalsFetcher() {}

	Type type() const { return m_type; }

	// if between two iterations the container changes, seal() must be called before the second begin_iter() or begin_chrom_iter(), i.e.:
	// begin_iter() / begin_chrom_iter()
	// ...
	// [container changes]
	// seal()
	// begin_iter() / begin_chrom_iter()
	// ...
	virtual void seal() = 0;

	// returns number of intervals in the intervals set
	virtual size_t size() const = 0;

	// Iterators: there are two modes of iterators
	// 1. If iterator is started with begin_iter():
	//    next()          - advances the iterator until end of all chromosomes is reached
	//    next_in_chrom() - does nothing and returns always false
	//    isend()         - returns true if end of all chromosomes is reached
	//    isend_chrom()   - always returns true
	//    cur_interval()  - returns a valid object as long as isend() returns false
	// 2. If iterator is started with begin_chrom_iter():
	//    next()          - advances the iterator until end of all chromosomes is reached
	//    next_in_chrom() - advances the iterator until end of current chromosome is reached
	//    isend()         - returns true if end of all chromosomes is reached
	//    isend_chrom()   - returns true if end of current chromosome is reached
	//    cur_interval()  - returns a valid object as long as isend_chrom() (and not isend()!) returns false
	virtual void begin_iter() = 0;

	// returns false if end reached (i.e. isend() returns true)
	// Note: if iterator is started by begin_chrom_iter next() continues iteration beyond the selected chromosome(s)
	virtual bool next() = 0;

	// returns false if end reached (i.e. isend_chrom() returns true)
	virtual bool next_in_chrom() = 0;

	virtual bool isend() const = 0;
	virtual bool isend_chrom() const = 0;

	// returns current iterator index
	virtual size_t iter_index() const = 0;

	// returns current iterator index within the current chromosome (or pair of chromosomes)
	virtual size_t iter_chrom_index() const = 0;

	// returns current interval
	virtual const Interval &cur_interval() const = 0;

protected:
	Type       m_type;
};

#endif


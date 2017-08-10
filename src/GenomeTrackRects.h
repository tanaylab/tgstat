/*
 * GenomeTrackRects.h
 *
 *  Created on: Jan 16, 2012
 *      Author: hoichman
 */

#ifndef GENOMETRACKRECTS_H_
#define GENOMETRACKRECTS_H_

#include <stdio.h>

#include "GenomeTrack2D.h"
#include "StatQuadTreeCached.h"
#include "StatQuadTreeCachedSerializer.h"

//---------------------------------------- GenomeTrackRects --------------------------------------

// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

template <class T>
class GenomeTrackRects : public GenomeTrack2D {
public:
	typedef StatQuadTreeCached<T, uint64_t>           QTree;
	typedef StatQuadTreeCachedSerializer<T, uint64_t> QTree_serializer;

	GenomeTrackRects(int64_t chunk_size, int64_t max_num_chunks) :
		GenomeTrack2D(typeid(T) == typeid(Rectangle_val<float>) ? RECTS : POINTS), m_qtree(chunk_size, max_num_chunks) { m_iqtree = NULL; }

	virtual ~GenomeTrackRects() { delete m_iqtree; }

	QTree &get_qtree() { return m_qtree; }
	void load();

	virtual void read_interval(const Rectangle &interval, const DiagonalBand &band);

	void write(const StatQuadTree<T, uint64_t> &qtree) { m_qtree.serialize(m_bfile, qtree); }

	// an alternative method to write() is through serializer. Call init_serializer() and then insert the objects to the serializer (serializer::begin() is called within init_serializer).
	// init_serializer() must be called anyway after init_write().
	void init_serializer(QTree_serializer &qtree_serializer, int64_t x1, int64_t y1, int64_t x2, int64_t y2, unsigned num_subtrees, bool check_overlaps, 
						 unsigned max_depth = 20, unsigned max_node_objs = 20);

	virtual bool begin_interval();
	virtual bool next_interval();
	virtual bool is_end_interval() { return !m_iqtree || m_iqtree->is_end(); }

protected:
	QTree                     m_qtree;
	typename QTree::Iterator *m_iqtree;
};

typedef GenomeTrackRects< Rectangle_val<float> > GenomeTrackRectsRects;
typedef GenomeTrackRects< Point_val<float> >     GenomeTrackRectsPoints;

//------------------------------------------- IMPLEMENTATION -------------------------------------------------------

template <class T>
void GenomeTrackRects<T>::read_interval(const Rectangle &interval, const DiagonalBand &band)
{
	if (m_bfile.opened()) {
		typename QTree::Stat result;

		load();
		if (band.is_non_empty_area())
			m_qtree.get_stat(interval, band, result);
		else
			m_qtree.get_stat(interval, result);
		m_last_occupied_area = result.occupied_area;
		m_last_weighted_sum = result.weighted_sum;
		m_last_min = result.min_val;
		m_last_max = result.max_val;
	} else {
		m_last_occupied_area = 0;
		m_last_weighted_sum = m_last_min = m_last_max = numeric_limits<double>::quiet_NaN();
	}
}

template <class T>
void GenomeTrackRects<T>::init_serializer(QTree_serializer &qtree_serializer, int64_t x1, int64_t y1, int64_t x2, int64_t y2, unsigned num_subtrees,
									      bool check_overlaps, unsigned max_depth, unsigned max_node_objs)
{
	qtree_serializer.begin(m_bfile, x1, y1, x2, y2, num_subtrees, check_overlaps, m_qtree.get_chunk_size(), m_qtree.get_max_num_chunks(), max_depth, max_node_objs);
}


template <class T>
bool GenomeTrackRects<T>::begin_interval()
{
	load();
	m_interval.chromid1() = m_chromid1;
	m_interval.chromid2() = m_chromid2;
	delete m_iqtree;
	m_iqtree = new typename QTree::Iterator(&m_qtree);
	bool retv = m_iqtree->begin();
	if (retv) {
		m_interval.start1() = (*m_iqtree)->get_x1();
		m_interval.end1() = (*m_iqtree)->get_x2();
		m_interval.start2() = (*m_iqtree)->get_y1();
		m_interval.end2() = (*m_iqtree)->get_y2();
	}
	return retv;
}

template <class T>
bool GenomeTrackRects<T>::next_interval()
{
	bool retv = m_iqtree->next();
	if (retv) {
		m_interval.start1() = (*m_iqtree)->get_x1();
		m_interval.end1() = (*m_iqtree)->get_x2();
		m_interval.start2() = (*m_iqtree)->get_y1();
		m_interval.end2() = (*m_iqtree)->get_y2();
	}
	return retv;
}

template <class T>
void GenomeTrackRects<T>::load()
{
	if (!m_loaded) {
		m_qtree.unserialize(m_bfile);
		m_loaded = true;
	}
}

#endif /* GENOMETRACKRECTS_H_ */

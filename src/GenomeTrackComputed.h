/*
 * GenomeTrackComputed.h
 *
 *  Created on: May 29, 2012
 *      Author: eitany
 */

#ifndef GENOMETRACKCOMPUTED_H_
#define GENOMETRACKCOMPUTED_H_

#include <stdio.h>

#include "GenomeTrack2D.h"
#include "StatQuadTreeCached.h"
#include "Computer2D.h"
#include "DiagonalBand.h"
#include "Rectangle.h"

//---------------------------------------- GenomeTrackComputed --------------------------------------

// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

//----------------------------------------- Computed_val -----------------------------------

template<typename T>
struct Computed_val : public Rectangle {
	T v;

	Computed_val() {}
	Computed_val(int64_t _x1, int64_t _y1, int64_t _x2, int64_t _y2, const T &_v = T()) : Rectangle(_x1, _y1, _x2, _y2), v(_v) {}
	Computed_val(const Rectangle &rect, const T &_v = T()) : Rectangle(rect), v(_v) {}

	double val(const Rectangle &rect, void* uptr) const {
        if (x1 == rect.x1 && x2 == rect.x2 && y1 == rect.y1 && y2 == rect.y2)
            return v;
        return ((Computer2D*)uptr)->compute(rect);
    }

	double val(const Rectangle &rect, const DiagonalBand &band, void* uptr) const {
        if (x1 == rect.x1 && x2 == rect.x2 && y1 == rect.y1 && y2 == rect.y2 && band.do_contain(rect))
            return v;
        return ((Computer2D*)uptr)->compute(rect, band);
    }

	char *debug_str() const {
		static char str[200];
		sprintf(str, "(%ld - %ld) (%ld - %ld) %g", x1, x2, y1, y2, (double)v);
		return str;
	}
};

typedef StatQuadTree<Computed_val<float>, uint64_t> ComputedQuadTree;
typedef StatQuadTreeCached<Computed_val<float>, uint64_t> ComputedQuadTreeCached;

class GenomeTrackComputed : public GenomeTrack2D {
public:
	typedef ComputedQuadTreeCached QTree;

	GenomeTrackComputed(const string &trackdb_path, int64_t chunk_size, int64_t max_num_chunks) :
		GenomeTrack2D(COMPUTED), m_qtree(chunk_size, max_num_chunks), m_iqtree(NULL), m_computer(NULL), m_trackdb_path(trackdb_path) {}

    virtual ~GenomeTrackComputed();

    ComputedQuadTreeCached &get_qtree() { return m_qtree; }

    // GenomeTrackComputed becomes the owner of the computer and destructs it at the end
	void set_computer(Computer2D* computer);

	Computer2D* get_computer() { return m_computer; }

	virtual void read_interval(const Rectangle &interval, const DiagonalBand &band);

	void load();

	void write(const ComputedQuadTree &qtree);

	virtual bool begin_interval();
	virtual bool next_interval();
	virtual bool is_end_interval() { return !m_iqtree || m_iqtree->is_end(); }

protected:
	ComputedQuadTreeCached            m_qtree;
	ComputedQuadTreeCached::Iterator *m_iqtree;

    // pointer to computer
    Computer2D* m_computer;

    // path to root of trackdb, needed to initialize computers
    string m_trackdb_path;
};

#endif /* GENOMETRACKCOMPUTED_H_ */

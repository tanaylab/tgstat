/*
 * GenomeTrackComputed.cpp
 *
 *  Created on: May 29, 2012
 *      Author: eitany
 */

#include "GenomeTrackComputed.h"

GenomeTrackComputed::~GenomeTrackComputed()
{
	delete m_iqtree;
	if (m_computer)
		delete m_computer;
}

void GenomeTrackComputed::set_computer(Computer2D* computer)
{
	delete m_computer;
	m_computer = computer;

    // set computer for quad tree
    m_qtree.set_uptr(m_computer);
}

void GenomeTrackComputed::read_interval(const Rectangle &interval, const DiagonalBand &band)
{
	if (m_bfile.opened()) {
		ComputedQuadTree::Stat result;

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


bool GenomeTrackComputed::begin_interval()
{
	load();
	m_interval.chromid1() = m_chromid1;
	m_interval.chromid2() = m_chromid2;
	delete m_iqtree;
	m_iqtree = new ComputedQuadTreeCached::Iterator(&m_qtree);
	bool retv = m_iqtree->begin();
	if (retv) {
        Rectangle &r = **m_iqtree;
        m_interval.start1() = r.x1;
        m_interval.end1() = r.x2;
        m_interval.start2() = r.y1;
        m_interval.end2() = r.y2;
	}
	return retv;
}

bool GenomeTrackComputed::next_interval()
{
	bool retv = m_iqtree->next();
	if (retv) {
        Rectangle &r = **m_iqtree;
        m_interval.start1() = r.x1;
        m_interval.end1() = r.x2;
        m_interval.start2() = r.y1;
        m_interval.end2() = r.y2;
	}
	return retv;
}

void GenomeTrackComputed::load()
{
	if (!m_loaded) {
        // read computer
        Computer2D* computer = Computer2D::unserializeComputer2D(m_bfile, m_trackdb_path.c_str(), m_chromid1, m_chromid2);
        set_computer(computer);

        // read tree
		m_qtree.unserialize(m_bfile);

		m_loaded = true;
	}
}

void GenomeTrackComputed::write(const ComputedQuadTree &qtree)
{
    if (m_computer == NULL)
        TGLError("cannot write track: m_computer not defined");

    // write computer
    Computer2D::serializeComputer2D(m_bfile, m_computer);

    // write tree
    m_qtree.serialize(m_bfile, qtree);
}

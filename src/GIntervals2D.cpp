#include "GIntervals2D.h"
#include "StatQuadTree.h"

//------------------------------------- GIntervals2D --------------------------------------------

void GIntervals2D::sort(bool (*cmp_function)(const GInterval2D &, const GInterval2D &))
{
	// do not sort automatically, check first that the intervals are not sorted yet
	for (iterator iinterv = begin() + 1; iinterv < end(); ++iinterv) {
		if (cmp_function(*iinterv, *(iinterv - 1))) {
			::sort(begin(), end(), cmp_function);
			break;
		}
	}
}

void GIntervals2D::verify_no_overlaps(const GenomeChromKey &chromkey, const char *error_prefix) const
{
	RectsQuadTree qtree;
	int chromid1 = -1;
	int chromid2 = -1;
	int chroms_start_idx = 0;

	for (const_iterator iinterv = begin(); iinterv != end(); ++iinterv) {
		if (iinterv != begin() && *iinterv < *(iinterv - 1))
			TGLError<GIntervalsFetcher2D>(UNSORTED_INTERVALS, "%sTo verify overlaps 2D intervals must be sorted", error_prefix);

		if (chromid1 != iinterv->chromid1() || chromid2 != iinterv->chromid2()) {
			chromid1 = iinterv->chromid1();
			chromid2 = iinterv->chromid2();
			qtree.reset(0, 0, chromkey.get_chrom_size(chromid1), chromkey.get_chrom_size(chromid2));
			chroms_start_idx = iinterv - begin();
		}

		if (qtree.do_intersect(*iinterv)) {
			vector<Rectangle> intersection;
			vector<uint64_t> intersected_objs_indices;

			qtree.intersect(*iinterv, intersection, intersected_objs_indices);
			const GInterval2D &interv = at(chroms_start_idx + intersected_objs_indices[0]);

			TGLError<GIntervalsFetcher2D>(OVERLAPPING_INTERVAL,
					"%sIntervals (%s, %ld, %ld, %s, %ld, %ld) and (%s, %ld, %ld, %s, %ld, %ld) overlap",
					error_prefix,
					chromkey.id2chrom(iinterv->chromid1()).c_str(), iinterv->start1(), iinterv->end1(),
					chromkey.id2chrom(iinterv->chromid2()).c_str(), iinterv->start2(), iinterv->end2(),
					chromkey.id2chrom(interv.chromid1()).c_str(), interv.start1(), interv.end1(),
					chromkey.id2chrom(interv.chromid2()).c_str(), interv.start2(), interv.end2());
		}

		qtree.insert(RectsQuadTree::ValueType(*iinterv));
	}
}

double GIntervals2D::surface() const
{
	double surface = 0;

	for (const_iterator iinterv = begin(); iinterv < end(); ++iinterv)
		surface += iinterv->surface();
	return surface;
}

double GIntervals2D::surface(int chromid1, int chromid2) const
{
	double surface = 0;

	for (const_iterator iinterv = begin(); iinterv < end(); ++iinterv) {
		if (iinterv->chromid1() == chromid1 && iinterv->chromid2() == chromid2) 
			surface += iinterv->surface();
	}
	return surface;
}

void GIntervals2D::begin_chrom_iter(int chromid1, int chromid2)
{
	build_chrom_map();
	m_cur_chromid1 = chromid1;
	m_cur_chromid2 = chromid2;
	m_iter_chrom_index = 0;

	if (chromid1 >= m_num_chroms || chromid2 >= m_num_chroms && chromid1 >= m_num_chroms - 1) 
		m_iinterval = end();
	else if (chromid2 >= m_num_chroms) 
		m_iinterval = m_chrom2itr[chroms2idx(chromid1 + 1, 0)];
	else
		m_iinterval = m_chrom2itr[chroms2idx(chromid1, chromid2)];
}

GIntervals2D::const_iterator GIntervals2D::get_chrom_begin() const
{
	build_chrom_map();
	return m_iinterval->chromid1() < m_num_chroms && m_iinterval->chromid2() < m_num_chroms ?
		m_chrom2itr[chroms2idx(m_iinterval->chromid1(), m_iinterval->chromid2())] :
		end();
}

GIntervals2D::const_iterator GIntervals2D::get_chrom_end() const
{
	build_chrom_map();
	if (m_iinterval->chromid1() < m_num_chroms && m_iinterval->chromid2() < m_num_chroms) {
		int idx = chroms2idx(m_iinterval->chromid1(), m_iinterval->chromid2()) + 1;
		return idx < m_chrom2itr.size() ? m_chrom2itr[idx] : end();
	}
	return end();
}

#include "GIntervals.h"

void GIntervals::sort(bool (*cmp_function)(const GInterval &, const GInterval &))
{
	// do not sort automatically, check first that the intervals are not sorted yet
	for (iterator iinterv = begin() + 1; iinterv < end(); ++iinterv) {
		if (cmp_function(*iinterv, *(iinterv - 1))) {
			::sort(begin(), end(), cmp_function);
			break;
		}
	}
}

// intervs are expected to be already sorted
void GIntervals::verify_no_overlaps(const GenomeChromKey &chromkey, const char *error_prefix) const
{
	for (const_iterator iinterv = begin() + 1; iinterv < end(); ++iinterv) {
		if (*iinterv < *(iinterv - 1))
			TGLError<GIntervalsFetcher1D>(UNSORTED_INTERVALS, "%sTo verify overlaps intervals must be sorted", error_prefix);

		if (iinterv->chromid == (iinterv - 1)->chromid && ((iinterv - 1)->end > iinterv->start))
			TGLError<GIntervalsFetcher1D>(OVERLAPPING_INTERVAL, "%sIntervals (%s, %ld, %ld) and (%s, %ld, %ld) overlap",
										  error_prefix,
										  chromkey.id2chrom((iinterv - 1)->chromid).c_str(), (iinterv - 1)->start, (iinterv - 1)->end,
										  chromkey.id2chrom(iinterv->chromid).c_str(), iinterv->start, iinterv->end);
	}
}

int64_t GIntervals::range() const
{
	int64_t range = 0;

	for (const_iterator iinterv = begin(); iinterv < end(); ++iinterv)
		range += iinterv->range();
	return range;
}

int64_t GIntervals::range(int chromid) const
{
	int64_t range = 0;

	for (const_iterator iinterv = begin(); iinterv < end(); ++iinterv) {
		if (iinterv->chromid == chromid) 
			range += iinterv->range();
	}
	return range;
}

// intervs are expected to be already sorted
void GIntervals::unify_overlaps(bool unify_touching_intervals)
{
	if (empty())
		return;

	size_t cur_idx = 0;

	for (size_t i = 1; i < size(); i++) {
		if (operator[](cur_idx).chromid != operator[](i).chromid || operator[](cur_idx).end < operator[](i).start || !unify_touching_intervals && operator[](cur_idx).end == operator[](i).start)
			operator[](++cur_idx) = operator[](i);
		// unite overlapping intervals
		else if (operator[](cur_idx).end < operator[](i).end)
			operator[](cur_idx).end = operator[](i).end;
	}
	erase(begin() + cur_idx + 1, end());
}

void GIntervals::unify(const GIntervals &intervs1, const GIntervals &intervs2, GIntervals &res_intervs)
{
	const_iterator iintervs[] = { intervs1.begin(), intervs2.begin() };
	const_iterator intervends[] = { intervs1.end(), intervs2.end() };
	int last_chromid[] = { -1, -1 };
	int idx = 0;

	res_intervs.clear();
	res_intervs.reserve(intervs1.size() + intervs2.size());

	// merge the two intervs into one vector sorted by chrom and start coord (complexity: O(n))
	while (iintervs[0] != intervends[0] && iintervs[1] != intervends[1]) {
		if (iintervs[0]->chromid == iintervs[1]->chromid)
			idx = iintervs[0]->start < iintervs[1]->start ? 0 : 1;
		else if (last_chromid[0] != iintervs[0]->chromid || last_chromid[1] != iintervs[1]->chromid) {
			idx = compare_by_start_coord(*iintervs[0], *iintervs[1]) ? 0 : 1;
			last_chromid[0] = iintervs[0]->chromid;
			last_chromid[1] = iintervs[1]->chromid;
		}

		res_intervs.push_back(*iintervs[idx]);
		++iintervs[idx];
	}

	for (int i = 0; i < 2; i++) {
		for (const_iterator iinterv = iintervs[i]; iinterv != intervends[i]; ++iinterv)
			res_intervs.push_back(*iinterv);
	}

	// having the whole set of merged (i.e. sorted) intervals let's unify the overlapping intervals (complexity: O(n))
	res_intervs.unify_overlaps();
}

void GIntervals::intersect(const GIntervals &_intervs1, const GIntervals &_intervs2, GIntervals &res_intervs)
{
	GIntervals intervs[] = { _intervs1, _intervs2 };  // yes, we want to copy the intervals (we are going to change them)
	iterator iintervs[2] = { intervs[0].begin(), intervs[1].begin() };
	int last_chromid[2] = { -1, -1 };
	int idx = 0;

	res_intervs.clear();
	while (iintervs[0] != intervs[0].end() && iintervs[1] != intervs[1].end()) {
		if (iintervs[0]->chromid == iintervs[1]->chromid) {
			if (iintervs[0]->start < iintervs[1]->start && iintervs[0]->end <= iintervs[1]->start)
				++iintervs[0];
			else if (iintervs[1]->start < iintervs[0]->start && iintervs[1]->end <= iintervs[0]->start)
				++iintervs[1];
			else { // intervals intersect
				int64_t start = max(iintervs[0]->start, iintervs[1]->start);
				int64_t end = min(iintervs[0]->end, iintervs[1]->end);

				res_intervs.push_back(GInterval(iintervs[0]->chromid, start, end, 0));
				for (int i = 0; i < 2; i++) {
					if (iintervs[i]->end == end)
						++iintervs[i];
					else
						iintervs[i]->start = end;
				}
			}
		} else {
			if (last_chromid[0] != iintervs[0]->chromid || last_chromid[1] != iintervs[1]->chromid) {
				idx = compare_by_start_coord(*iintervs[0], *iintervs[1]) ? 0 : 1;
				last_chromid[0] = iintervs[0]->chromid;
				last_chromid[1] = iintervs[1]->chromid;
			} else
				++iintervs[idx];
		}
	}
}

void GIntervals::diff(const GIntervals &_intervs1, const GIntervals &_intervs2, GIntervals &res_intervs)
{
	GIntervals intervs[] = { _intervs1, _intervs2 };  // yes, we want to copy the intervals (we are going to change them)
	iterator iintervs[] = { intervs[0].begin(), intervs[1].begin() };
	int last_chromid[] = { -1, -1 };
	int idx = 0;

	res_intervs.clear();
	while (iintervs[0] != intervs[0].end() && iintervs[1] != intervs[1].end()) {
		if (iintervs[0]->chromid == iintervs[1]->chromid) {
			if (iintervs[0]->start < iintervs[1]->start && iintervs[0]->end <= iintervs[1]->start) {
				res_intervs.push_back(*iintervs[0]);
				++iintervs[0];
			}
			else if (iintervs[1]->start < iintervs[0]->start && iintervs[1]->end <= iintervs[0]->start)
				++iintervs[1];
			else { // intervals intersect
				int64_t intersect_start = max(iintervs[0]->start, iintervs[1]->start);
				int64_t intersect_end = min(iintervs[0]->end, iintervs[1]->end);

				if (iintervs[0]->start < intersect_start)
					res_intervs.push_back(GInterval(iintervs[0]->chromid, iintervs[0]->start, intersect_start, 0));

				for (int i = 0; i < 2; i++) {
					if (iintervs[i]->end == intersect_end)
						++iintervs[i];
					else
						iintervs[i]->start = intersect_end;
				}
			}
		} else {
			if (last_chromid[0] != iintervs[0]->chromid || last_chromid[1] != iintervs[1]->chromid) {
				idx = compare_by_start_coord(*iintervs[0], *iintervs[1]) ? 0 : 1;
				last_chromid[0] = iintervs[0]->chromid;
				last_chromid[1] = iintervs[1]->chromid;
			} else {
				if (!idx)
					res_intervs.push_back(*iintervs[0]);
				++iintervs[idx];
			}
		}
	}

	for (const_iterator iinterv = iintervs[0]; iinterv != intervs[0].end(); ++iinterv)
		res_intervs.push_back(*iinterv);
}

const GInterval *GIntervals::containing_interval(const GInterval &interv)
{
	// run binary search
	GIntervals::const_iterator istart_interval = begin();
	GIntervals::const_iterator iend_interval = end();

	while (iend_interval - istart_interval > 1) {
		GIntervals::const_iterator imid_interval = istart_interval + (iend_interval - istart_interval) / 2;

		if (imid_interval->do_overlap(interv))
			return imid_interval->do_contain(interv) ? &*imid_interval : NULL;

		// is mid_interval < interval?
		if (GIntervals::compare_by_start_coord(*imid_interval, interv))
			istart_interval = imid_interval;
		else
			iend_interval = imid_interval;
	}

	if (iend_interval - istart_interval == 1 && istart_interval->do_overlap(interv))
		return istart_interval->do_contain(interv) ? &*istart_interval : NULL;

	return NULL;
}

void GIntervals::read(const GenomeChromKey &chromkey, istream &tab, int nostrand)
{
	string chrom;
	int64_t start, end;
	int strand;

	strand = 1;
	tab >> chrom;
	while(tab) {
		tab >> start >> end;
		if(!nostrand) {
			tab >> strand;
		}
		GInterval interval(chromkey.chrom2id(chrom.c_str()), start, end, (char)strand);
		interval.verify(chromkey);
		push_back(interval);
		tab >> chrom;
	}
}

void GIntervals::read_bed(const GenomeChromKey &chromkey, istream &bed)
{
	string chrom;
	int64_t start, end;
	char strandcode;
	float score;
	string name;
	int strand = 0;

	bed >> chrom;
	while(bed) {
		bed >> start >> end >> name >> score >> strandcode;
		try {
			strand = GInterval::char2strand(strandcode);
		} catch (TGLException &e) {
			TGLError<GInterval>(GInterval::BAD_STRAND, "Reading interval (%s, %ld, %ld, %c): %s", chrom.c_str(), start, end, strandcode, e.msg());
		}
		GInterval interval(chromkey.chrom2id(chrom.c_str()), start, end, (char)strand);
		interval.verify(chromkey);
		push_back(interval);
		while(bed.get() != '\n') {}
		bed >> chrom;
	}
}

void GIntervals::begin_chrom_iter(int chromid)
{
	build_chrom_map();
	m_cur_chromid = chromid;
	m_iter_chrom_index = 0;
	if (chromid < m_chrom2itr.size())
		m_iinterval = m_chrom2itr[chromid];
	else
		m_iinterval = end();
}

GIntervals::const_iterator GIntervals::get_chrom_begin() const
{
	build_chrom_map();
	return m_chrom2itr[m_iinterval->chromid];
}

GIntervals::const_iterator GIntervals::get_chrom_end() const
{
	build_chrom_map();
	return m_iinterval->chromid + 1 < m_chrom2itr.size() ? m_chrom2itr[m_iinterval->chromid + 1] : end();
}

void GIntervals::write(const GenomeChromKey &chromkey, ostream &tab)
{
	for(const_iterator i = begin(); i != end(); i++) {
		tab << chromkey.id2chrom(i->chromid)
		<< "\t" << i->start
		<< "\t" << i->end
		<< "\t" << (int)i->strand << "\n";
	}
}

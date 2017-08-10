#include <errno.h>
#include <limits>

#include "GenomeArraysCsv.h"
#include "GInterval.h"
#include "strutil.h"

void GenomeArraysCsv::init(const char *filename, const GenomeChromKey &chromkey)
{
	m_bfile.close();
	m_chroms_positions.clear();
	m_colnames.clear();
	m_intervals.clear();

	m_chromkey = (GenomeChromKey *)&chromkey;

	if (m_bfile.open(filename, "r"))
		TGLError<GenomeArraysCsv>(FILE_ERROR, "Opening a file %s: %s", filename, strerror(errno));

	long lineno = 0;

	lineno += split_line(m_bfile, m_fields, '\t');

	if (m_fields.size() < GInterval::NUM_COLS + 1)
		TGLError<GenomeArraysCsv>(FORMAT_ERROR, "File %s, line %ld: invalid format", filename, lineno);

	for (int i = 0; i < GInterval::NUM_COLS; ++i) {
		if (m_fields[i] != GInterval::COL_NAMES[i]) 
			TGLError<GenomeArraysCsv>(FORMAT_ERROR, "File %s, line %ld: invalid format", filename, lineno);
	}

	for (vector<string>::const_iterator ifield = m_fields.begin() + GInterval::NUM_COLS; ifield < m_fields.end(); ++ifield) 
		m_colnames.push_back(*ifield);

	m_chroms_positions.resize(m_chromkey->get_num_chroms());

	while (1) {
		Position pos(m_bfile.tell(), lineno);

		lineno += read_fields(pos);
		if (m_bfile.eof()) 
			return;

		try {
		    int chromid = m_chromkey->chrom2id(m_fields[GInterval::CHROM]);
			m_chroms_positions[chromid].push_back(pos);
		} catch (TGLException &) {
			// there might be unrecognized chromosomes, ignore them
		}
	}
}

const GIntervals &GenomeArraysCsv::get_intervals(int chromid)
{
	m_intervals.clear();
	const Positions &positions = m_chroms_positions[chromid];

	for (Positions::const_iterator ipos = positions.begin(); ipos != positions.end(); ++ipos) {
		long lineno = ipos->lineno + read_fields(*ipos);
		char *endptr;
		int64_t start, end;

		start = strtoll(m_fields[GInterval::START].c_str(), &endptr, 10);
		if (*endptr || start < 0) 
			TGLError<GenomeArraysCsv>(FORMAT_ERROR, "File %s, line %ld: invalid format of start coordinate", m_bfile.file_name().c_str(), lineno);

		end = strtoll(m_fields[GInterval::END].c_str(), &endptr, 10);
		if (*endptr) 
			TGLError<GenomeArraysCsv>(FORMAT_ERROR, "File %s, line %ld: invalid format of start coordinate", m_bfile.file_name().c_str(), lineno);

		if (start >= end) 
			TGLError<GenomeArraysCsv>(FORMAT_ERROR, "File %s, line %ld: start coordinate exceeds or equals the end coordinate", m_bfile.file_name().c_str(), lineno);

		if (end > m_chromkey->get_chrom_size(chromid)) 
			TGLError<GenomeArraysCsv>(FORMAT_ERROR, "File %s, line %ld: end coordinate exceeds chromosome's size", m_bfile.file_name().c_str(), lineno);

		m_intervals.push_back(GInterval(chromid, start, end, 0, (void *)&*ipos));
	}

	m_intervals.sort();

	for (GIntervals::const_iterator iinterv = m_intervals.begin() + 1; iinterv < m_intervals.end(); ++iinterv) {
		if (((iinterv - 1)->end > iinterv->start))
			TGLError<GenomeArraysCsv>(FORMAT_ERROR, "File %s, lines %ld and %ld: intervals overlap",
								 m_bfile.file_name().c_str(), ((Position *)(iinterv - 1)->udata)->lineno + 1, ((Position *)iinterv->udata)->lineno + 1);
	}

	return m_intervals;
}

void GenomeArraysCsv::get_sliced_vals(GIntervals::const_iterator iinterval, vector<float> &vals)
{
	const Position *position = (Position *)iinterval->udata;
	long lineno = position->lineno + read_fields(*position);
	char *endptr;

	vals.clear();
	for (vector<string>::const_iterator ifield = m_fields.begin() + GInterval::NUM_COLS; ifield < m_fields.end(); ++ifield) {
		if (ifield->empty()) 
			vals.push_back(numeric_limits<float>::quiet_NaN());
		else {
			vals.push_back(strtod(ifield->c_str(), &endptr));
			if (*endptr) 
				TGLError<GenomeArraysCsv>(FORMAT_ERROR, "File %s, line %ld, column %ld: invalid value", m_bfile.file_name().c_str(), lineno, ifield - m_fields.begin());
		}
	}
}

int GenomeArraysCsv::read_fields(const Position &pos)
{
	m_bfile.seek(pos.bytes, SEEK_SET);
	int total_numcols = m_colnames.size() + GInterval::NUM_COLS;
	int delta_lines = split_line(m_bfile, m_fields, '\t', total_numcols);

	if (m_fields.empty()) {
		if (m_bfile.error())
			TGLError<GenomeArraysCsv>(FILE_ERROR, "Reading a file %s: %s", m_bfile.file_name().c_str(), strerror(errno));
		return 0;
	}

	if (m_fields.size() != total_numcols) 
		TGLError<GenomeArraysCsv>(FORMAT_ERROR, "File %s, line %ld: expecting %ld columns, read %ld",
								  m_bfile.file_name().c_str(), pos.lineno + delta_lines, total_numcols, m_fields.size());

	return delta_lines;
}

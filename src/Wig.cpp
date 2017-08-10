#include <errno.h>

#include "port.h"
#include "strutil.h"
#include "Wig.h"

string Wig::FIELDS_STRS[NUM_FIELDS] = { "chrom=", "start=", "step=", "span=" };

#include <fstream>

void Wig::init(const GenomeChromKey &chromkey, const string &filename, bool ignore_unknown_chroms)
{
	m_chromkey = (GenomeChromKey *)&chromkey;
	m_ignore_unknown_chroms = ignore_unknown_chroms;
	m_chrom_fpos.resize(chromkey.get_num_chroms(), -1);
	m_chrom_lineno.resize(chromkey.get_num_chroms(), -1);

	// read the file to determine the positions of first chrom appearences
	int last_chromid = -1;

	if (m_bfile.open(filename.c_str(), "r"))
		TGLError<Wig>(FILE_ERROR, "Failed to open WIG file %s: %s", m_bfile.file_name().c_str(), strerror(errno));

	Rec rec;
	int64_t lineno = 1;
	int64_t last_lineno;
	long last_fpos;
	int last_rec_type = -1;
	int chrom_rec_type = -1;

	while (1) {
		last_fpos = m_bfile.tell();
		last_lineno = lineno;

		if (!read_record(rec, lineno))
			break;

		if (last_rec_type != -1 && rec.type == BEDGRAPH_REC && last_rec_type != BEDGRAPH_REC)
			TGLError<Wig>(INVALID_FORMAT, "Invalid format of WIG file %s: looks like a bizar mixture of WIG / BedGraph formats", m_bfile.file_name().c_str());

		if (rec.type == FIXED_STEP_REC || rec.type == VAR_STEP_REC || rec.type == BEDGRAPH_REC) {
			int chromid = -1;

			if (rec.type == FIXED_STEP_REC)
				chromid = rec.fixed_step.chromid;
			else if (rec.type == VAR_STEP_REC)
				chromid = rec.var_step.chromid;
			else
				chromid = rec.bedgraph.chromid;

			chrom_rec_type = rec.type;

			if (chromid == -1) // unknown chrom
				continue;

			if (chromid != last_chromid) {
				last_chromid = chromid;
				if (m_chrom_fpos[chromid] != -1)
					TGLError<Wig>(INVALID_FORMAT, "Invalid format of WIG file %s: file is not sorted by chromosomes", m_bfile.file_name().c_str());
				m_chrom_fpos[chromid] = last_fpos;
				m_chrom_lineno[chromid] = last_lineno;
			}
		} else if (rec.type == VAL_REC) {
			if (chrom_rec_type != FIXED_STEP_REC)
				TGLError<Wig>(INVALID_FORMAT, "Invalid format of WIG file %s, line %ld: value is not preceeded by fixedStep header",
						m_bfile.file_name().c_str(), lineno - 1);
		} else if (rec.type == COORD_VAL_REC) {
			if (chrom_rec_type != VAR_STEP_REC)
				TGLError<Wig>(INVALID_FORMAT, "Invalid format of WIG file %s, line %ld: value is not preceeded by variableStep header",
						m_bfile.file_name().c_str(), lineno - 1);
		}
	}
}

bool Wig::get_data(int chromid, GIntervals &intervals)
{
	intervals.clear();

	if (m_chrom_fpos[chromid] == -1)
		return false;

	int64_t lineno = m_chrom_lineno[chromid];
	int64_t span = -1;
	int64_t start = -1;
	int64_t step = 0;
	Rec rec;

	m_bfile.seek(m_chrom_fpos[chromid], SEEK_SET);

	while (1) {
		if (!read_record(rec, lineno))
			break;

		if (rec.type == FIXED_STEP_REC) {
			if (rec.fixed_step.chromid != chromid)
				break;

			start = rec.fixed_step.start;
			step = rec.fixed_step.step;
			span = rec.fixed_step.span;

			if (start >= (int64_t)m_chromkey->get_chrom_size(rec.fixed_step.chromid))
				TGLError<Wig>(BAD_START, "WIG file %s, line %ld: start coordinate %ld exceeds chromosome %s boundaries",
						m_bfile.file_name().c_str(), lineno - 1, start, m_chromkey->id2chrom(chromid).c_str());

			if (span > step)
				TGLError<Wig>(INVALID_FORMAT, "Invalid format of WIG file %s, line %ld: span exceeds step", m_bfile.file_name().c_str(), lineno - 1);
		} else if (rec.type == VAR_STEP_REC) {
			if (rec.var_step.chromid != chromid)
				break;

			start = -1;
			step = -1;
			span = rec.var_step.span;
		} else if (rec.type == VAL_REC) {
			if (start < 0)
				TGLError<Wig>(INVALID_FORMAT, "Invalid format of WIG file %s, line %ld: value is not preceeded by fixedStep header",
						m_bfile.file_name().c_str(), lineno - 1);

			if (start >= (int64_t)m_chromkey->get_chrom_size(chromid))
				TGLError<Wig>(INVALID_FORMAT, "[1] Invalid format of WIG file %s, line %ld: coordinate %ld exceeds chromosome %s boundaries",
						m_bfile.file_name().c_str(), lineno - 1, start, m_chromkey->id2chrom(chromid).c_str());

			int64_t end_coord = min(start + span, (int64_t)m_chromkey->get_chrom_size(chromid));
			intervals.push_back(GInterval(chromid, start, end_coord, 0, GInterval::cast2udata(rec.val)));
			start += step;
		} else if (rec.type == COORD_VAL_REC) {
			if (span < 0)
				TGLError<Wig>(INVALID_FORMAT, "Invalid format of WIG file %s, line %ld: value is not preceeded by variableStep header",
						m_bfile.file_name().c_str(), lineno - 1);

			if (rec.coord_val.coord >= (int64_t)m_chromkey->get_chrom_size(chromid))
				TGLError<Wig>(BAD_COORD, "[2] WIG file %s, line %ld: coordinate %ld exceeds chromosome %s boundaries",
						m_bfile.file_name().c_str(), lineno - 1, rec.coord_val.coord, m_chromkey->id2chrom(chromid).c_str());

			int64_t end_coord = min(rec.coord_val.coord + span, (int64_t)m_chromkey->get_chrom_size(chromid));
			intervals.push_back(GInterval(chromid, rec.coord_val.coord, end_coord, 0, GInterval::cast2udata(rec.coord_val.val)));
		} else if (rec.type == BEDGRAPH_REC) {
			if (rec.bedgraph.chromid != chromid)
				break;

			if (rec.bedgraph.start >= rec.bedgraph.end)
				TGLError<Wig>(BAD_COORD, "WIG file %s, line %ld: start coordinate exceeds end coordinate", m_bfile.file_name().c_str());

			if (rec.bedgraph.end > (int64_t)m_chromkey->get_chrom_size(chromid))
				TGLError<Wig>(BAD_COORD, "WIG file %s, line %ld: coordinate %ld exceeds chromosome %s boundaries",
						m_bfile.file_name().c_str(), lineno - 1, rec.bedgraph.end, m_chromkey->id2chrom(chromid).c_str());

			intervals.push_back(GInterval(chromid, rec.bedgraph.start, rec.bedgraph.end, 0, GInterval::cast2udata(rec.bedgraph.val)));
		}
	}

	for (GIntervals::const_iterator iinterv = intervals.begin() + 1; iinterv < intervals.end(); ++iinterv) {
		if ((iinterv - 1)->end > iinterv->start)
			TGLError<Wig>(INVALID_FORMAT, "Invalid format of WIG file %s: file contains overlapping intervals", m_bfile.file_name().c_str());
	}

	return true;
}

bool Wig::read_record(Rec &rec, int64_t &lineno)
{
	vector<string> fields;

	while (1) {
		lineno += split_line_by_space_chars(m_bfile, fields);

		if (m_bfile.error())
			TGLError<Wig>(FILE_ERROR, "Failed to read WIG file %s: %s", m_bfile.file_name().c_str(), strerror(errno));

		if (fields.empty())
			break;

		if (fields[0] == "track" || fields[0].size() && fields[0][0] == '#')
			continue;

		if (fields[0] == "variableStep") {
			if (fields.size() < 2 || fields.size() > 3 ||
					fields[1].compare(0, FIELDS_STRS[CHROM_FIELD].size(), FIELDS_STRS[CHROM_FIELD]) ||
					(fields.size() == 3 && fields[2].compare(0, FIELDS_STRS[SPAN_FIELD].size(), FIELDS_STRS[SPAN_FIELD])))
				TGLError<Wig>(INVALID_FORMAT, "Invalid format of WIG file %s, line %ld", m_bfile.file_name().c_str(), lineno - 1);

			rec.var_step.chromid = str2chromid(fields[1].c_str() + FIELDS_STRS[CHROM_FIELD].length(), lineno - 1);
			rec.var_step.span = fields.size() == 3 ? str2span(fields[2].c_str() + FIELDS_STRS[SPAN_FIELD].length(), lineno - 1) : 1;
			rec.type = VAR_STEP_REC;
			return true;
		}

		if (fields[0] == "fixedStep") {
			if (fields.size() < 4 || fields.size() > 5 ||
					fields[1].compare(0, FIELDS_STRS[CHROM_FIELD].size(), FIELDS_STRS[CHROM_FIELD]) ||
					fields[2].compare(0, FIELDS_STRS[START_FIELD].size(), FIELDS_STRS[START_FIELD]) ||
					fields[3].compare(0, FIELDS_STRS[STEP_FIELD].size(), FIELDS_STRS[STEP_FIELD]) ||
					(fields.size() == 5 && fields[4].compare(0, FIELDS_STRS[SPAN_FIELD].size(), FIELDS_STRS[SPAN_FIELD])))
				TGLError<Wig>(INVALID_FORMAT, "Invalid format of WIG file %s, line %ld", m_bfile.file_name().c_str(), lineno - 1);

			rec.fixed_step.chromid = str2chromid(fields[1].c_str() + FIELDS_STRS[CHROM_FIELD].length(), lineno - 1);
			rec.fixed_step.start = str2start(fields[2].c_str() + FIELDS_STRS[START_FIELD].length(), lineno - 1);
			rec.fixed_step.step = str2step(fields[3].c_str() + FIELDS_STRS[STEP_FIELD].length(), lineno - 1);
			rec.fixed_step.span = fields.size() == 5 ? str2span(fields[4].c_str() + FIELDS_STRS[SPAN_FIELD].length(), lineno - 1) : 1;
			rec.type = FIXED_STEP_REC;
			return true;
		}

		if (fields.size() == 1) {
			rec.val = str2val(fields[0].c_str(), lineno - 1);
			rec.type = VAL_REC;
			return true;
		}

		if (fields.size() == 2) {
			rec.coord_val.coord = str2coord(fields[0].c_str(), 1, lineno - 1);
			rec.coord_val.val = str2val(fields[1].c_str(), lineno - 1);
			rec.type = COORD_VAL_REC;
			return true;
		}

		if (fields.size() == 4 && fields[0].size() > 3 && !fields[0].compare(0, 3, "chr")) {
			rec.bedgraph.chromid = str2chromid(fields[0].c_str(), lineno - 1);
			rec.bedgraph.start = str2coord(fields[1].c_str(), 0, lineno - 1);
			rec.bedgraph.end = str2coord(fields[2].c_str(), 0, lineno - 1);
			rec.bedgraph.val = str2val(fields[3].c_str(), lineno - 1);
			rec.type = BEDGRAPH_REC;
			return true;
		}

		TGLError<Wig>(INVALID_FORMAT, "Invalid format of WIG file %s, line %ld", m_bfile.file_name().c_str(), lineno - 1);
	}

	return false;
}

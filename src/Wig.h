#ifndef WIG_H_
#define WIG_H_

#include <vector>
#include "BufferedFile.h"
#include "GenomeChromKey.h"
#include "GIntervals.h"

using namespace std;

//----------------------------------------------------------------------------
// Wig class parses file in WIG or BedGraph formats and returns its data
//
// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

class Wig {
public:
	enum Errors { FILE_ERROR, INVALID_FORMAT, BAD_CHROM, BAD_START, BAD_STEP, BAD_SPAN, BAD_VAL, BAD_COORD };

	Wig() {}
	Wig(const GenomeChromKey &chromkey, const string &filename, bool ignore_unknown_chroms) { init(chromkey, filename, ignore_unknown_chroms); }

	void init(const GenomeChromKey &chromkey, const string &filename, bool ignore_unknown_chroms);

	// returns true if data for chromosome exists; (float)udata of each interval contains the value
	bool get_data(int chromid, GIntervals &intervals);

private:
	enum { CHROM_FIELD, START_FIELD, STEP_FIELD, SPAN_FIELD, NUM_FIELDS };
	enum RecordType { FIXED_STEP_REC, VAR_STEP_REC, VAL_REC, COORD_VAL_REC, BEDGRAPH_REC };

	struct FixedStepRec {
		int     chromid;
		int64_t start;
		int64_t step;
		int64_t span;
	};

	struct VarStepRec {
		int     chromid;
		int64_t span;
	};

	struct CoordValRec {
		int64_t coord;
		float   val;
	};

	struct BedGraph {
		int     chromid;
		int64_t start;
		int64_t end;
		float   val;
	};

	struct Rec {
		RecordType type;
		union {
			FixedStepRec fixed_step;
			VarStepRec   var_step;
			CoordValRec  coord_val;
			float        val;
			BedGraph     bedgraph;
		};
	};

	static string FIELDS_STRS[NUM_FIELDS];

	GenomeChromKey *m_chromkey;
	vector<float>   m_chrom_data;
	vector<long>    m_chrom_fpos;
	vector<int64_t> m_chrom_lineno;
	BufferedFile    m_bfile;
	bool            m_ignore_unknown_chroms;

	bool read_record(Rec &rec, int64_t &lineno);

	int     str2chromid(const char *str, int64_t lineno);
	int64_t str2span(const char *str, int64_t lineno);
	int64_t str2step(const char *str, int64_t lineno);
	int64_t str2start(const char *str, int64_t lineno);
	int64_t str2coord(const char *str, int64_t base, int64_t lineno);
	float   str2val(const char *str, int64_t lineno);
};


//-------------------------------------IMPLEMENTATION -----------------------------------------

inline int Wig::str2chromid(const char *str, int64_t lineno)
{
	try {
		return m_chromkey->chrom2id(str);
	} catch (TGLException &e) {
		if (m_ignore_unknown_chroms)
			return -1;
		TGLError<Wig>(BAD_CHROM, "WIG file %s, line %ld: %s\n", m_bfile.file_name().c_str(), lineno, e.msg());
	}
	return -1;
}

inline int64_t Wig::str2span(const char *str, int64_t lineno)
{
	char *endptr;
	int64_t span = strtoll(str, &endptr, 10);
	if (*endptr || span < 1)
		TGLError<Wig>(BAD_SPAN, "WIG file %s, line %ld: invalid value of span", m_bfile.file_name().c_str(), lineno);
	return span;
}

inline int64_t Wig::str2step(const char *str, int64_t lineno)
{
	char *endptr;
	int64_t step = strtoll(str, &endptr, 10);
	if (*endptr || step < 1)
		TGLError<Wig>(BAD_STEP, "WIG file %s, line %ld: invalid value of step", m_bfile.file_name().c_str(), lineno);
	return step;
}

inline int64_t Wig::str2start(const char *str, int64_t lineno)
{
	char *endptr;
	int64_t start = strtoll(str, &endptr, 10);
	if (*endptr || start < 1)
		TGLError<Wig>(BAD_START, "WIG file %s, line %ld: invalid value of start", m_bfile.file_name().c_str(), lineno);
	return start - 1;
}

inline int64_t Wig::str2coord(const char *str, int64_t base, int64_t lineno)
{
	char *endptr;
	int64_t coord = strtoll(str, &endptr, 10);
	if (*endptr || coord < base)
		TGLError<Wig>(BAD_COORD, "WIG file %s, line %ld: invalid coordinate", m_bfile.file_name().c_str(), lineno);
	return coord - base;
}

inline float Wig::str2val(const char *str, int64_t lineno)
{
	char *endptr;
	float val = strtod(str, &endptr);
	if (*endptr)
		TGLError<Wig>(BAD_VAL, "WIG file %s, line %ld: invalid value", m_bfile.file_name().c_str(), lineno);
	return val;
}

#endif /* WIG_H_ */

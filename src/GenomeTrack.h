/*
 * GenomeTrack.h
 *
 *  Created on: Mar 10, 2010
 *      Author: hoichman
 */

#ifndef GENOMETRACK_H_
#define GENOMETRACK_H_

#include <map>
#include <string>

#include "BufferedFile.h"
#include "GenomeChromKey.h"

using namespace std;

// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

class GenomeTrack {
public:
	typedef map<string, string> TrackAttrs;

	enum Type { FIXED_BIN, SPARSE, ARRAYS, RECTS, POINTS, COMPUTED, OLD_RECTS1, OLD_RECTS2, OLD_COMPUTED1, OLD_COMPUTED2, OLD_COMPUTED3, NUM_TYPES };
	enum Error { NOT_2D, BAD_FORMAT, OBSOLETE_FORMAT, MISMATCH_FORMAT, FILE_ERROR, BAD_ATTRS };

	static const char *TYPE_NAMES[NUM_TYPES];
	static const int   FORMAT_SIGNATURES[NUM_TYPES];

	virtual ~GenomeTrack() {}

	Type get_type() const { return m_type; }

	bool is_1d() const { return is_1d(m_type); }
	bool is_2d() const { return is_2d(m_type); }

	static bool is_1d(Type type) { return IS_1D_TRACK[type]; }
	static bool is_2d(Type type) { return !IS_1D_TRACK[type]; }

	// returns true if the track file is opened
	bool opened() const { return m_bfile.opened(); }

	const string &file_name() const { return m_bfile.file_name(); }

	static Type get_type(const char *track_dir, const GenomeChromKey &chromkey, bool return_obsolete_types = false);

	static void load_attrs(const char *track, const char *filename, TrackAttrs &attrs);

	static void save_attrs(const char *track, const char *filename, const TrackAttrs &attrs);

	static const string &get_1d_filename(const GenomeChromKey &chromkey, int chromid) { return chromkey.id2chrom(chromid); }

	static const string get_2d_filename(const GenomeChromKey &chromkey, int chromid1, int chromid2) {
		return chromkey.id2chrom(chromid1) + "-" + chromkey.id2chrom(chromid2);
	}

	static const int get_chromid_1d(const GenomeChromKey &chromkey, const string &filename) { return chromkey.chrom2id(filename); }

	static const pair<int, int> get_chromid_2d(const GenomeChromKey &chromkey, const string &filename);

protected:
	static const bool IS_1D_TRACK[NUM_TYPES];

	BufferedFile m_bfile;
	Type         m_type;

	GenomeTrack(Type type) : m_type(type) {}

	void read_type(const char *filename, const char *mode = "rb");

	void write_type(const char *filename, const char *mode = "wb");

	static Type s_read_type(const char *filename, const char *mode = "rb");

	static Type s_read_type(BufferedFile &bfile, const char *filename, const char *mode = "rb");
};

#endif /* GENOMETRACK_H_ */

/*
 * GenomeTrack.cpp
 *
 *  Created on: Mar 10, 2010
 *      Author: hoichman
 */

#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "GenomeTrack.h"
#include "TGLException.h"

const char *GenomeTrack::TYPE_NAMES[GenomeTrack::NUM_TYPES] = { "dense", "sparse", "array", "rectangles", "points", "computed", "obsolete rectangles", "obsolete rectangles", "obsolete computed", "obsolete computed", "obsolete computed" };
const bool  GenomeTrack::IS_1D_TRACK[NUM_TYPES] =             { true,    true,     true,     false,        false,    false,      false,                 false,                 false,               false,               false };
const int   GenomeTrack::FORMAT_SIGNATURES[NUM_TYPES] =       { 0,       -1,       -8,       -9,           -10,      -11,        -4,                    -6,                    -3,                  -5,                  -7 };

const pair<int, int> GenomeTrack::get_chromid_2d(const GenomeChromKey &chromkey, const string &filename)
{
	size_t pos = filename.find_first_of("-");

	if (pos == string::npos)
		TGLError<GenomeTrack>(NOT_2D, "File %s does not belong to 2D track", filename.c_str());

	string chrom1(filename, 0, pos);
	string chrom2(filename, pos + 1);

	return pair<int, int>(chromkey.chrom2id(chrom1), chromkey.chrom2id(chrom2));
}

GenomeTrack::Type GenomeTrack::get_type(const char *track_dir, const GenomeChromKey &chromkey, bool return_obsolete_types)
{
	if (access(track_dir, F_OK))
		TGLError<GenomeTrack>(FILE_ERROR, "Accessing directory %s: %s\n", track_dir, strerror(errno));

	// first try to access it as a 1D track
	Type type;
	bool is_1d = false;

	for (size_t chromid = 0; chromid < chromkey.get_num_chroms(); chromid++) {
		try {
			type = s_read_type((string(track_dir) + "/" + get_1d_filename(chromkey, chromid)).c_str());
			is_1d = true;
			break;
		} catch (TGLException &err) {}
	}

	if (is_1d) {
		if (type != FIXED_BIN && type != SPARSE && type != ARRAYS)
			TGLError<GenomeTrack>(BAD_FORMAT, "Invalid format of track file at %s", track_dir);
		return type;
	}

	for (size_t chromid1 = 0; chromid1 < chromkey.get_num_chroms(); chromid1++) {
		for (size_t chromid2 = 0; chromid2 < chromkey.get_num_chroms(); chromid2++) {
			bool is_2d = false;

			try {
				type = s_read_type((string(track_dir) + "/" + get_2d_filename(chromkey, chromid1, chromid2)).c_str());
				is_2d = true;
			} catch (TGLException &) {}

			if (is_2d) {
				if (type == OLD_RECTS1 || type == OLD_RECTS2 ||type == OLD_COMPUTED1 || type == OLD_COMPUTED2 || type == OLD_COMPUTED3) {
					if (return_obsolete_types)
						return type;
					TGLError<GenomeTrack>(OBSOLETE_FORMAT, "Track file at %s is in obsolete format and requires conversion", track_dir);
				}

				if (type != RECTS && type != POINTS && type != COMPUTED)
					TGLError<GenomeTrack>(BAD_FORMAT, "Invalid format of track file at %s", track_dir);
				return type;
			}
		}
	}

	TGLError<GenomeTrack>(BAD_FORMAT, "Invalid format of track at %s", track_dir);
	return NUM_TYPES;
}

void GenomeTrack::read_type(const char *filename, const char *mode)
{
	Type type = s_read_type(m_bfile, filename, mode);

	if (type != m_type)
		TGLError<GenomeTrack>(MISMATCH_FORMAT, "Track file %s is in %s format while expected to be in %s format", filename, TYPE_NAMES[type], TYPE_NAMES[m_type]);
}

GenomeTrack::Type GenomeTrack::s_read_type(const char *filename, const char *mode)
{
	BufferedFile bfile;
	return s_read_type(bfile, filename, mode);
}

GenomeTrack::Type GenomeTrack::s_read_type(BufferedFile &bfile, const char *filename, const char *mode)
{
	if (bfile.open(filename, mode))
		TGLError<GenomeTrack>(FILE_ERROR, "Opening a track file %s: %s", filename, strerror(errno));

	int format_signature;

	if (bfile.read(&format_signature, sizeof(format_signature)) != sizeof(format_signature)) {
		if (bfile.error())
			TGLError<GenomeTrack>(FILE_ERROR, "Reading a track file %s: %s", filename, strerror(errno));
		TGLError<GenomeTrack>(BAD_FORMAT, "Invalid format of track file %s", filename);
	}

	if (format_signature > 0)
		return FIXED_BIN;

	for (int type = FIXED_BIN + 1; type < NUM_TYPES; ++type) {
		if (format_signature == FORMAT_SIGNATURES[type])
			return (Type)type;
	}

	TGLError<GenomeTrack>(BAD_FORMAT, "Invalid format of genome track file %s", filename);
	return NUM_TYPES;
}

void GenomeTrack::write_type(const char *filename, const char *mode)
{
	umask(07);

	if (m_bfile.open(filename, mode))
		TGLError<GenomeTrack>(FILE_ERROR, "Opening a track file %s: %s", filename, strerror(errno));

	if (m_bfile.write(&FORMAT_SIGNATURES[m_type], sizeof(FORMAT_SIGNATURES[m_type])) != sizeof(FORMAT_SIGNATURES[m_type])) {
		if (m_bfile.error())
			TGLError<GenomeTrack>(FILE_ERROR, "Failed to write a %s track file %s: %s", TYPE_NAMES[m_type], filename, strerror(errno));
		TGLError<GenomeTrack>(FILE_ERROR, "Failed to write a %s track file %s", TYPE_NAMES[m_type], filename);
	}
}

void GenomeTrack::load_attrs(const char *, const char *filename, TrackAttrs &attrs)
{
	BufferedFile bfile;
	int c;
	int idx = 0;
	string name;
	string val;

	attrs.clear();

	if (bfile.open(filename, "rb")) {
		if (errno == ENOENT)   // no file = no attributes
			return; 
		TGLError<GenomeTrack>(FILE_ERROR, "Failed to read attributes file %s: %s", filename, strerror(errno));
	}

	while ((c = bfile.getc()) >= 0) {
		if (c) {
			if (idx) 
				val.push_back((char)c);
			else
				name.push_back((char)c);
		} else {
			if (idx) {
				if (name.empty() || val.empty())
					TGLError<GenomeTrack>(BAD_FORMAT, "Invalid format of attributes file %s", filename); 

				if (attrs.find(name) != attrs.end()) // duplicated attributes
					TGLError<GenomeTrack>(BAD_FORMAT, "Invalid format of attributes file %s", filename); 

				attrs[name] = val;
				name.clear();
				val.clear();
			}
			idx = 1 - idx;
		}
	}

	if (bfile.error()) 
		TGLError<GenomeTrack>(FILE_ERROR, "Failed to read attributes file %s: %s", filename, strerror(errno));

	if (idx) 
		TGLError<GenomeTrack>(BAD_FORMAT, "Invalid format of attributes file %s", filename); 
}

void GenomeTrack::save_attrs(const char *track, const char *filename, const TrackAttrs &attrs)
{
	bool empty_attrs = true;

	for (TrackAttrs::const_iterator iattr = attrs.begin(); iattr != attrs.end(); ++iattr) { 
		if (!iattr->second.empty()) {
			empty_attrs = false;
			break;
		}
	}

	if (empty_attrs) {
		if (unlink(filename) && errno != ENOENT)
			TGLError<GenomeTrack>(FILE_ERROR, "Failed accessing attributes file %s: %s", filename, strerror(errno));
		return;
	}

	for (TrackAttrs::const_iterator iattr = attrs.begin(); iattr != attrs.end(); ++iattr) {
		if (iattr->first.empty())
				TGLError<GenomeTrack>(BAD_ATTRS, "Track %s: attribute name is an empty string", track); 
	}

	BufferedFile bfile;

	if (bfile.open(filename, "wb"))
		TGLError<GenomeTrack>(FILE_ERROR, "Failed to write attributes file %s: %s", filename, strerror(errno));

	for (TrackAttrs::const_iterator iattr = attrs.begin(); iattr != attrs.end(); ++iattr) {
		if (!iattr->second.empty())  {
			bfile.write(iattr->first.c_str(), iattr->first.length() + 1);
			bfile.write(iattr->second.c_str(), iattr->second.length() + 1);
		}
	}

	if (bfile.error())
		TGLError<GenomeTrack>(FILE_ERROR, "Failed to write attributes file %s: %s", filename, strerror(errno));
}


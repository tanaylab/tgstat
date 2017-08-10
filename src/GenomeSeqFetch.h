/*
 * GenomeSeqFetch.h
 *
 *  Created on: Nov 21, 2010
 *      Author: hoichman
 */

#ifndef GENOMESEQFETCH_H_
#define GENOMESEQFETCH_H_

#include <string>
#include <vector>

#include "BufferedFile.h"
#include "GenomeChromKey.h"
#include "GenomeUtils.h"
#include "GInterval.h"

// -------------------- GenomeSeqFetch  -----------------------
// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

class GenomeSeqFetch {
public:
	enum Errors { FILE_READ_FAILED };

	GenomeSeqFetch() : m_cur_chromid(-1) {}

	void set_seqdir(const std::string &dir) { m_seqdir = dir; }
	void read_interval(const GInterval &interval, const GenomeChromKey &chromkey, std::vector<char> &result);

private:
	std::string  m_seqdir;
	int          m_cur_chromid;
	BufferedFile m_bfile;
};

#endif /* GENOMESEQFETCH_H_ */

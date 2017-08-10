#include <errno.h>

#include "port.h"
BASE_CC_FILE

#include "GenomeSeqFetch.h"

void GenomeSeqFetch::read_interval(const GInterval &interval, const GenomeChromKey &chromkey, vector<char> &result)
{
	if (m_cur_chromid != interval.chromid) {
		char filename[PATH_MAX];

		m_cur_chromid = interval.chromid;
		sprintf(filename, "%s/%s.seq", m_seqdir.c_str(), chromkey.id2chrom(interval.chromid).c_str());
		m_bfile.close();
		m_bfile.open(filename, "rb");

		if (m_bfile.error())
			TGLError<GenomeSeqFetch>(FILE_READ_FAILED, "Reading sequence file %s failed: %s", filename, strerror(errno));
	}

	interval.verify(chromkey, false);
	result.clear();

	int64_t size = min(interval.end, m_bfile.file_size()) - interval.start;

	if (size < 0)
		return;

	if (!size)
		size++;

	result.resize(size);
	m_bfile.seek(interval.start, SEEK_SET);
	if (m_bfile.read(&*result.begin(), result.size()) != result.size()) {
		if (m_bfile.error())
			TGLError<GenomeSeqFetch>(FILE_READ_FAILED, "Reading sequence file %s failed: %s", m_bfile.file_name().c_str(), strerror(errno));

		TGLError<GenomeSeqFetch>(FILE_READ_FAILED, "Reading sequence file %s failed", m_bfile.file_name().c_str());
	}

	// if strand == -1 make the sequence reverse-complementary
	if (interval.strand == -1) {
		for (vector<char>::iterator i = result.begin(); i != result.end(); ++i)
			*i = basepair2complementary(*i);
		reverse(result.begin(), result.end());
	}
}

#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "TGLException.h"
#include "GenomeTrack2D.h"

void GenomeTrack2D::init_read(const char *filename, int chromid1, int chromid2)
{
	m_bfile.close();
	m_loaded = false;

	if (!access(filename, R_OK) || errno != ENOENT)
		read_type(filename);

	m_chromid1 = chromid1;
	m_chromid2 = chromid2;
}

void GenomeTrack2D::init_write(const char *filename, int chromid1, int chromid2)
{
	m_bfile.close();
	m_loaded = false;
	write_type(filename);
	m_chromid1 = chromid1;
	m_chromid2 = chromid2;
}

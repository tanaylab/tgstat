#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "BufferedFile.h"
#include "TGLException.h"

int64_t BufferedFile::file_size(const char *path)
{
	struct stat st;

	if (stat(path, &st))
		TGLError("Cannot stat file %s: %s\n", path, strerror(errno));
	return (int64_t)st.st_size;
}


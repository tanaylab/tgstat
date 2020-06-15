#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/sendfile.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "FileUtils.h"
#include "TGLException.h"

struct FD {
    int fd{-1};
    ~FD() {
        if (fd != -1)
            close(fd);
    }
};

void FileUtils::copy_file(const char *src, const char *tgt)
{
    FD srcfd;
    FD tgtfd;
    struct stat srcstat;

    if ((srcfd.fd = open(src, O_RDONLY, 0)) == -1)
        TGLError(errno, "Error opening file %s for reading: %s", src, strerror(errno));
    if (fstat(srcfd.fd, &srcstat) == -1)
        TGLError(errno, "Error trying to stat file %s: %s", src, strerror(errno));
    if ((tgtfd.fd = creat(tgt, srcstat.st_mode)) == -1)
        TGLError(errno, "Error opeining file %s for writing: %s", tgt, strerror(errno));
    if (sendfile(tgtfd.fd, srcfd.fd, NULL, srcstat.st_size) == -1)
        TGLError(errno, "Error copying file %s to %s: %s\n", src, tgt, strerror(errno));
}

void FileUtils::move_file(const char *src, const char *tgt)
{
    if (rename(src, tgt) == -1) {
        if (errno == EXDEV) {
            FileUtils::copy_file(src, tgt);
            if (unlink(src) == -1) {
                auto olderrno = errno;
                unlink(tgt);
                TGLError(olderrno, "Error removing file %s: %s", src, strerror(olderrno));
            }
        } else
            TGLError(errno, "Error moving file %s to %s: %s\n", src, tgt);
    }
}

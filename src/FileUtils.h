#ifndef FILEUTILS_H_INCLUDED
#define FILEUTILS_H_INCLUDED

namespace FileUtils {
    // Makes a fast copy of a file while preserving permissions.
    // Throws TGLException on error, TGLException::code contains errno.
    void copy_file(const char *src, const char *tgt);

    // Renames a file or copies it, if the target is located in a different file system.
    // Throws TGLException on error, TGLException::code contains errno.
    void move_file(const char *src, const char *tgt);
}

#endif

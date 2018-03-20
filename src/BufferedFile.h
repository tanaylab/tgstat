#ifndef BUFFEREDFILE_H_
#define BUFFEREDFILE_H_

#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <string>
#include <vector>

#ifndef _FILE_OFFSET_BITS
#define _FILE_OFFSET_BITS 64
#endif

class BufferedFile {
public:
	BufferedFile() { init(1024); }
	BufferedFile(unsigned bufsize) { init(1024); }
	~BufferedFile();

	// returns 0 on success, -1 on failure
	int open(const char *path, const char *mode);

	// see fclose for return value
	int close();

	// see fgetc for return value
	int getc();

	// see fread for return value
	size_t read(void *ptr, size_t size);

	// see fwrite for return value
	size_t write(const void *ptr, size_t size);

	// see ftell for return value
	int64_t tell() const { return m_virt_pos; }

	// see fseek for return value
	int seek(int64_t offset, int whence);

	int error() const { return !m_fp || ferror(m_fp); }

	int eof() const { return m_eof; }

	bool opened() const { return m_fp != NULL; }

	const std::string &file_name() const { return m_filename; }

	int64_t file_size() const { return m_file_size; }

	// in case of error throws an exception
	static int64_t file_size(const char *path);

protected:
	FILE       *m_fp{NULL};
	int         m_eof;
	std::string m_filename;
	char       *m_buf{NULL};
	unsigned    m_bufsize;
	int64_t     m_file_size;
	int64_t     m_virt_pos;
	int64_t     m_phys_pos;
	int64_t     m_sbuf_pos;
	int64_t     m_ebuf_pos;

	void init(unsigned bufsize);

private:
	// restrict copying of class (otherwise need reference counters for m_buf and m_fp)
	BufferedFile(const BufferedFile &) = delete;
	BufferedFile &operator=(const BufferedFile &obj) = delete;
};


class BufferedFiles : public std::vector<BufferedFile *> {
public:
	BufferedFiles() : std::vector<BufferedFile *>() {}
	BufferedFiles(size_type n) : std::vector<BufferedFile *>(n, NULL) {}

	~BufferedFiles() {
		for (iterator ifile = begin(); ifile != end(); ++ifile)
			delete *ifile;
	}

private:
	BufferedFiles(const BufferedFiles &) = delete;
	BufferedFiles &operator=(const BufferedFiles &) = delete;
};

//--------------------------- IMPLEMENTATION ------------------------------------

inline BufferedFile::~BufferedFile()
{
	close();
	delete []m_buf;
}

inline void BufferedFile::init(unsigned bufsize) {
	m_fp = NULL;
	m_eof = true;
	m_buf = NULL;
	m_bufsize = bufsize;
	m_file_size = 0;
	m_virt_pos = -1;
	m_phys_pos = 0;
	m_sbuf_pos = 0;
	m_ebuf_pos = 0;
	m_buf = new char[m_bufsize];
}

inline int BufferedFile::open(const char *path, const char *mode)
{
	close();
	m_filename = path;
	m_fp = fopen(path, mode);
	if (m_fp) {
		m_eof = false;
		m_virt_pos = m_phys_pos = 0;
		m_sbuf_pos = m_ebuf_pos = 0;

		fseeko(m_fp, 0, SEEK_END);
		m_file_size = ftello(m_fp);
		fseeko(m_fp, 0, SEEK_SET);
		return 0;
	}
	return -1;
}

inline int BufferedFile::close()
{
	if (m_fp) {
		int retv = fclose(m_fp);
		m_fp = NULL;
		m_eof = true;
		m_phys_pos = -1;
		return retv;
	}
	return 0;
}

inline int BufferedFile::getc()
{
	// is the new read already cached?
	if (m_virt_pos >= m_sbuf_pos && m_virt_pos + 1 <= m_ebuf_pos)
		return m_buf[m_virt_pos++ - m_sbuf_pos];

	char c;
	if (!read(&c, 1))
		return -1;
	return c;
}

inline size_t BufferedFile::read(void *ptr, size_t size)
{
	// is the new read already cached?
	if (m_virt_pos >= m_sbuf_pos && m_virt_pos + (long)size <= m_ebuf_pos) {
		memcpy(ptr, m_buf + m_virt_pos - m_sbuf_pos, size);
		m_virt_pos += size;
		return size;
	}

	// we must perform a new read
	if (m_phys_pos != m_virt_pos)
		fseeko(m_fp, m_virt_pos, SEEK_SET);

	// We could check whether the buffer overlaps the reading section,
	// but that would make the code slower and more complicated and would not
	// give any real benefits considering that the buffer size is relatively small.
	// In some cases (when the buffer overlaps the central part of the reading section
	// that would be even harmful as it would result in two freads instead of one.

	// if size is smaller than bufsize => cache the read
	if (size <= m_bufsize) {
		size_t bytes_read = fread(m_buf, 1, m_bufsize, m_fp);

		m_phys_pos = m_virt_pos + bytes_read;
		m_sbuf_pos = m_virt_pos;
		m_ebuf_pos = m_phys_pos;

		if (bytes_read > size)
			bytes_read = size;

		m_virt_pos += bytes_read;
		memcpy(ptr, m_buf, bytes_read);
		if (!bytes_read && feof(m_fp))
			m_eof = true;
		return bytes_read;
	}

	// make the read without caching
	size_t bytes_read = fread(ptr, 1, size, m_fp);
	m_virt_pos += bytes_read;
	m_phys_pos = m_virt_pos;
	if (!bytes_read && feof(m_fp))
		m_eof = true;
	return bytes_read;
}

inline size_t BufferedFile::write(const void *ptr, size_t size)
{
	if (m_phys_pos != m_virt_pos) {
		fseeko(m_fp, m_virt_pos, SEEK_SET);
		m_phys_pos = m_virt_pos;
	}

	size_t retv = fwrite(ptr, 1, size, m_fp);

	if (retv) {
		// if the cached buffer overlaps the written region => trash the buffer
		if (std::max(m_virt_pos, m_sbuf_pos) < std::min(m_virt_pos + (long)retv, m_ebuf_pos))
			m_sbuf_pos = m_ebuf_pos = 0;

		m_virt_pos += retv;
		m_phys_pos = m_virt_pos;

		if (m_file_size < m_virt_pos)
			m_file_size = m_virt_pos;
	}

	return retv;
}

inline int BufferedFile::seek(int64_t offset, int whence)
{
	switch (whence) {
	case SEEK_END:
		if (offset < 0 || m_file_size - 1 < offset) {
			errno = EINVAL;
			return -1;
		}
		m_virt_pos = m_file_size - offset - 1;
		break;
	case SEEK_CUR:
		if (m_virt_pos + offset < 0 || m_virt_pos + offset > m_file_size - 1) {
			errno = EINVAL;
			return -1;
		}
		m_virt_pos += offset;
		break;
	case SEEK_SET:
		if (offset < 0 || offset > m_file_size) {
			errno = EINVAL;
			return -1;
		}
		m_virt_pos = offset;
		break;
	default:
		errno = EINVAL;
		return -1;
	}

	m_eof = m_virt_pos == m_file_size;
	return 0;
}

#endif

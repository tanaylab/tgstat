#ifndef HASHFUNC_H_
#define HASHFUNC_H_

#ifdef __GXX_EXPERIMENTAL_CXX0X__
#include <functional>
#else
#include <ext/hash_fun.h>
#endif

#include <sys/param.h>
#ifdef _BSD
#include <sys/endian.h>
#elif defined(__linux__)
#include <byteswap.h>
#endif

namespace __gnu_cxx
{
	template<> struct hash<std::string>
	{
		size_t operator()(const std::string &s) const { return hash<const char *>()(s.c_str()); }
	};

	// This hash works for a pair of ids that can reach up to sizeof_t.
	// Ideally we would like the hash to be:
	//      v1 | bit_reverse(v2)
	// Yet bit_reverse requires quite some computational effort. So our hash is:
	//      little_edian(v1) | big_endian(v2)
	template<> struct hash< std::pair<size_t, size_t> >
	{
		size_t operator()(const std::pair<size_t, size_t> &v) const {
#if (__WORDSIZE == 64)
			return v.first ^ bswap_64(v.second);
#else
			return v.first ^ bswap_32(v.second);
#endif
		}
	};
}

#endif /* HASHFUNC_H_ */

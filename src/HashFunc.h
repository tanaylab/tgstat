#ifndef HASHFUNC_H_
#define HASHFUNC_H_

#include <functional>

#include <sys/param.h>
#ifdef _BSD
#include <sys/endian.h>
#elif defined(__linux__)
#include <byteswap.h>
#elif defined(__APPLE__)
#include <libkern/OSByteOrder.h>
#define bswap_64 OSSwapInt64
#define bswap_32 OSSwapInt32
#endif

namespace std
{
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


    // This hash works for a pair of ids that can reach up to sizeof_t.
    // Ideally we would like the hash to be:
    //      v1 | bit_reverse(v2)
    // Yet bit_reverse requires quite some computational effort. So our hash is:
    //      little_edian(v1) | big_endian(v2)
    template<> struct hash< std::pair<unsigned, unsigned> >
    {
        size_t operator()(const std::pair<size_t, size_t> &v) const {
#if (__WORDSIZE == 64)
            return v.first | (v.second << 32);
#else
            return v.first ^ bswap_32(v.second);
#endif
        }
    };
}

#endif /* HASHFUNC_H_ */

#ifndef HASHFUNC_H_
#define HASHFUNC_H_

#include <cstdint>
#include <functional>

#ifndef BSWAP_8
#define	BSWAP_8(x)	((x) & 0xff)
#endif

#ifndef BSWAP_16
#define	BSWAP_16(x)	((BSWAP_8(x) << 8) | BSWAP_8((x) >> 8))
#endif

#ifndef BSWAP_32
#define	BSWAP_32(x)	((BSWAP_16(x) << 16) | BSWAP_16((x) >> 16))
#endif

#ifndef BSWAP_64
#define	BSWAP_64(x)	((BSWAP_32(x) << 32) | BSWAP_32((x) >> 32))
#endif

namespace std
{
	// This hash works for a pair of ids that can reach up to sizeof_t.
	// Ideally we would like the hash to be:
	//      v1 | bit_reverse(v2)
	// Yet bit_reverse requires quite some computational effort. So our hash is:
	//      little_edian(v1) | big_endian(v2)
	template<> struct hash< std::pair<uint64_t, uint64_t> >
	{
		uint64_t operator()(const std::pair<uint64_t, uint64_t> &v) const {
#if (__WORDSIZE == 64)
			return v.first ^ BSWAP_64(v.second);
#else
			return v.first ^ BSWAP_32(v.second);
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
            return v.first ^ BSWAP_32(v.second);
#endif
        }
    };
}

#endif /* HASHFUNC_H_ */

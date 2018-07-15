#ifndef RANDOM_SHUFFLE_H_INCLUDED
#define RANDOM_SHUFFLE_H_INCLUDED

#include <algorithm>
#include <stdlib.h>

// tgs_random_shuffle replaces std::random_shuffle.
// 
// std::random_shuffle is implemented differently on Linux and OSX which results in different shuffling.
// The results achieved on Linux are not reproducible on OSX even if the same seed is used.

template<class RandomIt>
void tgs_random_shuffle(RandomIt first, RandomIt last, double (*rnd_func)())
{
    typename std::iterator_traits<RandomIt>::difference_type i, n;

    n = last - first;
    for (i = n - 1; i > 0; --i) {
        using std::swap;
        swap(first[i], first[(size_t)(rnd_func() * (i + 1))]);
    }
}


#endif

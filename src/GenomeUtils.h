#ifndef GENOMEUTILS_H_
#define GENOMEUTILS_H_

#include <string>

extern char s_complementary_basepair[256];
static struct Complementary_basepair_initializer {
	Complementary_basepair_initializer();
	static bool s_initialized;
} __complementary_basepair_initializer;

inline char basepair2complementary(char c) { return s_complementary_basepair[(int)c]; }

std::string seq2reverse_complementary(const std::string &seq);

#endif

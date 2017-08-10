#include "GenomeUtils.h"

using namespace std;

char s_complementary_basepair[256];

bool Complementary_basepair_initializer::s_initialized = false;

Complementary_basepair_initializer::Complementary_basepair_initializer()
{
	if (!s_initialized) {
		s_initialized = true;

		for (unsigned i = 0; i < sizeof(s_complementary_basepair) / sizeof(s_complementary_basepair[0]); i++)
			s_complementary_basepair[i] = (char)i;

		s_complementary_basepair[(int)'a'] = 't';
		s_complementary_basepair[(int)'c'] = 'g';
		s_complementary_basepair[(int)'g'] = 'c';
		s_complementary_basepair[(int)'t'] = 'a';
		s_complementary_basepair[(int)'A'] = 'T';
		s_complementary_basepair[(int)'C'] = 'G';
		s_complementary_basepair[(int)'G'] = 'C';
		s_complementary_basepair[(int)'T'] = 'A';
	}
}

std::string seq2reverse_complementary(const std::string &seq)
{
	string res;
	res.resize(seq.length());
	string::iterator ires = res.begin();
	for (string::const_reverse_iterator iseq = seq.rbegin(); iseq != seq.rend(); ++iseq)
		*ires++ = basepair2complementary(*iseq);
	return res;
}

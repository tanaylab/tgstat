#ifndef GENOMECHROMKEY_H_
#define GENOMECHROMKEY_H_

#include <stdint.h>
#include <ext/hash_map>
#include "HashFunc.h"
#include "TGLException.h"

using namespace std;
using namespace __gnu_cxx;

// -------------------- GenomeChromKey  -----------------------
// GenomeChromKey manages the mapping between chromosome id and chromosome name.
// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

class GenomeChromKey {
public:
	enum Errors { CHROM_EXISTS, CHROM_NOEXISTS, ID_NOEXISTS, FILE_READ_FAILED, BAD_FILE_FORMAT, NUM_ERRORS };

	GenomeChromKey() : m_id(0) {}

	int           chrom2id(const string &chrom) const;
	int           chrom2id(const char *chrom) const;
	const string &id2chrom(int id) const;
	uint64_t      get_chrom_size(int id) const;
	size_t        get_num_chroms() const { return m_id2chrom.size(); }

	void read_chroms_sizes_file(const char *fname);

	// returns id of the new chromosome
	int add_chrom(const string &chrom, uint64_t size);

private:
	struct Chrom {
		string   name;
		uint64_t size;

		Chrom(const string &_name, uint64_t _size) : name(_name), size(_size) {}
		bool operator<(const Chrom &chrom) const { return name < chrom.name; }
	};

	typedef __gnu_cxx::hash_map<string, int> Name2id;
	typedef std::vector<Chrom> Id2chrom;

	Name2id    m_name2id;
	Id2chrom   m_id2chrom;
	int        m_id;
};


// --------------------- implementation -----------------------

inline int GenomeChromKey::add_chrom(const string &name, uint64_t size)
{
	if (m_name2id.find(name) != m_name2id.end())
		TGLError<GenomeChromKey>(CHROM_EXISTS, "Chromosome %s already exists", name.c_str());
	m_name2id[name] = m_id;
	m_id2chrom.push_back(Chrom(name, size));
	return m_id++;
}

inline int GenomeChromKey::chrom2id(const string &name) const
{
	Name2id::const_iterator iname2id = m_name2id.find(name);
	if (iname2id == m_name2id.end())
		TGLError<GenomeChromKey>(CHROM_NOEXISTS, "Chromosome \"%s\" does not exist", name.c_str());
	return iname2id->second;
}

inline int GenomeChromKey::chrom2id(const char *name) const
{
	Name2id::const_iterator iname2id = m_name2id.find(name);
	if (iname2id == m_name2id.end())
		TGLError<GenomeChromKey>(CHROM_NOEXISTS, "Chromosome \"%s\" does not exist", name);
	return iname2id->second;
}

inline const string &GenomeChromKey::id2chrom(int id) const
{
	if (id >= (int)m_id2chrom.size())
		TGLError<GenomeChromKey>(ID_NOEXISTS, "Id %d cannot be mapped to any chromosome", id);
	return m_id2chrom[id].name;
}

inline uint64_t GenomeChromKey::get_chrom_size(int id) const
{
	if (id >= (int)m_id2chrom.size())
		TGLError<GenomeChromKey>(ID_NOEXISTS, "Id %d cannot be mapped to any chromosome", id);
	return m_id2chrom[id].size;
}

#endif /* GENOMECHROMKEY_H_ */

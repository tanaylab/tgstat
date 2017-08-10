#ifndef seqpack_PssmSet_h
#define seqpack_PssmSet_h 1

/*=================================================
=================================================*/

#include "DnaPSSM.h"

class PssmSet  {

protected:

	vector<string> m_pssm_name;

	vector<DnaPSSM> m_pssms;

	vector<string> m_pssm_bid;

public:

	const string &get_pssm_name(int i) const {
		return(m_pssm_name[i]);
	}

	const DnaPSSM &get_pssm(int id) const {
		return(m_pssms[id]);
	}

	const string &get_pssm_bid(int id ) const {
		return(m_pssm_bid[id]);
	}

	int id_space_size() const {
		return(m_pssms.size());
	}
	int size() const {
		return(m_pssms.size());
	}

	DnaPSSM &get_pssm(int id);

	void read(const string &pssm_key, const string &pssm_data, float prior = 0.025, bool logodds = false);
	void read(istream &pssm_key, istream &pssm_data, float prior = 0.025, bool logodds = false);
	void read_old(const string &pssm_key, const string &pssm_data, int max_range = 1000000, int with_bic_id = 0, float prior = 0.025, bool logodds = false);
	void read_old(istream &pssm_key, istream &pssm_data, int max_range = 1000000, int with_bic = 0, float prior = 0.025, bool logodds = false);
	void write(ostream &pssm_key, ostream &pssm_data);
};

#endif //PssmSet_h

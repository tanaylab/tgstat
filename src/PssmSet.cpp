#include "port.h"
BASE_CC_FILE

#include <sstream>

#include "TGLException.h"

#include "PssmSet.h"

DnaPSSM &PssmSet::get_pssm(int id)
{
	if (id < 0 || id >= (int)m_pssms.size())
		TGLError<PssmSet>("Pssm id is out of range");
	return m_pssms[id];
}

void PssmSet::read(const string &key_fn, const string &data_fn, float prior, bool logodds)
{
	ifstream key(key_fn.c_str());
	ifstream data(data_fn.c_str());

	if (!key)
		TGLError<PssmSet>("Cannot open PSSM file %s", key_fn.c_str());
	if (!data)
		TGLError<PssmSet>("Cannot open PSSM file %s", data_fn.c_str());

	read(key, data, prior, logodds);
}

//
void PssmSet::read_old(const string &key_fn, const string &data_fn, int max_range, int with_set_id, float prior, bool logodds)
{
	ifstream key(key_fn.c_str());
	ifstream data(data_fn.c_str());

	if (!key)
		TGLError<PssmSet>("Cannot open PSSM file %s", key_fn.c_str());
	if (!data)
		TGLError<PssmSet>("Cannot open PSSM file %s", data_fn.c_str());

	read_old(key, data, max_range, with_set_id, prior, logodds);
}

void PssmSet::read(istream &pssm_key, istream &pssm_data, float prior, bool logodds)
{
	int base_id = m_pssms.size();
	int id;
	string name;
	int bid;
	pssm_key >> id;
	id += base_id;
	while(pssm_key) {
		pssm_key >> name >> bid;
	    char c = pssm_key.get();
	    while(pssm_key && c != '\n') { c = pssm_key.get(); }
		if((int)m_pssm_name.size() <= id)  {
			m_pssm_name.resize(id + 1);
			m_pssms.resize(id + 1);
			m_pssm_bid.resize(id + 1);
		}
		m_pssms[id].set_bidirect(bid);
		m_pssm_name[id] =  name;

		pssm_key >> id;
		id += base_id;
	}

	m_pssms.resize(m_pssm_name.size());

	int pos;
	float pa,pc,pg,pt;
	pssm_data >> id;
	id += base_id;
	while(pssm_data) {
		pssm_data >> pos >> pa >> pc >> pg >> pt;

		if ((int)id >= (int)m_pssms.size()) {
			ostringstream str;
			str << "Read data on pssm id %d" << id << " pos " << pos << " but max id is " << m_pssms.size();
			TGLError<PssmSet>("%s", str.str().c_str());
		}

		m_pssms[id].resize(pos + 1);
		if(!logodds) {
			m_pssms[id][pos] = DnaProbVec(pa + prior, pc + prior, pg + prior, pt + prior);
		} else {
			m_pssms[id][pos] = DnaProbVec();
			m_pssms[id][pos].reset_log_odds(pa, pc, pg, pt);
		}

		pssm_data >> id;
		id += base_id;
	}
}

void PssmSet::read_old(istream &pssm_key, istream &pssm_data, int max_range, int with_set_id, float prior, bool logodds)
{
	int base_id = m_pssms.size();
	cerr << "read with set id " << with_set_id << endl;
	int id;
	string name;
	int fr, to, bid;
	pssm_key >> id;
	id += base_id;
	string set_id;
	float score;
	while(pssm_key) {
		if(with_set_id) {
			pssm_key >> set_id;
		}
		pssm_key >> name >> fr >> to;
		if(with_set_id) {
			pssm_key >> score;
		}
	    pssm_key >> bid;
	    char c = pssm_key.get();
	    while(pssm_key && c != '\n') { c = pssm_key.get(); }
		if((int)m_pssm_name.size() <= id)  {
			m_pssm_name.resize(id + 1);
			m_pssms.resize(id + 1);
			m_pssm_bid.resize(id + 1);
		}
		cerr << "read id " << id << " set " << set_id << " sc " << score << " direct " << bid << endl;

		m_pssms[id].set_range(max_range + fr, max_range + to);
		m_pssms[id].set_bidirect(bid);
		m_pssm_bid[id] = set_id;

		m_pssm_name[id] =  name;

		pssm_key >> id;
		id += base_id;
	}

	cerr << "done reading " << m_pssm_name.size() << " pssms " << endl;
	cerr << "logodds is " << logodds << endl;

	m_pssms.resize(m_pssm_name.size());

	int pos;
	float pa,pc,pg,pt;
	pssm_data >> id;
	id += base_id;
	while(pssm_data) {
		pssm_data >> pos >> pa >> pc >> pg >> pt;

		if ((int)id >= (int)m_pssms.size()) {
			ostringstream str;
			str << "Read data on pssm id %d" << id << " pos " << pos << " but max id is " << m_pssms.size();
			TGLError<PssmSet>("%s", str.str().c_str());
		}

		m_pssms[id].resize(pos + 1);
		if(!logodds) {
			m_pssms[id][pos] = DnaProbVec(pa + prior, pc + prior, pg + prior, pt + prior);
		} else {
			m_pssms[id][pos] = DnaProbVec();
			m_pssms[id][pos].reset_log_odds(pa, pc, pg, pt);
		}

		pssm_data >> id;
		id += base_id;
	}
}

void PssmSet::write(ostream &pssm_key, ostream &pssm_data)
{

}

#ifndef seqpack_DnaPSSM_h
#define seqpack_DnaPSSM_h 1

/*=================================================
=================================================*/

#include <string>
#include <list>
#include "RaList.h"

class DnaProbVec {
private:
	float m_p[4];
	float m_logp[4];

public:
	int encode(char c) const {
		switch(c) {
			case 'A' : return(0);
				   break;
			case 'C' : return(1);
				   break;
			case 'G' : return(2);
				   break;
			case 'T' : return(3);
				   break;
			default:   return(-1);
		}
	}
	DnaProbVec() {
		m_p[0] = 0.25;
		m_p[1] = 0.25;
		m_p[2] = 0.25;
		m_p[3] = 0.25;
		normalize();
	}

	DnaProbVec(float ap, float cp, float gp, float tp) {
		m_p[0] = ap;
		m_p[1] = cp;
		m_p[2] = gp;
		m_p[3] = tp;
		normalize();
	}

	void direct_incr_weight(int code, float weight = 1) {
		m_p[code] += weight;
	}
	void direct_incr_log_weight(int code, float lweight) {
		log_sum_log(m_logp[code], lweight);
	}

	void incr_weight(char c, float weight = 1) {
		if(c!='*')
			m_p[encode(c)] += weight;
	}
	void incr_log_weight(char c, float lweight) {
		if(c!='*')
    		log_sum_log(m_logp[encode(c)], lweight);
	}
	void set_weight(char c, float w) {
		int ec = encode(c);
		if(ec != -1) {
    		m_p[ec] = w;
		} else {
			cerr << "Set weight of PSSM with wrong character " << c;
			exit(1);
		}
	}

	void reset_log_odds(float pa, float pc, float pg, float pt) {
		m_logp[0] = pa;
		m_logp[1] = pc;
		m_logp[2] = pg;
		m_logp[3] = pt;
		m_p[0] = 0;
		m_p[1] = 0;
		m_p[2] = 0;
		m_p[3] = 0;
	}
	void reset(const vector<float> &v) {
		m_p[0] = v[0];
		m_p[1] = v[1];
		m_p[2] = v[2];
		m_p[3] = v[3];
		m_logp[0] = -_REAL(MAX);
		m_logp[1] = -_REAL(MAX);
		m_logp[2] = -_REAL(MAX);
		m_logp[3] = -_REAL(MAX);
		if(m_p[0] != 0) {
			m_logp[0] = log(m_p[0]);
		}
		if(m_p[1] != 0) {
			m_logp[1] = log(m_p[1]);
		}
		if(m_p[2] != 0) {
			m_logp[2] = log(m_p[2]);
		}
		if(m_p[3] != 0) {
			m_logp[3] = log(m_p[3]);
		}
	}

	float dot(DnaProbVec &arg) {
		return(m_p[0] * arg.m_p[0] +
			m_p[1] * arg.m_p[1] +
			m_p[2] * arg.m_p[2] +
			m_p[3] * arg.m_p[3]);
	}


	void normalize();
	void normalize_log();

	float get_direct_prob(int ci) const{
		return(m_p[ci]);
	}
	void set_direct_prob(int ci, float p) {
		m_p[ci] = p;
		m_logp[ci] = log(p);
	}

	float get_prob(char c) const {
		return(m_p[encode(c)]);
	}
	float get_log_prob(char c) const {
		return(m_logp[encode(c)]);
	}

	float get_avg_log_prob() const {
		return((m_logp[0]+m_logp[1]+m_logp[2]+m_logp[3])/4);
	}

	float get_max_log_prob() const;

	float get_entropy() const;

	void print(ostream &out);
};

ostream &operator<<(ostream &out, DnaProbVec &pvec);
ostream &operator<<(ostream &out, const DnaProbVec &pvec);

class DnaPSSM  {

public:

	static const float CONSENSUS_SINGLE_THRESH;
	static const float CONSENSUS_DOUBLE_THRESH;

protected:

	vector<DnaProbVec> m_chars;

	int m_min_range;
	int m_max_range;

	bool m_bidirect;

public:

	bool is_bidirect() const {
		return(m_bidirect);
	}

	void set_bidirect(bool dir) {
		m_bidirect = dir;
	}

	void set_range(int fr, int to) {
		m_min_range = fr;
		m_max_range = to;
	}

	int get_min_range() const {
		return(m_min_range);
	}
	int get_max_range() const {
		return(m_max_range);
	}

	DnaPSSM() :
		m_min_range(0),
		m_max_range(1000000),
		m_bidirect(false)
	{}

	DnaPSSM(const DnaPSSM &other) :
		m_chars(other.m_chars),
		m_min_range(other.m_min_range),
		m_max_range(other.m_max_range),
		m_bidirect(other.m_bidirect)
	{}

	const DnaPSSM &operator=(const DnaPSSM &other);

	void init_from_seed(const string &seed, float prior);

	DnaProbVec &operator[](int id) {
		return(m_chars[id]);
	}
	const DnaProbVec &operator[](int id) const {
		return(m_chars[id]);
	}

	uint length() const {
		return(m_chars.size());
	}
	uint size() const {
		return(m_chars.size());
	}

	void resize(int new_size);

	void push_back(DnaProbVec &pv) {
		m_chars.push_back(pv);
	}

	string::const_iterator max_like_match(const std::string &target,
					float &best_p, int &best_dir) const;

    void integrate_like_seg(const char *min_i, const char *max_i, float &energy) const;
	void integrate_like(const string &target, float &energy, vector<float> *spat_dist = 0) const;
	void integrate_energy(const string &target, float &energy, vector<float> &spat_func, int spat_bin_size) const;
	float get_max_ll() const;

	void calc_like(const std::string &target, float &logp) const;
	void calc_like_rc(const std::string &target, float &logp) const;
	void calc_like(string::const_iterator &i, float &logp) const;
	void calc_like_rc(string::const_iterator &i, float &logp) const;

	void update_like_vec(const string &str, vector<float> &old_like,
				vector<float> &deltas, vector<int1> &dirs);

	void normalize();
	void normalize_logs();

	void count(string::const_iterator seq, float weight = 1, int dir = 1);
	void count_weighted(const string &target, vector<float> &wgts,
					vector<int1> &dirs, float thresh_wgt);
	void count_log_weighted(const string &target, vector<float> &wgts,
					vector<int1> &dirs, float thresh_wgt);

	void reset_prior(const vector<float> &prior);

	float dot_product(DnaPSSM &arg);
	float log_dot_product(DnaPSSM &arg);

	void write_tab(ostream &pssmd, int id) const;
	void like_thresh_match(const string &target, float thresh,
		list<int> &poss, list<float> &vals, list<int> &dirs);

	string get_consensus() const;
	void permut_randomize();
};

ostream &operator<<(ostream &out, const DnaPSSM &pat);

#endif // DnaPSSM

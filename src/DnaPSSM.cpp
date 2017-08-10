#include "port.h"
BASE_CC_FILE
#include "DnaPSSM.h"
#include "Random.h"

#include <algorithm>


void DnaProbVec::normalize()
{
	float sum = m_p[0] + m_p[1] + m_p[2] + m_p[3];

	m_p[0] /= sum;
	m_p[1] /= sum;
	m_p[2] /= sum;
	m_p[3] /= sum;

	if(m_p[0] != 0) {
		m_logp[0] = log(m_p[0]);
	} else {
		m_logp[0] = -_REAL(MAX)/100;
	}
	if(m_p[1] != 0) {
		m_logp[1] = log(m_p[1]);
	} else {
		m_logp[1] = -_REAL(MAX)/100;
	}
	if(m_p[2] != 0) {
		m_logp[2] = log(m_p[2]);
	} else {
		m_logp[2] = -_REAL(MAX)/100;
	}
	if(m_p[3] != 0) {
		m_logp[3] = log(m_p[3]);
	} else {
		m_logp[3] = -_REAL(MAX)/100;
	}
}
void DnaProbVec::normalize_log()
{
	float sum = m_logp[0];
       	log_sum_log(sum, m_logp[1]);
       	log_sum_log(sum, m_logp[2]);
       	log_sum_log(sum, m_logp[3]);

	cerr << "normalize, sum = " << sum << " 0 " << m_logp[0] << " 1 " << m_logp[1] << endl;

	m_logp[0] -= sum;
	m_logp[1] -= sum;
	m_logp[2] -= sum;
	m_logp[3] -= sum;

	m_p[0] = exp(m_logp[0]);
	m_p[1] = exp(m_logp[1]);
	m_p[2] = exp(m_logp[2]);
	m_p[3] = exp(m_logp[3]);
}

float DnaProbVec::get_entropy() const
{
	return(-m_p[0]*m_logp[0]
		- m_p[1]*m_logp[1]
		- m_p[2]*m_logp[2]
		- m_p[3]*m_logp[3]);
}

float DnaProbVec::get_max_log_prob() const
{
	float logp = m_logp[0];
	if(m_logp[1] > logp) {
		logp = m_logp[1];
	}
	if(m_logp[2] > logp) {
		logp = m_logp[2];
	}
	if(m_logp[3] > logp) {
		logp = m_logp[3];
	}

	return(logp);
}

ostream &operator<<(ostream &out, DnaProbVec &pvec)
{
	out << int(pvec.get_prob('A')*1000)/1000.0 << "\t"
	<< int(pvec.get_prob('C')*1000)/1000.0 << "\t"
	<< int(pvec.get_prob('G')*1000)/1000.0 << "\t"
	<< int(pvec.get_prob('T')*1000)/1000.0 << endl;
	return(out);
}

ostream &operator<<(ostream &out, const DnaProbVec &pvec)
{
	out << int(pvec.get_prob('A')*1000)/1000.0 << "\t"
	<< int(pvec.get_prob('C')*1000)/1000.0 << "\t"
	<< int(pvec.get_prob('G')*1000)/1000.0 << "\t"
	<< int(pvec.get_prob('T')*1000)/1000.0 << endl;
	return(out);
}


const DnaPSSM &DnaPSSM::operator=(const DnaPSSM &other)
{
	m_chars = other.m_chars;
	m_min_range = other.m_min_range;
	m_max_range = other.m_max_range;
	m_bidirect = other.m_bidirect;
	return(*this);
}
void DnaPSSM::resize(int sz)
{
	m_chars.resize(sz);
}

void DnaPSSM::init_from_seed(const string &seed, float prior)
{
	m_chars.resize(seed.length());
	vector<DnaProbVec>::iterator p = m_chars.begin();
	vector<float> back(4, prior);
	for(string::const_iterator i = seed.begin(); i != seed.end(); i++) {
		p->reset(back);
		switch(*i) {
			case 'A': p->set_direct_prob(0, 1 - 3*prior); break;
			case 'C': p->set_direct_prob(1, 1 - 3*prior); break;
			case 'G': p->set_direct_prob(2, 1 - 3*prior); break;
			case 'T': p->set_direct_prob(3, 1 - 3*prior); break;
		}
		p->normalize();
		p++;
	}
}

float DnaPSSM::get_max_ll() const
{
	float logp = 0;
	for(vector<DnaProbVec>::const_iterator p = m_chars.begin();
	    p < m_chars.end();
	    p++) {
		logp += p->get_max_log_prob();
	}
	return(logp);

}

void DnaPSSM::calc_like(const string &target, float &logp) const
{

	string::const_iterator i = target.begin();
	logp = 0;
	for(vector<DnaProbVec>::const_iterator p = m_chars.begin();
	    p < m_chars.end();
	    p++) {
		if(*i != 'A' && *i != 'C' && *i != 'G' && *i != 'T') {
			logp = -_REAL(MAX);
			return;
		}
		logp += p->get_log_prob(*i);
		i++;
	}
}
void DnaPSSM::calc_like_rc(const string &target, float &logp) const
{

	string::const_iterator i = target.begin();
	logp = 0;
	for(vector<DnaProbVec>::const_reverse_iterator p = m_chars.rbegin();
	    p != m_chars.rend();
	    p++) {
		char c;
		switch(*i) {
			case 'A': c = 'T';
				  break;
			case 'T': c = 'A';
				  break;
			case 'C': c = 'G';
				  break;
			case 'G': c = 'C';
				  break;
			default:  logp = -_REAL(MAX);
				  return;
		}
		logp += p->get_log_prob(c);
		i++;
	}
}
void DnaPSSM::calc_like(string::const_iterator &j, float &logp) const
{

	logp = 0;
	string::const_iterator i = j;
	for(vector<DnaProbVec>::const_iterator p = m_chars.begin();
	    p < m_chars.end();
	    p++) {
		if(*i != 'A' && *i != 'C' && *i != 'G' && *i != 'T') {
			logp = -_REAL(MAX);
			return;
		}
		logp += p->get_log_prob(*i);
		i++;
	}
}

static const float c_log_quarter = -1.38629;

void DnaPSSM::calc_like_rc(string::const_iterator &j, float &logp) const
{

	logp = 0;
	string::const_iterator i = j;
	for(vector<DnaProbVec>::const_reverse_iterator p = m_chars.rbegin();
	    p != m_chars.rend();
	    p++) {
		char c;
		switch(*i) {
			case 'A': c = 'T';
				  break;
			case 'T': c = 'A';
				  break;
			case 'C': c = 'G';
				  break;
			case 'G': c = 'C';
				  break;
			default:  logp = -_REAL(MAX);
				  return;
		}
		logp += p->get_log_prob(c);
		i++;
	}
}

string::const_iterator DnaPSSM::max_like_match(const string &target,
				float &best_logp, int &best_dir) const
{
	if(target.length() < m_chars.size()) {
		best_logp = -_REAL(MAX);
		return(target.begin());
	}

	string::const_iterator max_i = target.begin() + m_max_range;
	if(max_i > target.end() - m_chars.size()) {
		max_i = target.end() - m_chars.size();
	}
	string::const_iterator best_pos;
	best_logp = -_REAL(MAX)/100;
	for(string::const_iterator i = target.begin() + m_min_range;
	    i < max_i;
	    i++) {
		string::const_iterator j = i;
		float logp = 0;
		for(vector<DnaProbVec>::const_iterator p = m_chars.begin();
		    p < m_chars.end();
		    p++) {
			if(!(*j)) {
				logp = -_REAL(MAX);
				break;
			}
			if(*j == 'N' || *j =='*') {
				logp += p->get_avg_log_prob();
			} else {
				logp += p->get_log_prob(*j);
			}
			if(logp < best_logp) {
				break;
			}
			j++;
		}
		if(logp > best_logp) {
			best_logp = logp;
			best_pos = i;
			best_dir = 1;
		}
		if(m_bidirect) {
			logp = 0;
			j = i;
			for(vector<DnaProbVec>::const_reverse_iterator p = m_chars.rbegin();
			    p != m_chars.rend();
			    p++) {
				if(!(*j)) {
					logp = -_REAL(MAX);
					break;
				}
				char c = 0;
				switch(*j) {
					case 'A': logp += p->get_log_prob('T');
						  break;
					case 'T': logp += p->get_log_prob('A');
						  break;
					case 'C': logp += p->get_log_prob('G');
						  break;
					case 'G': logp += p->get_log_prob('C');
						  break;
					case '*': logp += p->get_avg_log_prob();
						  break;
					case 'N': logp += p->get_avg_log_prob();
						  break;
					default:  break;
				}
				j++;
			}
			if(logp > best_logp) {
				best_logp = logp;
				best_pos = i;
				best_dir = -1;
			}
		}
	}
	return(best_pos);
}

void DnaPSSM::update_like_vec(const string &target,
			vector<float> &likes, vector<float> &deltas,
			vector<int1> &dirs)
{
	if(target.length() < m_chars.size()) {
		return;
	}

	string::const_iterator max_i = target.begin() + m_max_range;
	if(max_i > target.end() - m_chars.size()) {
		max_i = target.end() - m_chars.size();
	}
	vector<float>::iterator delta = deltas.begin() + m_min_range;
	vector<float>::iterator like = likes.begin() + m_min_range;
	vector<int1>::iterator dir = dirs.begin() + m_min_range;
	for(string::const_iterator i = target.begin() + m_min_range;
	    i < max_i;
	    i++) {
		string::const_iterator j = i;
		float logp = 0;
		for(vector<DnaProbVec>::const_iterator p = m_chars.begin();
		    p < m_chars.end();
		    p++) {
			if(!(*j)) {
				logp = -_REAL(MAX);
				break;
			}
			if(*j == 'N' || *j =='*') {
				logp += c_log_quarter;
			} else {
				logp += p->get_log_prob(*j);
			}
			j++;
		}
		*dir = 1;
		if(m_bidirect) {
			float rlogp = 0;
			j = i;
			for(vector<DnaProbVec>::reverse_iterator p =
							m_chars.rbegin();
			    p != m_chars.rend();
			    p++) {
				if(!(*j)) {
					rlogp = -_REAL(MAX);
					break;
				}
				char c = 0;
				switch(*j) {
					case 'A': rlogp += p->get_log_prob('T');
						  break;
					case 'T': rlogp += p->get_log_prob('A');
						  break;
					case 'C': rlogp += p->get_log_prob('G');
						  break;
					case 'G': rlogp += p->get_log_prob('C');
						  break;
					case '*': rlogp += p->get_avg_log_prob();
						  break;
					case 'N': rlogp += p->get_avg_log_prob();
						  break;
					default:  break;
				}
				j++;
			}
			if(rlogp > logp) {
				logp = rlogp;
				*dir = -1;
			}
		}
		if(logp == -_REAL(MAX)) {
			*delta = -_REAL(MAX);
			*like = -_REAL(MAX);
		} else {
			*delta = -(*like);
			*delta += logp;
			*like = logp;
		}
		like++;
		delta++;
		dir++;
	}
}

void DnaPSSM::integrate_like_seg(const char *min_i, const char *max_i, float &energy) const
{
	energy = -_REAL(MAX)/100;
	for(const char *i = min_i;
	    i < max_i;
	    i++) {
		const char *j = i;
		float logp = 0;
		for(vector<DnaProbVec>::const_iterator p = m_chars.begin();
		    p < m_chars.end();
		    p++) {
			if(!(*j)) {
				logp = -_REAL(MAX);
				break;
			}
			if(*j == 'N' || *j == 'n' || *j =='*') {
				logp += p->get_avg_log_prob();
			} else if(*j > 'Z') {
				logp += p->get_log_prob(*j-('a'-'A'));
			} else {
				logp += p->get_log_prob(*j);
			}
			j++;
		}
   		log_sum_log(energy, logp);
//   		cerr << " +e at " << i-max_i << " is " << logp << endl;
		if(m_bidirect) {
			logp = 0;
			j = i;
			for(vector<DnaProbVec>::const_reverse_iterator
							p = m_chars.rbegin();
			    p != m_chars.rend();
			    p++) {
				if(!(*j)) {
					logp = -_REAL(MAX);
					break;
				}
				char c = 0;
				switch(*j) {
					case 'A': logp += p->get_log_prob('T');
						  break;
					case 'T': logp += p->get_log_prob('A');
						  break;
					case 'C': logp += p->get_log_prob('G');
						  break;
					case 'G': logp += p->get_log_prob('C');
						  break;
					case 'a': logp += p->get_log_prob('T');
						  break;
					case 't': logp += p->get_log_prob('A');
						  break;
					case 'c': logp += p->get_log_prob('G');
						  break;
					case 'g': logp += p->get_log_prob('C');
						  break;
					case '*': logp += p->get_avg_log_prob();
						  break;
					case 'N': logp += p->get_avg_log_prob();
						  break;
					case 'n': logp += p->get_avg_log_prob();
						  break;
					default:  break;
				}
				j++;
			}
 //      		cerr << "-e at " << i-max_i << " is " << logp << endl;
       		log_sum_log(energy, logp);
		}
	}
}
void DnaPSSM::integrate_like(const string &target, float &energy, vector<float> *spat_dist) const
{
	if(target.length() < m_chars.size()) {
		energy = -_REAL(MAX);
		return;
	}

	string::const_iterator max_i = target.begin() + m_max_range;
	if(max_i > target.end() - m_chars.size()) {
		max_i = target.end() - m_chars.size();
	}
	energy = -_REAL(MAX)/100;
	for(string::const_iterator i = target.begin() + m_min_range;
	    i < max_i;
	    i++) {
		string::const_iterator j = i;
		float logp = 0;
		for(vector<DnaProbVec>::const_iterator p = m_chars.begin();
		    p < m_chars.end();
		    p++) {
			if(!(*j)) {
				logp = -_REAL(MAX);
				break;
			}
			if(*j == 'N' || *j =='*') {
				logp += c_log_quarter;
			} else {
				logp += p->get_log_prob(*j);
			}
			j++;
		}
   		log_sum_log(energy, logp);
		if(spat_dist) {
			log_sum_log((*spat_dist)[i - target.begin()], logp);
		}
		if(m_bidirect) {
			logp = 0;
			j = i;
			for(vector<DnaProbVec>::const_reverse_iterator
							p = m_chars.rbegin();
			    p != m_chars.rend();
			    p++) {
				if(!(*j)) {
					logp = -_REAL(MAX);
					break;
				}
				char c = 0;
				switch(*j) {
					case 'A': logp += p->get_log_prob('T');
						  break;
					case 'T': logp += p->get_log_prob('A');
						  break;
					case 'C': logp += p->get_log_prob('G');
						  break;
					case 'G': logp += p->get_log_prob('C');
						  break;
					case '*': logp += p->get_avg_log_prob();
						  break;
					case 'N': logp += p->get_avg_log_prob();
						  break;
					default:  break;
				}
				j++;
			}
   			log_sum_log(energy, logp);
			if(spat_dist) {
				log_sum_log((*spat_dist)[i - target.begin()], logp);
			}
		}
	}
}

void DnaPSSM::count(string::const_iterator seq, float weight, int dir)
{
	if(dir == 1) {
		for(vector<DnaProbVec>::iterator i = m_chars.begin();
		    i != m_chars.end();
		    i++) {
			i->incr_weight(*seq, weight);
			++seq;
		}
	} else {
		for(vector<DnaProbVec>::reverse_iterator i = m_chars.rbegin();
		    i != m_chars.rend();
		    i++) {
			char c;
			switch(*seq) {
				case 'A': c = 'T';
					  break;
				case 'T': c = 'A';
					  break;
				case 'C': c = 'G';
					  break;
				case 'G': c = 'C';
					  break;
				default: break;
			}
			i->incr_weight(c, weight);
			++seq;
		}
	}
}

void DnaPSSM::count_weighted(const string &target, vector<float> &wgts,
				vector<int1> &dirs, float thresh_wgt)
{
	//iterate on the correct range
	if(target.length() < m_chars.size()) {
		return;
	}

	string::const_iterator max_i = target.begin() + m_max_range;
	if(max_i > target.end() - m_chars.size()) {
		max_i = target.end() - m_chars.size();
	}
	vector<float>::iterator wgt = wgts.begin() + m_min_range;
	vector<int1>::iterator dir = dirs.begin() + m_min_range;
	for(string::const_iterator i = target.begin() + m_min_range;
	    i < max_i;
	    i++) {
		//if logp is very small - ignore it
		if(*wgt < thresh_wgt) {
			wgt++;
			dir++;
			continue;
		}
		string::const_iterator j = i;
		if(*dir == 1) {
			for(vector<DnaProbVec>::iterator p = m_chars.begin();
			    p < m_chars.end();
			    p++) {
				if((*j) && *j != 'N' && *j !='*') {
					p->incr_weight(*j, *wgt);
				}
				j++;
			}
		} else {
			for(vector<DnaProbVec>::reverse_iterator p =
							m_chars.rbegin();
			    p != m_chars.rend();
			    p++) {
				switch(*j) {
					case 'A': p->direct_incr_weight(3, *wgt);
						  break;
					case 'T': p->direct_incr_weight(0, *wgt);
						  break;
					case 'C': p->direct_incr_weight(2, *wgt);
						  break;
					case 'G': p->direct_incr_weight(1, *wgt);
						  break;
					default:  break;
				}
				j++;
			}
		}
		wgt++;
		dir++;
	}
}
void DnaPSSM::count_log_weighted(const string &target, vector<float> &wgts,
				vector<int1> &dirs, float thresh_wgt)
{
	//iterate on the correct range
	if(target.length() < m_chars.size()) {
		return;
	}

	string::const_iterator max_i = target.begin() + m_max_range;
	if(max_i > target.end() - m_chars.size()) {
		max_i = target.end() - m_chars.size();
	}
	vector<float>::iterator wgt = wgts.begin() + m_min_range;
	vector<int1>::iterator dir = dirs.begin() + m_min_range;
	for(string::const_iterator i = target.begin() + m_min_range;
	    i < max_i;
	    i++) {
		//if logp is very small - ignore it
		if(*wgt < thresh_wgt) {
			wgt++;
			dir++;
			continue;
		}
		string::const_iterator j = i;
		if(*dir == 1) {
			for(vector<DnaProbVec>::iterator p = m_chars.begin();
			    p < m_chars.end();
			    p++) {
				if((*j) && *j != 'N' && *j !='*') {
					p->incr_log_weight(*j, *wgt);
				}
				j++;
			}
		} else {
			for(vector<DnaProbVec>::reverse_iterator p =
							m_chars.rbegin();
			    p != m_chars.rend();
			    p++) {
				switch(*j) {
					case 'A': p->direct_incr_log_weight(3, *wgt);
						  break;
					case 'T': p->direct_incr_log_weight(0, *wgt);
						  break;
					case 'C': p->direct_incr_log_weight(2, *wgt);
						  break;
					case 'G': p->direct_incr_log_weight(1, *wgt);
						  break;
					default:  break;
				}
				j++;
			}
		}
		wgt++;
		dir++;
	}
}

void DnaPSSM::normalize()
{
	for(vector<DnaProbVec>::iterator i = m_chars.begin();
	    i != m_chars.end();
	    i++) {
		i->normalize();
	}
}
void DnaPSSM::normalize_logs()
{
	for(vector<DnaProbVec>::iterator i = m_chars.begin();
	    i != m_chars.end();
	    i++) {
		i->normalize_log();
	}
}

void DnaPSSM::reset_prior(const vector<float> &prior)
{
	for(vector<DnaProbVec>::iterator i = m_chars.begin();
	    i != m_chars.end();
	    i++) {
		i->reset(prior);
	}
}

//currently assuming same length profiles
float DnaPSSM::dot_product(DnaPSSM &arg)
{
	ASSERT(arg.size() == size(), "dot product support equal sized profiles, extend the code if you ment something else");
	vector<DnaProbVec>::iterator j = arg.m_chars.begin();
	float prod = 1;
	for(vector<DnaProbVec>::iterator i = m_chars.begin();
	    i != m_chars.end();
	    i++) {
		prod *= (*i).dot(*j);
		j++;
	}
	return(prod);
}
float DnaPSSM::log_dot_product(DnaPSSM &arg)
{
	ASSERT(arg.size() == size(), "dot product support equal sized profiles, extend the code if you ment something else");
	vector<DnaProbVec>::iterator j = arg.m_chars.begin();
	float prod = 1;
	for(vector<DnaProbVec>::iterator i = m_chars.begin();
	    i != m_chars.end();
	    i++) {
		prod *= (*i).dot(*j);
		j++;
	}
	return(log(prod));
}

void DnaPSSM::permut_randomize()
{
	int max_i = m_chars.size();
	for(int count = 0; count < max_i*2; count++) {
		int i = int(Random::fraction() * max_i);
		int j = int(Random::fraction() * max_i);
		DnaProbVec tmp = m_chars[i];
		m_chars[i] = m_chars[j];
		m_chars[j] = tmp;
	}
}

void DnaPSSM::write_tab(ostream &pssmd, int id) const
{
	int pos = 0;
	for(vector<DnaProbVec>::const_iterator i = m_chars.begin();
	    i != m_chars.end();
	    i++) {
		pssmd << id << "\t" << pos << "\t" << *i;
		pos++;
	}
}

const float DnaPSSM::CONSENSUS_SINGLE_THRESH = 0.6;
const float DnaPSSM::CONSENSUS_DOUBLE_THRESH = 0.85;

string DnaPSSM::get_consensus() const
{
	string output;
	vector<int> ps(4);
	for(vector<DnaProbVec>::const_iterator i = m_chars.begin();
	    i != m_chars.end();
	    i++) {
		ps[0] = int(1000 * i->get_prob('A'))*4;
		ps[1] = int(1000 * i->get_prob('C'))*4 + 1;
		ps[2] = int(1000 * i->get_prob('G'))*4 + 2;
		ps[3] = int(1000 * i->get_prob('T'))*4 + 3;

		sort(ps.begin(), ps.end());

		if(ps[3] > CONSENSUS_SINGLE_THRESH*4000) {
			int code = ps[3] % 4;
			switch(code) {
				case 0: output += 'A';
					break;
				case 1: output += 'C';
					break;
				case 2: output += 'G';
					break;
				case 3: output += 'T';
					break;
				default:
					break;
			}
			continue;
		}
		if(ps[3] + ps[2] >= 4000*CONSENSUS_DOUBLE_THRESH) {
			int code = (ps[3]%4) * 4 + ps[2]%4;
			switch(code) {
				case 1: output += 'M';	//AC
					break;
				case 2: output += 'R';	//AG
					break;
				case 3: output += 'W';	//AT
					break;
				case 4: output += 'M';	//CA
					break;
				case 6: output += 'S';	//CG
					break;
				case 7: output += 'Y';	//CT
					break;
				case 8: output += 'R';	//GA
					break;
				case 9: output += 'S';	//GC
					break;
				case 11: output += 'K';	//GT
					break;
				case 12: output += 'W';	//TA
					break;
				case 13: output += 'Y';	//TC
					break;
				case 14: output += 'K';	//TG
					break;
				default:
					output += 'e';
					break;
			}
			continue;
		}
		output += "*";
	}
	return(output);
}

ostream &operator<<(ostream &out, const DnaPSSM &pssm)
{
	cerr << "[" << pssm.get_min_range() << "," << pssm.get_max_range() << "] dir=" << pssm.is_bidirect() << endl;
	for(int i = 0; i < pssm.size(); i++) {
		out << pssm[i];
	}
	out << endl;
	for(int i = 0; i < pssm.size(); i++) {
		out << pssm[i].get_log_prob('A') << "\t" << pssm[i].get_log_prob('C') << "\t" << pssm[i].get_log_prob('G') << "\t" << pssm[i].get_log_prob('T') << endl;
	}
	return(out);
}

void DnaPSSM::integrate_energy(const string &target, float &energy, vector<float> &spat_func, int spat_bin_size) const
{
	if(target.length() < m_chars.size()) {
		energy = -_REAL(MAX);
		return;
	}

	string::const_iterator max_i = target.begin() + m_max_range;
	if(max_i > target.end() - m_chars.size()) {
		max_i = target.end() - m_chars.size();
	}
	energy = -_REAL(MAX)/100;
	int pos = 0;
	for(string::const_iterator i = target.begin() + m_min_range;
	    i < max_i;
	    i++) {
		int spat_bin = int(pos/spat_bin_size);
		pos++;
		string::const_iterator j = i;
		float logp = 0;
		for(vector<DnaProbVec>::const_iterator p = m_chars.begin();
		    p < m_chars.end();
		    p++) {
			if(!(*j)) {
				logp = -_REAL(MAX);
				break;
			}
			if(*j == 'N' || *j =='*') {
				logp += p->get_avg_log_prob();
			} else {
				logp += p->get_log_prob(*j);
			}
			j++;
		}
		logp += log(spat_func[spat_bin]);
		log_sum_log(energy, logp);
		if(m_bidirect) {
			logp = 0;
			j = i;
			for(vector<DnaProbVec>::const_reverse_iterator
							p = m_chars.rbegin();
			    p != m_chars.rend();
			    p++) {
				if(!(*j)) {
					logp = -_REAL(MAX);
					break;
				}
				char c = 0;
				switch(*j) {
					case 'A': logp += p->get_log_prob('T');
						  break;
					case 'T': logp += p->get_log_prob('A');
						  break;
					case 'C': logp += p->get_log_prob('G');
						  break;
					case 'G': logp += p->get_log_prob('C');
						  break;
					case '*': logp += c_log_quarter;
						  break;
					case 'N': logp += c_log_quarter;
						  break;
					default:  break;
				}
				j++;
			}
			logp += log(spat_func[spat_bin]);
			log_sum_log(energy, logp);
		}
	}
}
void DnaPSSM::like_thresh_match(const string &target, float thresh,
		list<int> &poss, list<float> &vals, list<int> &dirs)
{
	if(target.length() < m_chars.size()) {
		return;
	}

	string::const_iterator max_i = target.begin() + m_max_range;
	if(max_i > target.end() - m_chars.size()) {
		max_i = target.end() - m_chars.size();
	}
	for(string::const_iterator i = target.begin() + m_min_range;
	    i < max_i;
	    i++) {
		string::const_iterator j = i;
		float logp = 0;
		for(vector<DnaProbVec>::const_iterator p = m_chars.begin();
		    p != m_chars.end();
		    p++) {
			if(!(*j)) {
				logp = -_REAL(MAX);
				break;
			}
			if(*j == 'N' || *j =='*') {
				logp += c_log_quarter;
			} else {
				logp += p->get_log_prob(*j);
			}
			if(logp < thresh) {
				break;
			}
			j++;
		}
		if(logp > thresh) {
			poss.push_back(i - target.begin());
			dirs.push_back(1);
			vals.push_back(logp);
		}
		if(m_bidirect) {
			logp = 0;
			j = i;
			for(vector<DnaProbVec>::reverse_iterator p = m_chars.rbegin();
			    p != m_chars.rend();
			    p++) {
				if(!(*j)) {
					logp = -_REAL(MAX);
					break;
				}
				char c = 0;
				switch(*j) {
					case 'A': logp += p->get_log_prob('T');
						  break;
					case 'T': logp += p->get_log_prob('A');
						  break;
					case 'C': logp += p->get_log_prob('G');
						  break;
					case 'G': logp += p->get_log_prob('C');
						  break;
					case '*': logp += c_log_quarter;
						  break;
					case 'N': logp += c_log_quarter;
						  break;
					default:  break;
				}
				j++;
			}
			if(logp > thresh) {
				poss.push_back(i - target.begin());
				dirs.push_back(-1);
				vals.push_back(logp);
			}
		}
	}
}

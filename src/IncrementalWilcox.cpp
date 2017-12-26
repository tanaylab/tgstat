#include <cmath>

#include "IncrementalWilcox.h"

using namespace std;

const unsigned IncrementalWilcox::MIN_RELIABLE_WINSIZE = 10;

IncrementalWilcox::IncrementalWilcox(bool one_tailed)
{
	m_one_tailed = one_tailed;
	reset();
}

void IncrementalWilcox::reset()
{
	m_z = 1;
	m_pval = -2;
	m_std_dev_u = 0;
	m_mean_u = 0;

	for (int i = 0; i < 2; i++) {
		u[i] = 0;
		n[i] = 0;
		m_v2cnt[i].clear();
	}
}

void IncrementalWilcox::update(double old_v1, double new_v1, double old_v2, double new_v2)
{
	double old_vs[] = { old_v1, old_v2 };
	double new_vs[] = { new_v1, new_v2 };
	double old_ns[] = { n[0], n[1] };
	bool vals_differ = false;

	for (int i = 0; i < 2; i++) {
		double old_v = old_vs[i];
		double new_v = new_vs[i];

		// Basiacally we want to say: "if (old_v == new_v)" however if one of the values is NaN
		// we will get an exception in naive comparison. We also can't use "if (islessgreater(old_v, new_v))"
		// since islessgreater(NaN, NaN) returns false.
		if ((std::isnan(old_v) && std::isnan(new_v)) || (!std::isnan(old_v) && !std::isnan(new_v) && old_v == new_v))
			continue;

		int j = 1 - i;
		map<double, unsigned>::iterator old_itr = m_v2cnt[j].end();
		map<double, unsigned>::iterator new_itr = m_v2cnt[j].end();

		vals_differ = true;

		if (!std::isnan(old_v)) {
			old_itr = m_v2cnt[i].find(old_v);
			if (old_itr->second == 1)
				m_v2cnt[i].erase(old_itr);
			else
				old_itr->second--;

			old_itr = m_v2cnt[j].lower_bound(old_v);
			n[i]--;
		}

		if (!std::isnan(new_v)) {
			new_itr = m_v2cnt[i].find(new_v);
			if (new_itr == m_v2cnt[i].end())
				m_v2cnt[i][new_v] = 1;
			else
				new_itr->second++;

			new_itr = m_v2cnt[j].lower_bound(new_v);
			n[i]++;
		}

		map<double, unsigned>::iterator itr;

		// update U:
		//   For each x inserted: Uy = Uy + |y, y>x| - 0.5 * |y, y=x|
		//   For each x deleted:  Uy = Uy - |y, y>x| - 0.5 * |y, y=x|
		// Having Uy we can calculate Ux:
		//    Ux = n1*n2 - Uy
		//
		// Since we simultaneously delete x1 and insert x2 (assume x1 < x2):
		//   Uy = Uy - |y, x1<y<x2| - 0.5 * |y, y=x1 OR y=x2|
		// Similarly if x1 > x2
		//   Uy = Uy + |y, x1<y<x2| + 0.5 * |y, y=x1 OR y=x2|

		if (old_itr != m_v2cnt[j].end() && (new_itr == m_v2cnt[j].end() || old_v < new_v)) {
			// old_v < new_v

			if (old_itr->first == old_v) {
				u[j] -= 0.5 * old_itr->second;
				if (old_itr != new_itr)
					++old_itr;
			}

			for (itr = old_itr; itr != new_itr; ++itr)
				u[j] -= itr->second;

			if (new_itr != m_v2cnt[j].end() && new_itr->first == new_v)
				u[j] -= 0.5 * itr->second;

		} else if (new_itr != m_v2cnt[j].end()) {
			// new_v < old_v

			if (new_itr->first == new_v) {
				u[j] += 0.5 * new_itr->second;
				if (new_itr != old_itr)
					++new_itr;
			}

			for (itr = new_itr; itr != old_itr; ++itr)
				u[j] += itr->second;

			if (old_itr != m_v2cnt[j].end() && old_itr->first == old_v)
				u[j] += 0.5 * itr->second;
		}
		u[i] = n[0] * n[1] - u[1 - i];
	}

	if (vals_differ) {
		if (n[0] < MIN_RELIABLE_WINSIZE || n[1] < MIN_RELIABLE_WINSIZE) {
			m_pval = -1;
			m_z = 1;
		}
		else {
			double U = min(u[0], u[1]);

			if (old_ns[0] != n[0] || old_ns[1] != n[1]) {
				double n_product = n[0] * n[1];
				m_mean_u = 0.5 * n_product;
				m_std_dev_u = 1. / sqrt(n_product * (n[0] + n[1] + 1) * 0.08333333333333333333333); // 0.08333333333333 = 1/12
			}

			m_z = (U - m_mean_u) * m_std_dev_u;
			m_pval = -2;  // -2 = not calculated yet
		}
	}
}

double IncrementalWilcox::pval() const
{
	if (m_pval == -2) {
		// use erfl instead of erf as the latter returns -1 when z is low enough
		m_pval = erfl(m_z * 0.707106781) + 1;  // 0.707106781 = 1 / sqrt(2)
		if (m_one_tailed)
			m_pval *= 0.5;
	}

	return m_pval;
}

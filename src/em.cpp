/***************************************************************************
 *   Copyright (C) 2006 by Reed A. Cartwright                              *
 *   reed@scit.us                                                          *
 *                                                                         *
 *   Permission is hereby granted, free of charge, to any person obtaining *
 *   a copy of this software and associated documentation files (the       *
 *   "Software"), to deal in the Software without restriction, including   *
 *   without limitation the rights to use, copy, modify, merge, publish,   *
 *   distribute, sublicense, and/or sell copies of the Software, and to    *
 *   permit persons to whom the Software is furnished to do so, subject to *
 *   the following conditions:                                             *
 *                                                                         *
 *   The above copyright notice and this permission notice shall be        *
 *   included in all copies or substantial portions of the Software.       *
 *                                                                         *
 *   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,       *
 *   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF    *
 *   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.*
 *   IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR     *
 *   OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, *
 *   ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR *
 *   OTHER DEALINGS IN THE SOFTWARE.                                       *
 ***************************************************************************/

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#define _USE_MATH_DEFINES

#include <algorithm>
#include <cfloat>
#include <iostream>

#include "em_k2p.h"

#include "series.h"

inline double pow(double p, size_t x){return std::pow(p, (int)x);}
inline double log(size_t x) {return log((double)x);}
inline double sq(double s) { return s*s; }
inline double csch(double x) { return 1.0/sinh(x); }

double g_underflow = 2.50521045e-293;
double g_overflow = 1.596672248e+293;

using namespace std;

// z = Inverse[x = -Zeta'[z]/Zeta[z]] --- better if x < 1 and Z -> infinity
const double zeta_est_point = 1.0;
const double zeta_est_data_upper[] = {
	 1.680417359204037500, -0.493398380631477840,  0.377608050241197260, -0.301902134467186440,
	 0.249853398007347650, -0.212439584794023220,  0.184494444464681050, -0.162933942397337080,
	 0.145841311210660050, -0.131979173168241980,  0.120520348279367420, -0.110894137358286860,
	 0.102695535442322360, -0.095629845948996370,  0.089477841440473970, -0.084073215603702570,
	 0.079287602481201300, -0.075020382065199870,  0.071191590224201210, -0.067736889174703800,
	 0.064603935050871820, -0.061749711362303580,  0.059138542234139880, -0.056740591939189365,
	 0.054530717559993086, -0.052487581632782734,  0.050592958627466676, -0.048831187665996256,
	 0.047188736747310720, -0.045653852866719885,  0.044216278938531035, -0.042867023050507410
};

// z = Inverse[x = -Zeta[z]/Zeta'[z]] --- better if x > 1 and Z -> 1
const double zeta_est_data_lower[] = {
	 1.680417359204037500,   0.493398380631477840,    -0.115790330390280470,    0.04008441461626977,
	-0.016427235302097935,   0.007404979234440228,    -0.003548973529314551,    0.0017751835647363359,
	-0.0009163421000515879,  0.00048457153952771037,  -0.000261186186659685,    0.00014297530344206137,
	-0.00007927275153616899, 0.000044427372473812556, -0.000025127499617590096, 0.000014324210323756111,
	-8.221891812351182e-6,   4.7477719420294596e-6,   -2.7562873430247877e-6,   1.6077701677387667e-6,
	-9.418361102383064e-7,   5.538560303121795e-7,    -3.268378316914802e-7,    1.9348452404924703e-7,
	-1.1487368890742988e-7,  6.838357015243887e-8,    -4.080830139851217e-8,    2.4407869107849063e-8,
	-1.4629306380255232e-8,  8.785483531996262e-9,    -5.285642575287723e-9,    3.185429137389365e-9
};

taylor_series<double> zeta_est_up(zeta_est_point, zeta_est_data_upper);
taylor_series<double> zeta_est_lo(zeta_est_point, zeta_est_data_lower);

inline double zeta_est(double s) {
	return (s > 1.0) ? zeta_est_lo(1.0/s) : zeta_est_up(s);
}

extern const int g_pupy[];
extern const int g_atgc[];

/***************************************************************************
 * class em_k2p_zeta                                                       *
 ***************************************************************************/

void em_k2p_zeta::preallocate(size_t maxa, size_t maxd)
{
	sz_height = maxa+1;
	sz_width = maxd+1;
	
	size_t sz_max = std::max(sz_height, sz_width);
	cache_log.resize(sz_max, 0.0);
	for(size_t u=1;u<sz_max;++u)
		cache_log[u] = log(u);

}

void em_k2p_zeta::expectation_setup(const params_type &params)
{

}

double em_k2p_zeta::expectation(ex_type &ex, const sequence &seq_a, const sequence &seq_d, size_t index) const
{
	size_t sz_anc = seq_a.size();
	size_t sz_dec = seq_d.size();
	
	const model_type &model = get_model();
	
	aln_table table_w(sz_anc+1,sz_dec+1,ublas::zero_vector<double>(eSize));
	
	table_w(0,0)[eW] = prob_scale;
	double dt,dp;
	
	for(size_t d = 1; d <= sz_dec; ++d)
	{
		for(size_t k = d; k > 0; --k)
		{
			dp = model.p_indel_size[k];
			dt = table_w(0,d-k)[eW]*dp;
			noalias(table_w(0,d)) += table_w(0,d-k)*dp;
			table_w(0,d)[eG] += dt;
			table_w(0,d)[eL] += cache_log[k]*dt;
		}
	}
		
	for(size_t a = 1; a <= sz_anc; ++a)
	{
		for(size_t k = a; k > 0; --k)
		{
			dp = model.p_indel_size[k];
			dt = table_w(a-k,0)[eW]*dp;
			noalias(table_w(a,0)) += table_w(a-k,0)*dp;
			table_w(a,0)[eG] += dt;
			table_w(a,0)[eL] += cache_log[k]*dt;
		}
	}

	int nuc_a, nuc_d;
		
	for(size_t a=1;a<=sz_anc;++a)
	{
		for(size_t d=1;d<= sz_dec;++d)
		{
			nuc_a = seq_a[a-1];
			nuc_d = seq_d[d-1];
				
			ex_type &W = table_w(a,d);
			noalias(W) = table_w(a-1,d-1)*model.p_substitution[nuc_a][nuc_d];
			if(nuc_a == nN || nuc_d == nN)
			{
				W[eM] += model.p_match*W[eW];
				W[eS] += model.p_ts*W[eW];
				W[eV] += model.p_tv*W[eW];
			}
			else if(nuc_a == nuc_d)
				W[eM] += W[eW];
			else if(g_pupy[nuc_a] == g_pupy[nuc_d])
				W[eS] += W[eW];
			else
				W[eV] += W[eW];
						
			for(size_t k = d; k > 0; --k)
			{
				dp = model.p_indel_size[k];
				dt = table_w(a,d-k)[eW]*dp;
				noalias(W) += table_w(a,d-k)*dp;
				W[eG] += dt;
				W[eL] += cache_log[k]*dt;
			}			
			for(size_t k = a; k > 0; --k)
			{		
				dp = model.p_indel_size[k];
				dt = table_w(a-k,d)[eW]*dp;
				noalias(W) += table_w(a-k,d)*dp;
				W[eG] += dt;
				W[eL] += cache_log[k]*dt;
			}
		}
	}
	
	ex = table_w(sz_anc,sz_dec);
	double w = ex[eW];
	ex /= w;
	return log(w) + log(model.p_end) - log(prob_scale)
		+ static_cast<double>(sz_anc+sz_dec)*log(model.nuc_scale)
		+ static_cast<double>(get_seq_info()[index][nN])*log(model.amb_scale);
}

void em_k2p_zeta::maximization(params_type &p, const exvec_type &ev)
{	
	ex_type ex;
	ex.clear();
	ex = for_each(ev.begin(), ev.end(), sum_elements<ex_type>(ex)).sum;
	
	double W = ex[eW];
	double M = ex[eM];
	double S = ex[eS];
	double V = ex[eV];
	double G = ex[eG];
	double L = ex[eL];
	
	double N = M+S+V;
	double A = G+N;
	double P = S/N;
	double Q = V/N;
	double e1 = log(1.0-2.0*Q);
	double e2 = log(1.0-2.0*P-Q);
	
	double t = -0.25*(e1+2.0*e2);
	double k = e2/e1-0.5;
	
	double r = -0.5/t * log(1.0 - G/A);
	double z = zeta_est(L/G);
	double a = A/W;
	
	if(t == 0.0)
		t = k = 0.0;
	else if(!std::isfinite(k))
		k = DBL_MAX;	
	if(!std::isfinite(r))
		r = -0.5 * log(1.0 - G/A);
	if(r == 0.0) {
		r = 0.0;
		z = 2.0;
	}
	
	
	p[pT] = t;
	p[pK] = k;
	p[pR] = r;
	p[pZ] = z;
	p[pA] = a;
}

std::string em_k2p_zeta::state() const
{
	const params_type &pm = get_params();
	const exvec_type &ev = get_exvec();
	ex_type ex;
	ex.clear();
	ex = for_each(ev.begin(), ev.end(), sum_elements<ex_type>(ex)).sum;
	
	double llrs = (get_loglh() - log(0.25)*static_cast<double>(total_nucs));
	
	ostringstream ostr;
	ostr << std::setprecision(DBL_DIG);
	ostr << "Expectations:" << endl;
	ostr << "  adj-log-likelihood = " << llrs << endl;
	ostr << "  num-seq-pairs = " << ex[eW] << endl;
	ostr << "  matches = " << ex[eM] << endl;
	ostr << "  transitions = " << ex[eS] << endl;
	ostr << "  transversions = " << ex[eV] << endl;
	ostr << "  num-gaps = " << ex[eG] << endl;
	ostr << "  avg-ln-gap-size = " << ex[eL]/ex[eG] << endl;
	ostr << endl;
	ostr << "Estimated Parameters:" << endl;
	ostr << "  model = " << "zeta" << endl;
	ostr << "  scale = " << scale << endl;
	ostr << "  branch-length = " << pm[pT] << endl;
	ostr << "  ratio = " << pm[pK] << endl;
	ostr << "  indel-rate = " << pm[pR] << endl;
	ostr << "  indel-slope = " << pm[pZ] << endl;
	ostr << "  avgaln = " << pm[pA] << endl;
	return ostr.str();
}

/***************************************************************************
 * class em_k2p_geo                                                        *
 ***************************************************************************/

void em_k2p_geo::preallocate(size_t maxa, size_t maxd)
{
	sz_height = maxa+1;
	sz_width = maxd+1;
	
	size_t sz_max = std::max(sz_height, sz_width);
	cache_size.resize(sz_max, 0.0);
	for(size_t u=1;u<sz_max;++u)
		cache_size[u] = double(u);
}

void em_k2p_geo::expectation_setup(const params_type &params)
{

}

double em_k2p_geo::expectation(ex_type &ex, const sequence &seq_a, const sequence &seq_d, size_t index) const
{
	size_t sz_anc = seq_a.size();
	size_t sz_dec = seq_d.size();
	
	const model_type &model = get_model();
	
	aln_table table_w(sz_anc+1,sz_dec+1,ublas::zero_vector<double>(eSize));
	
	table_w(0,0)[eW] = prob_scale;
	
	ex_type row_cache;
	vector<ex_type> col_cache(sz_dec+1, ublas::zero_vector<double>(eSize));

	row_cache.clear();
	for(size_t d = 1; d <= sz_dec; ++d)
	{
		row_cache *= model.p_extend;
		noalias(row_cache) += table_w(0,d-1)*model.p_open;
		row_cache[eL] += row_cache[eW];
		row_cache[eG] += table_w(0,d-1)[eW]*model.p_open;
		noalias(table_w(0,d)) += row_cache;
	}
	row_cache.clear();
	for(size_t a = 1; a <= sz_anc; ++a)
	{
		row_cache *= model.p_extend;
		noalias(row_cache) += table_w(a-1,0)*model.p_open;
		row_cache[eL] += row_cache[eW];
		row_cache[eG] += table_w(a-1,0)[eW]*model.p_open;
		noalias(table_w(a,0)) += row_cache;
	}
	
	int nuc_a, nuc_d;
		
	for(size_t a=1;a<=sz_anc;++a)
	{
		row_cache.clear();
		for(size_t d=1;d<= sz_dec;++d)
		{
			nuc_a = seq_a[a-1];
			nuc_d = seq_d[d-1];
				
			ex_type &W = table_w(a,d);
			noalias(W) = table_w(a-1,d-1)*model.p_substitution[nuc_a][nuc_d];
			if(nuc_a == nN || nuc_d == nN)
			{
				W[eM] += model.p_match*W[eW];
				W[eS] += model.p_ts*W[eW];
				W[eV] += model.p_tv*W[eW];
			}
			else if(nuc_a == nuc_d)
				W[eM] += W[eW];
			else if(g_pupy[nuc_a] == g_pupy[nuc_d])
				W[eS] += W[eW];
			else
				W[eV] += W[eW];
			
			row_cache *= model.p_extend;
			noalias(row_cache) += table_w(a,d-1)*model.p_open;
			row_cache[eL] += row_cache[eW];
			row_cache[eG] += table_w(a,d-1)[eW]*model.p_open;
			noalias(W) += row_cache;
			
			col_cache[d] *= model.p_extend;
			noalias(col_cache[d]) += table_w(a-1,d)*model.p_open;
			col_cache[d][eL] += col_cache[d][eW];
			col_cache[d][eG] += table_w(a-1,d)[eW]*model.p_open;
			noalias(W) += col_cache[d];
		}
	}
	
	ex = table_w(sz_anc,sz_dec);
	double w = ex[eW];
	ex /= w;
	return log(w) + log(model.p_end) - log(prob_scale)
		+ static_cast<double>(sz_anc+sz_dec)*log(model.nuc_scale)
		+ static_cast<double>(get_seq_info()[index][nN])*log(model.amb_scale);
}

void em_k2p_geo::maximization(params_type &p, const exvec_type &ev)
{	
	ex_type ex;
	ex.clear();
	ex = for_each(ev.begin(), ev.end(), sum_elements<ex_type>(ex)).sum;
	
	double W = ex[eW];
	double M = ex[eM];
	double S = ex[eS];
	double V = ex[eV];
	double G = ex[eG];
	double L = ex[eL];

	double N = M+S+V;
	double A = G+N;
	double P = S/N;
	double Q = V/N;
	double e1 = log(1.0-2.0*Q);
	double e2 = log(1.0-2.0*P-Q);
	
	double t = -0.25*(e1+2.0*e2);
	double k = e2/e1-0.5;
	double r = -0.5/t * log(1.0 - G/A);
	double q = L/G;
	
	if(t == 0.0)
		t = k = 0.0;
	else if(!std::isfinite(k))
		k = DBL_MAX;	
	if(!std::isfinite(r))
		r = -0.5 * log(1.0 - G/A);
	if(r == 0.0) {
		r = 0.0;
		q = 1.0;
	}
	
	double a = A/W;
	
	p[pT] = t;
	p[pK] = k;
	p[pR] = r;
	p[pQ] = q;
	p[pA] = a;
}

std::string em_k2p_geo::state() const
{
	const params_type &pm = get_params();
	const exvec_type &ev = get_exvec();
	ex_type ex;
	ex.clear();
	ex = for_each(ev.begin(), ev.end(), sum_elements<ex_type>(ex)).sum;
	
	double llrs = (get_loglh() - log(0.25)*static_cast<double>(total_nucs));

	ostringstream ostr;
	ostr << std::setprecision(DBL_DIG);
	ostr << "Expectations:" << endl;
	ostr << "  adj-log-likelihood = " << llrs << endl;
	ostr << "  matches = " << ex[eM] << endl;
	ostr << "  transitions = " << ex[eS] << endl;
	ostr << "  transversions = " << ex[eV] << endl;
	ostr << "  num-gaps = " << ex[eG] << endl;
	ostr << "  avg-gap-size = " << ex[eL]/ex[eG] << endl;
	ostr << endl;
	ostr << "Estimated Parameters:" << endl;
	ostr << "  model = " << "geo" << endl;
	ostr << "  scale = " << scale << endl;	
	ostr << "  branch-length = " << pm[pT] << endl;
	ostr << "  ratio = " << pm[pK] << endl;
	ostr << "  indel-rate = " << pm[pR] << endl;
	ostr << "  indel-mean = " << pm[pQ] << endl;
	ostr << "  avgaln = " << pm[pA] << endl;
	return ostr.str();
}

/***************************************************************************
 * class em_k2p_geo_guard                                                  *
 ***************************************************************************/

void em_k2p_geo_guard::preallocate(size_t maxa, size_t maxd)
{
	sz_height = maxa+1;
	sz_width = maxd+1;
	
	size_t sz_max = std::max(sz_height, sz_width);
	cache_size.resize(sz_max, 0.0);
	for(size_t u=1;u<sz_max;++u)
		cache_size[u] = double(u);
}

void em_k2p_geo_guard::expectation_setup(const params_type &params)
{

}

double em_k2p_geo_guard::expectation(ex_type &ex, const sequence &seq_a, const sequence &seq_d, size_t index) const
{
	size_t sz_anc = seq_a.size();
	size_t sz_dec = seq_d.size();
	
	const model_type &model = get_model();
	big_prob p_extend(model.p_extend);
	big_prob p_open(model.p_open);
	
	aln_table table_w(sz_anc+1,sz_dec+1);
	
	table_w(0,0).p = prob_scale;
	
	ex_type row_cache, cache_temp;
	vector<ex_type> col_cache(sz_dec+1);

	row_cache.clear();
	for(size_t d = 1; d <= sz_dec; ++d) {
		row_cache *= p_extend;
		cache_temp = table_w(0,d-1);
		cache_temp.p.q *= model.p_open;
		cache_temp.stat[eG] += 1.0f;
		row_cache += cache_temp;
		row_cache.stat[eL] += 1.0f;
		table_w(0,d) = row_cache;
	}
	row_cache.clear();
	for(size_t a = 1; a <= sz_anc; ++a) {
		row_cache *= p_extend;
		cache_temp = table_w(a-1,0);
		cache_temp.p.q *= model.p_open;
		cache_temp.stat[eG] += 1.0f;
		row_cache += cache_temp;
		row_cache.stat[eL] += 1.0f;
		table_w(a,0) = row_cache;
	}
	
	int nuc_a, nuc_d;
		
	for(size_t a=1;a<=sz_anc;++a)
	{
		row_cache.clear();
		for(size_t d=1;d<= sz_dec;++d)
		{
			nuc_a = seq_a[a-1];
			nuc_d = seq_d[d-1];
				
			ex_type &W = table_w(a,d);
			W = table_w(a-1,d-1);
			W.p.q *= model.p_substitution[nuc_a][nuc_d];
			if(nuc_a == nN || nuc_d == nN) {
				W.stat[eM] += model.p_match;
				W.stat[eS] += model.p_ts;
				W.stat[eV] += model.p_tv;
			} else if(nuc_a == nuc_d)
				W.stat[eM] += 1.0f;
			else if(g_pupy[nuc_a] == g_pupy[nuc_d])
				W.stat[eS] += 1.0f;
			else
				W.stat[eV] += 1.0f;
			
			row_cache.p.q *= model.p_extend;
			//row_cache.stat[eL] += 1.0;
			cache_temp = table_w(a,d-1);
			cache_temp.p.q *= model.p_open;
			//cache_temp.stat[eL] += 1.0;
			cache_temp.stat[eG] += 1.0f;
			row_cache += cache_temp;
			row_cache.stat[eL] += 1.0f;
			
			col_cache[d].p.q *= model.p_extend;
			//col_cache[d].stat[eL] += 1.0;
			cache_temp = table_w(a-1,d);
			cache_temp.p.q *= model.p_open;
			//cache_temp.stat[eL] += 1.0;
			cache_temp.stat[eG] += 1.0f;
			col_cache[d] += cache_temp;
			col_cache[d].stat[eL] += 1.0f;
			
			//W.plus2(row_cache, col_cache[d]);
			W += row_cache;
			W += col_cache[d];
		}
	}
	
	ex = table_w(sz_anc,sz_dec);
	big_prob w = ex.p;
	return w.log() + log(model.p_end) - log(prob_scale)
		+ static_cast<double>(sz_anc+sz_dec)*log(model.nuc_scale)
		+ static_cast<double>(get_seq_info()[index][nN])*log(model.amb_scale);
}

void em_k2p_geo_guard::maximization(params_type &p, const exvec_type &ev)
{	
	ex_type ex;
	ex.clear();
	ex.p = 1.0;
	// will call operator *=, which will sum our summary statistics
	ex = for_each(ev.begin(), ev.end(), prod_elements<ex_type>(ex)).prod;
	
	double W = ev.size();
	double M = ex.stat[eM];
	double S = ex.stat[eS];
	double V = ex.stat[eV];
	double G = ex.stat[eG];
	double L = ex.stat[eL];

	double N = M+S+V;
	double A = G+N;
	double P = S/N;
	double Q = V/N;
	double e1 = (1.0-2.0*Q);
	double e2 = (1.0-2.0*P-Q);

	double t = -0.25*(log(e1)+2.0*log(e2));
	double k = log(e2)/log(e1)-0.5;

	double r = -0.5/t * log(1.0 - G/A);
	double q = L/G;
	double a = A/W;
	
	p[pT] = t;
	p[pK] = k;
	p[pR] = r;
	p[pQ] = q;
	p[pA] = a;
}

std::string em_k2p_geo_guard::state() const
{
	const params_type &pm = get_params();
	const exvec_type &ev = get_exvec();
	ex_type ex;
	ex.clear();
	ex.p = 1.0;
	// will call operator *=, which will sum our summary statistics
	ex = for_each(ev.begin(), ev.end(), prod_elements<ex_type>(ex)).prod;
	
	double llrs = (get_loglh() - log(0.25)*static_cast<double>(total_nucs));

	ostringstream ostr;
	ostr << std::setprecision(DBL_DIG);
	ostr << "Expectations:" << endl;
	ostr << "  adj-log-likelihood = " << llrs << endl;
	//ostr << "  num-seq-pairs = " << ex[eW] << endl;	
	ostr << "  matches = " << ex.stat[eM] << endl;
	ostr << "  transitions = " << ex.stat[eS] << endl;
	ostr << "  transversions = " << ex.stat[eV] << endl;
	ostr << "  num-gaps = " << ex.stat[eG] << endl;
	ostr << "  avg-gap-size = " << ex.stat[eL]/ex.stat[eG] << endl;
	ostr << endl;
	ostr << "Estimated Parameters:" << endl;
	ostr << "  model = " << "geo-guard" << endl;
	ostr << "  scale = " << scale << endl;	
	ostr << "  branch-length = " << pm[pT] << endl;
	ostr << "  ratio = " << pm[pK] << endl;
	ostr << "  indel-rate = " << pm[pR] << endl;
	ostr << "  indel-mean = " << pm[pQ] << endl;
	ostr << "  avgaln = " << pm[pA] << endl;
	return ostr.str();
}


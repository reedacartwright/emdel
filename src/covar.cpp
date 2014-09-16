/***************************************************************************
 *   Copyright (C) 2007 by Reed A. Cartwright                              *
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

#include <boost/math/special_functions/zeta.hpp>

#include "covar.h"
#include "covar_k2p.h"
#include "ccvector.h"
#include "series.h"
#include "invert_matrix.h"
#include "table.h"

using namespace std;
namespace ublas = boost::numeric::ublas;

inline double pow(double p, size_t x){return std::pow(p, (int)x);}
inline double log(size_t x) {return log((double)x);}
inline double zeta(double z) { return boost::math::zeta<double>(z); }
inline double sq(double s) { return s*s; }
inline double csch(double x) { return 1.0/sinh(x); }

extern const int g_pupy[5];
extern const int g_atgc[5];

// zeta_ratio(z) = zeta'(z)/zeta(z)
const double zeta_ratio_point = 1.6804173592040375;
const double zeta_ratio_data[] = {
	-1.0,                2.0267597934151023, -3.1437505990361300, 4.6584860638091640,
	-6.8551150550990680, 10.076942606572075, -14.810472757463808, 21.766886187422337,
	-31.990529063023680, 47.016048846364384, -69.098840622593630, 101.55361290858988,
	-149.25194318992214, 219.35352061052055, -322.38084118972284, 473.79867199375957,
	-696.33536767600560, 1023.3944773113722, -1504.0687358544299, 2210.5090581667810,
	-3248.7546478130210, 4774.6498584537580, -7017.2369853102955, 10313.136327854983,
	-15157.074093346790, 22276.142559146083, -32738.939208142758, 48115.967009485650,
	-70715.372496922540, 103929.40676828961, -152743.61443374804, 224485.18158388760
};
taylor_series<double> zeta_ratio(zeta_ratio_point, zeta_ratio_data);

// zeta_ratio2(x) = zeta''(z)/zeta(z)
const double zeta_ratio2_data[] = {
	3.0267597934151023,   -10.341020784902463,   24.370714649703780,   -49.480686977316730,
	92.861375476624520,   -166.09430161114878,   287.63926958655040,   -486.72025090333090,
	809.35805827060480,   -1327.7000268719173,   2154.4096216995576,   -3464.8099653836620,
	5530.8904169767300,   -8773.4349769199200,   13841.793616644849,   -21735.763734672070,
	33991.540864093986,   -52965.035843404560,   82263.000139490190,   -127398.30500095984,
	196784.82990677876,   -303246.42196821440,   466303.35285683510,   -715633.74879524200,
	1.0963093267557952e6, -1.6767085619459874e6, 2.5604673613858290e6, -3.9045142800481310e6,
	5.9462695510516780e6, -9.0446375027637000e6, 1.3741749244271364e7, -2.0855904711249072e7
};
taylor_series<double> zeta_ratio2(zeta_ratio_point, zeta_ratio2_data);

/***************************************************************************
 * class covar_k2p_zeta                                                    *
 ***************************************************************************/

void covar_k2p_zeta::preallocate(size_t maxa, size_t maxd)
{
	sz_height = maxa+1;
	sz_width = maxd+1;
		
	size_t sz_max = std::max(sz_height, sz_width);
	
	cache_log.resize(sz_max, 0.0);
	for(size_t u=1;u<sz_max;++u)
		cache_log[u] = log(u);
}

void covar_k2p_zeta::expectation_setup(const params_type &params)
{

}

double covar_k2p_zeta::expectation(ex_type &ex, const sequence &seq_a, const sequence &seq_d, size_t index) const
{
	size_t sz_anc = seq_a.size();
	size_t sz_dec = seq_d.size();
	
	const model_type &model = get_model();
		
	typedef ublas::cc_vector<double, 6> ex_s;
	typedef ublas::cc_vector<double, 5> ex_t;
	
	table<ex_s> table_s(sz_anc+1,sz_dec+1,ublas::zero_vector<double>(6));
	table<ex_t> table_t(sz_anc+1,sz_dec+1,ublas::zero_vector<double>(5));
	
	table_s(0,0)[e1W] = prob_scale;
	double dp;

	// First Cycle e1 and e2
	
	for(size_t d = 1; d <= sz_dec; ++d)
	{
		ex_s &WS = table_s(0,d);
		ex_t &WT = table_t(0,d);
		for(size_t k = d; k > 0; --k)
		{
			dp = model.p_indel_size[k];
			ex_s &WSK = table_s(0,d-k);
			noalias(WT) += table_t(0,d-k)*dp;
			noalias(WS) += WSK*dp;
			WS[e1G] += WSK[e1W]*dp;
			WS[e1L] += WSK[e1W]*dp*cache_log[k];
			WT[e2GG] += 2.0*WSK[e1G]*dp + WSK[e1W]*dp;
			WT[e2LL] += 2.0*WSK[e1L]*dp*cache_log[k] + WSK[e1W]*dp*cache_log[k]*cache_log[k];
		}
	}
		
	for(size_t a = 1; a <= sz_anc; ++a)
	{
		ex_s &WS = table_s(a,0);
		ex_t &WT = table_t(a,0);
		for(size_t k = a; k > 0; --k)
		{
			dp = model.p_indel_size[k];
			ex_s &WSK = table_s(a-k,0);
			noalias(WT) += table_t(a-k,0)*dp;
			noalias(WS) += WSK*dp;
			WS[e1G] += WSK[e1W]*dp;
			WS[e1L] += WSK[e1W]*dp*cache_log[k];
			WT[e2GG] += 2.0*WSK[e1G]*dp + WSK[e1W]*dp;
			WT[e2LL] += 2.0*WSK[e1L]*dp*cache_log[k] + WSK[e1W]*dp*cache_log[k]*cache_log[k];
			
		}
	}

	int nuc_a, nuc_d;
		
	for(size_t a=1;a<=sz_anc;++a)
	{
		for(size_t d=1;d<= sz_dec;++d)
		{
			nuc_a = seq_a[a-1];
			nuc_d = seq_d[d-1];
			dp = model.p_substitution[nuc_a][nuc_d];
				
			ex_s &WS = table_s(a,d);
			ex_s &WSP = table_s(a-1,d-1);
			ex_t &WT = table_t(a,d);
			noalias(WS) = WSP*dp;
			noalias(WT) = table_t(a-1,d-1)*dp;
			
			if(nuc_a == nN || nuc_d == nN)
			{
				WS[e1M] += WS[e1W]*model.p_match;
				WS[e1S] += WS[e1W]*model.p_ts;
				WS[e1V] += WS[e1W]*model.p_tv;
				
				WT[e2MM] += 2.0*WSP[e1M]*dp*model.p_match + WS[e1W]*model.p_match*model.p_match;
				WT[e2SS] += 2.0*WSP[e1S]*dp*model.p_ts + WS[e1W]*model.p_ts*model.p_ts;
				WT[e2VV] += 2.0*WSP[e1V]*dp*model.p_tv + WS[e1W]*model.p_tv*model.p_tv;
			}
			else if(nuc_a == nuc_d)
			{
				WS[e1M] += WS[e1W];
				WT[e2MM] += 2.0*WSP[e1M]*dp + WS[e1W];
			}
			else if(g_pupy[nuc_a] == g_pupy[nuc_d])
			{
				WS[e1S] += WS[e1W];
				WT[e2SS] += 2.0*WSP[e1S]*dp + WS[e1W];
			}
			else
			{
				WS[e1V] += WS[e1W];
				WT[e2VV] += 2.0*WSP[e1V]*dp + WS[e1W];
			}
						
			for(size_t k = d; k > 0; --k)
			{
				dp = model.p_indel_size[k];
				ex_s &WSK = table_s(a,d-k);
				noalias(WT) += table_t(a,d-k)*dp;	
				noalias(WS) += WSK*dp;
				WS[e1G] += WSK[e1W]*dp;
				WS[e1L] += WSK[e1W]*dp*cache_log[k];
				WT[e2GG] += 2.0*WSK[e1G]*dp + WSK[e1W]*dp;
				WT[e2LL] += 2.0*WSK[e1L]*dp*cache_log[k] + WSK[e1W]*dp*cache_log[k]*cache_log[k];
			}			
			for(size_t k = a; k > 0; --k)
			{		
				dp = model.p_indel_size[k];
				ex_s &WSK = table_s(a-k,d);
				noalias(WT) += table_t(a-k,d)*dp;
				noalias(WS) += WSK*dp;
				WS[e1G] += WSK[e1W]*dp;
				WS[e1L] += WSK[e1W]*dp*cache_log[k];
				WT[e2GG] += 2.0*WSK[e1G]*dp + WSK[e1W]*dp;
				WT[e2LL] += 2.0*WSK[e1L]*dp*cache_log[k] + WSK[e1W]*dp*cache_log[k]*cache_log[k];
			}
		}
	}
	ex_s ex1 = table_s(sz_anc,sz_dec);
	ex_t ex2 = table_t(sz_anc,sz_dec);
		
	double w = ex1[e1W];
	ex1 /= w;
	ex2 /= w;
	
	ex[eM] = ex1[e1M];   ex[eS] = ex1[e1S];   ex[eV] = ex1[e1V];
	ex[eG] = ex1[e1G];   ex[eL] = ex1[e1L];
	ex[eMM] = ex2[e2MM]; ex[eSS] = ex2[e2SS]; ex[eVV] = ex2[e2VV];
	ex[eGG] = ex2[e2GG]; ex[eLL] = ex2[e2LL];
	
	output() << "." << std::flush;
	
	// Second Cycle: e3: e3MS, e3MV, e3MG, e3ML, e3SV
	table_t(0,0).clear();
	for(size_t d = 1; d <= sz_dec; ++d)
	{
		ex_t &WT = table_t(0,d);
		WT.clear();
		for(size_t k = d; k > 0; --k)
		{
			dp = model.p_indel_size[k];
			ex_s &WSK = table_s(0,d-k);
			noalias(WT) += table_t(0,d-k)*dp;
			WT[e3MG] += WSK[e1M]*dp;
			WT[e3ML] += WSK[e1M]*dp*cache_log[k];
		}
	}
	for(size_t a = 1; a <= sz_anc; ++a)
	{
		ex_t &WT = table_t(a,0);
		WT.clear();
		for(size_t k = a; k > 0; --k)
		{
			dp = model.p_indel_size[k];
			ex_s &WSK = table_s(a-k,0);
			noalias(WT) += table_t(a-k,0)*dp;
			WT[e3MG] += WSK[e1M]*dp;
			WT[e3ML] += WSK[e1M]*dp*cache_log[k];
		}
	}
	for(size_t a=1;a<=sz_anc;++a)
	{
		for(size_t d=1;d<= sz_dec;++d)
		{
			nuc_a = seq_a[a-1];
			nuc_d = seq_d[d-1];
			dp = model.p_substitution[nuc_a][nuc_d];
				
			ex_t &WT = table_t(a,d);
			ex_s &WSP = table_s(a-1,d-1);
			noalias(WT) = table_t(a-1,d-1)*dp;
			
			if(nuc_a == nN || nuc_d == nN)
			{
				WT[e3MS] += WSP[e1S]*dp*model.p_match
					+ WSP[e1M]*dp*model.p_ts
					+ WSP[e1W]*dp*model.p_ts*model.p_match;
				WT[e3MV] += WSP[e1V]*dp*model.p_match
					+ WSP[e1M]*dp*model.p_tv
					+ WSP[e1W]*dp*model.p_tv*model.p_match;
				WT[e3SV] += WSP[e1V]*dp*model.p_ts
					+ WSP[e1S]*dp*model.p_tv
					+ WSP[e1W]*dp*model.p_ts*model.p_tv;
				
				WT[e3MG] += WSP[e1G]*dp*model.p_match;
				WT[e3ML] += WSP[e1L]*dp*model.p_match;
			}
			else if(nuc_a == nuc_d)
			{
				WT[e3MS] += WSP[e1S]*dp;
				WT[e3MV] += WSP[e1V]*dp;
				WT[e3MG] += WSP[e1G]*dp;
				WT[e3ML] += WSP[e1L]*dp;	
			}
			else if(g_pupy[nuc_a] == g_pupy[nuc_d])
			{
				WT[e3MS] += WSP[e1M]*dp;
				WT[e3SV] += WSP[e1V]*dp;
			}
			else
			{
				WT[e3MV] += WSP[e1M]*dp;
				WT[e3SV] += WSP[e1S]*dp;
			}
						
			for(size_t k = d; k > 0; --k)
			{
				dp = model.p_indel_size[k];
				ex_s &WSK = table_s(a,d-k);
				noalias(WT) += table_t(a,d-k)*dp;	
				WT[e3MG] += WSK[e1M]*dp;
				WT[e3ML] += WSK[e1M]*dp*cache_log[k];
				
			}			
			for(size_t k = a; k > 0; --k)
			{		
				dp = model.p_indel_size[k];
				ex_s &WSK = table_s(a-k,d);
				noalias(WT) += table_t(a-k,d)*dp;
				WT[e3MG] += WSK[e1M]*dp;
				WT[e3ML] += WSK[e1M]*dp*cache_log[k];
			}
		}
	}	
	
	ex2 = table_t(sz_anc,sz_dec);
		
	ex2 /= w;
	ex[eMS] = ex2[e3MS];
	ex[eMV] = ex2[e3MV];
	ex[eMG] = ex2[e3MG];
	ex[eML] = ex2[e3ML];
	ex[eSV] = ex2[e3SV];
	
	output() << "." << std::flush;
		
	// Third Cycle: e4: e4SG, e4SL, e4VG, e4VL, e4GL
	table_t(0,0).clear();
	for(size_t d = 1; d <= sz_dec; ++d)
	{
		ex_t &WT = table_t(0,d);
		WT.clear();
		for(size_t k = d; k > 0; --k)
		{
			dp = model.p_indel_size[k];
			noalias(WT) += table_t(0,d-k)*dp;
			ex_s &WSK = table_s(0,d-k);
			WT[e4SG] += WSK[e1S]*dp;
			WT[e4VG] += WSK[e1V]*dp;
			WT[e4SL] += WSK[e1S]*dp*cache_log[k];
			WT[e4VL] += WSK[e1V]*dp*cache_log[k];
			WT[e4GL] += WSK[e1G]*dp*cache_log[k] + WSK[e1L]*dp + WSK[e1W]*dp*cache_log[k];
		}
	}	
	for(size_t a = 1; a <= sz_anc; ++a)
	{
		ex_t &WT = table_t(a,0);
		WT.clear();
		for(size_t k = a; k > 0; --k)
		{
			dp = model.p_indel_size[k];
			noalias(WT) += table_t(a-k,0)*dp;
			ex_s &WSK = table_s(a-k,0);
			WT[e4SG] += WSK[e1S]*dp;
			WT[e4VG] += WSK[e1V]*dp;
			WT[e4SL] += WSK[e1S]*dp*cache_log[k];
			WT[e4VL] += WSK[e1V]*dp*cache_log[k];
			WT[e4GL] += WSK[e1G]*dp*cache_log[k] + WSK[e1L]*dp + WSK[e1W]*dp*cache_log[k];
		}
	}
	for(size_t a=1;a<=sz_anc;++a)
	{
		for(size_t d=1;d<= sz_dec;++d)
		{
			nuc_a = seq_a[a-1];
			nuc_d = seq_d[d-1];
				
			dp = model.p_substitution[nuc_a][nuc_d];
			ex_t &WT = table_t(a,d);
			noalias(WT) = table_t(a-1,d-1)*dp;
			ex_s &WSP = table_s(a-1,d-1);
			
			if(nuc_a == nN || nuc_d == nN)
			{
				WT[e4SG] += WSP[e1G]*dp*model.p_ts;
				WT[e4SL] += WSP[e1L]*dp*model.p_tv;
				WT[e4VG] += WSP[e1G]*dp*model.p_ts;
				WT[e4VL] += WSP[e1L]*dp*model.p_tv;
				
			}
			else if(nuc_a == nuc_d)
			{
			}
			else if(g_pupy[nuc_a] == g_pupy[nuc_d])
			{
				WT[e4SG] += WSP[e1G]*dp;
				WT[e4SL] += WSP[e1L]*dp;
			}
			else
			{
				WT[e4VG] += WSP[e1G]*dp;
				WT[e4VL] += WSP[e1L]*dp;
			}
			for(size_t k = d; k > 0; --k)
			{
				dp = model.p_indel_size[k];
				ex_s &WSK = table_s(a,d-k);
				noalias(WT) += table_t(a,d-k)*dp;
				WT[e4SG] += WSK[e1S]*dp;
				WT[e4VG] += WSK[e1V]*dp;
				WT[e4SL] += WSK[e1S]*dp*cache_log[k];
				WT[e4VL] += WSK[e1V]*dp*cache_log[k];
				WT[e4GL] += WSK[e1G]*dp*cache_log[k] + WSK[e1L]*dp + WSK[e1W]*dp*cache_log[k];
			}
			for(size_t k = a; k > 0; --k)
			{
				dp = model.p_indel_size[k];
				ex_s &WSK = table_s(a-k,d);
				noalias(WT) += table_t(a-k,d)*dp;
				WT[e4SG] += WSK[e1S]*dp;
				WT[e4VG] += WSK[e1V]*dp;
				WT[e4SL] += WSK[e1S]*dp*cache_log[k];
				WT[e4VL] += WSK[e1V]*dp*cache_log[k];
				WT[e4GL] += WSK[e1G]*dp*cache_log[k] + WSK[e1L]*dp + WSK[e1W]*dp*cache_log[k];
			}
		}
	}
	ex2 = table_t(sz_anc,sz_dec);
	
	ex2 /= w;
	ex[eSG] = ex2[e4SG];
	ex[eSL] = ex2[e4SL];
	ex[eVG] = ex2[e4VG];
	ex[eVL] = ex2[e4VL];
	ex[eGL] = ex2[e4GL];
	
	output() << "." << std::flush;
	
	return log(w) + log(model.p_end) - log(prob_scale)
		+ static_cast<double>(sz_anc+sz_dec)*log(model.nuc_scale)
		+ static_cast<double>(get_seq_info()[index][nN])*log(model.amb_scale);
	
}

void covar_k2p_zeta::covariance(covar_type &v, const params_type &p, const exvec_type &ev)
{
	ex_type ex;
	ex.clear();
	
	for(exvec_type::const_iterator cit=ev.begin();cit != ev.end();++cit)
		add_expectation(ex, *cit);
	
	double W = get_seqs().size();
	double M = ex[eM], S = ex[eS], V = ex[eV], G = ex[eG], L = ex[eL];
	double MM = ex[eMM], SS = ex[eSS], VV = ex[eVV], GG = ex[eGG], LL = ex[eLL];
	double MS = ex[eMS], MV = ex[eMV], MG = ex[eMG], ML = ex[eML], SV = ex[eSV];
	double SG = ex[eSG], SL = ex[eSL], VG = ex[eVG], VL = ex[eVL], GL = ex[eGL];
		
	double t = p[pT], k = p[pK], r = p[pR],
		z = p[pZ], a = p[pA];
	
	double ek = 1.0+k;
	double ek2 = 1.0 + 2.0*k;
	double ea = exp(-t*ek2/ek);
	double eb = exp(-2.0*t/ek);
	double er = exp(-2.0*t*r);
	
	double da = 4.0/(1.0 - 2.0*ea + eb);
	double db = 4.0/(1.0 + 2.0*ea + eb);
	double dc = eb-1.0;
	double de = -2.0*eb/(dc*ek);
	double df = er/(1.0-er);
	
	double tt = ((8.0*de*(V + 2.0*ek*r*(MV + SV - df*VG)))/ek + 4.0*V*sq(de) +
		(2.0*da*(-2.0*eb*(S + 2.0*ek*r*(MS - df*SG + SS)) + 2.0*eb*ek*(de - 2.0*r)*SV +
		ea*ek2*(4.0*ek*MS*r + S + 2.0*k*S -
		2.0*ek*(2.0*df*r*SG + de*SV - 2.0*r*(SS + SV)))))/sq(ek) -
		(2.0*db*(ea*ek2*(M + 2.0*k*M +
		2.0*ek*(-(de*MV) + 2.0*(-(df*MG) + MM + MS + MV)*r)) +
		eb*(2.0*M + da*eb*MS +
		2.0*ek*(-(de*MV) + 2.0*(-(df*MG) + MM + MS + MV)*r)) -
		da*MS*sq(ea + 2.0*ea*k)))/sq(ek) +
		((S - SS)*sq(da)*sq(ea - eb + 2.0*ea*k))/sq(ek) +
		((M - MM)*sq(db)*sq(ea + eb + 2.0*ea*k))/sq(ek) +
		16.0*(-MM - 2.0*MS - 2.0*MV - SS - 2.0*SV +
		df*(G + df*G - df*GG + 2.0*(MG + SG + VG)))*sq(r) -
		(4.0*VV*sq(de)*sq(eb + dc*ek*r))/sq(eb))/4.0
		;
	
	double tk = (-2.0*da*dc*(dc*(ea + eb)*ek*S +
		(-(dc*ea*(2.0*ek*MS*r + S + 2.0*k*S + 2.0*ek*r*(-(df*SG) + SS))) -
		2.0*dc*eb*(S + ek*r*(MS - df*SG + SS)) +
		eb*(-2.0*eb + dc*de*ek - 2.0*dc*ek*r)*SV +
		ea*(dc*de*ek + 2.0*eb*ek2 - 2.0*dc*ek*r)*SV)*t) +
		8.0*eb*(2.0*dc*ek*MV*r*t - (dc*ek + 2.0*t)*V + 2.0*eb*t*VV +
		2.0*dc*ek*r*t*(SV - df*VG + VV)) +
		(ea + eb)*(ea - eb + 2.0*ea*k)*(S - SS)*t*sq(da)*sq(dc) +
		(ea - eb)*(ea + eb + 2.0*ea*k)*(M - MM)*t*sq(db)*sq(dc) +
		2.0*db*dc*(eb*(-(dc*M*(ek - 2.0*t)) +
		(da*dc*eb*MS + 2.0*eb*MV - dc*de*ek*MV +
		2.0*dc*ek*(-(df*MG) + MM + MS + MV)*r)*t) +
		ea*((dc*de*ek*MV + 2.0*eb*ek2*MV -
		2.0*dc*ek*(-(df*MG) + MM + MS + MV)*r)*t + dc*M*(ek - t - 2.0*k*t))
		+ da*dc*ek2*MS*t*sq(ea)))/(4.0*sq(ek)*ek*sq(dc))
		;

	
	double tr = (-2.0*k*M - 2.0*k*S + db*ea*MM*t + 2.0*db*ea*k*MM*t - da*ea*MS*t + db*ea*MS*t -
		2.0*da*ea*k*MS*t + 2.0*db*ea*k*MS*t - 2.0*de*MV*t + db*ea*MV*t -
		2.0*de*k*MV*t + 2.0*db*ea*k*MV*t + 4.0*MM*r*t + 4.0*k*MM*r*t + 8.0*MS*r*t +
		8.0*k*MS*r*t + 8.0*MV*r*t + 8.0*k*MV*r*t - da*ea*SS*t - 2.0*da*ea*k*SS*t +
		4.0*r*SS*t + 4.0*k*r*SS*t - 2.0*de*SV*t - da*ea*SV*t - 2.0*de*k*SV*t -
		2.0*da*ea*k*SV*t + 8.0*r*SV*t + 8.0*k*r*SV*t - 2.0*k*V - 2.0*(M + S + V) +
		dc*df*(2.0*ek*G*(-1.0 + 2.0*r*t) +
		t*(db*(ea + eb + 2.0*ea*k)*MG + 8.0*ek*MG*r -
		(da*(ea - eb + 2.0*ea*k) - 8.0*ek*r)*SG - 2.0*ek*(de - 4.0*r)*VG)) +
		4.0*ek*r*t*VV + eb*(2.0*ek*M + 2.0*S - db*(-1.0 + ea)*MM*t + 2.0*V +
		2.0*k*(S - db*ea*(MM + MS + MV)*t + V +
		t*(de*(MV + SV) + da*ea*(MS + SS + SV) -
		2.0*r*(MM + 2.0*MS + 2.0*MV + SS + 2.0*SV + VV))) +
		t*(-(db*(-1.0 + ea)*(MS + MV)) + da*(1.0 + ea)*(MS + SS + SV) +
		2.0*(de*(MV + SV) - 2.0*VV - 2.0*r*(MM + 2.0*MS + 2.0*MV + SS + 2.0*SV + VV)))
		) + 4.0*dc*ek*(G - GG)*r*t*sq(df) -
		(db*(MM + MS + MV) + da*(MS + SS + SV))*t*sq(eb))/(dc*ek)
		;

	double tz = -(db*(ea + eb + 2.0*ea*k)*ML - 4.0*df*ek*GL*r + 4.0*ML*r + 4.0*k*ML*r - da*ea*SL +
		da*eb*SL - 2.0*da*ea*k*SL + 4.0*r*SL + 4.0*k*r*SL - 2.0*ek*(de - 2.0*r)*VL +
		(db*(ea + eb + 2.0*ea*k)*MG - 4.0*df*ek*GG*r + 4.0*MG*r + 4.0*k*MG*r -
		da*ea*SG + da*eb*SG - 2.0*da*ea*k*SG + 4.0*r*SG + 4.0*k*r*SG -
		2.0*ek*(de - 2.0*r)*VG)*zeta_ratio(z))/(2.0*ek)
		;
	
	double ta = (db*(ea + eb + 2.0*ea*k)*(MG + MM + MS + MV - a*M*W) -
		da*(ea - eb + 2.0*ea*k)*(MS + SG + SS + SV - a*S*W) +
		2.0*ek*(-(de*(MV + SV + VG + VV - a*V*W)) +
		2.0*r*(MG + MM + 2.0*MS + 2.0*MV + SG + SS + 2.0*SV + VG -
		df*(GG + MG + SG + VG) + VV - a*(-(df*G) + M + S + V)*W)))/
		(2.0*a*(1.0 + a)*ek)
		;

	double kk = (t*(2.0*da*dc*(2.0*dc*(ea + eb)*ek*S +
		(dc*(ea - 2.0*eb)*S - 4.0*eb*(ea + eb)*SV)*t) +
		16.0*eb*((dc*ek + t)*V - eb*t*VV) +
		2.0*db*dc*(-2.0*dc*(ea - eb)*ek*M +
		t*(4.0*(ea - eb)*eb*MV -
		dc*(ea*M + eb*(2.0*M + da*eb*MS) - da*MS*sq(ea)))) +
		(M - MM)*t*sq(db)*sq(dc)*sq(ea - eb) +
		(S - SS)*t*sq(da)*sq(dc)*sq(ea + eb)))/(4.0*sq(ek*ek*dc))
		;
	
	double kr = ((db*dc*(ea - eb)*(df*MG - MM - MS - MV) +
		da*dc*(ea + eb)*(MS - df*SG + SS + SV) + 4.0*eb*(MV + SV - df*VG + VV))*
		sq(t))/(dc*sq(ek))
		;
	
   double kz = (t*(dc*(db*(-ea + eb)*ML + da*(ea + eb)*SL) + 4.0*eb*VL +
		(dc*(db*(-ea + eb)*MG + da*(ea + eb)*SG) + 4.0*eb*VG)*zeta_ratio(z)))/
		(2.0*dc*sq(ek))
		;
	
	double ka = (t*(db*dc*(ea - eb)*(MG + MM + MS + MV - a*M*W) -
		da*dc*(ea + eb)*(MS + SG + SS + SV - a*S*W) -
		4.0*eb*(MV + SV + VG + VV - a*V*W)))/(2.0*a*(1.0 + a)*dc*sq(ek))
		;

   double rr = 4.0*(-MM - 2.0*MS - 2.0*MV - SS - 2.0*SV +
		df*(G + df*G - df*GG + 2.0*(MG + SG + VG)) - VV)*sq(t)
		;

   double rz = (2.0*t*(ML + SL + VL - er*(GL + ML + SL + VL) +
		(MG + SG + VG - er*(GG + MG + SG + VG))*zeta_ratio(z)))/(-1.0 + er)
		;
	
	double ra = (2.0*t*(-(df*(-1.0 + er)*(GG + MG + SG + VG)) +
		(-1.0 + er)*(MG + MM + 2.0*MS + 2.0*MV + SG + SS + 2.0*SV + VG + VV) -
		a*(-M - S - V + er*(G + M + S + V))*W))/(a*(1.0 + a)*(-1.0 + er))
		;


	double zz = -LL - (G + GG)*sq(zeta_ratio(z)) - 2.0*GL*zeta_ratio(z) + G*zeta_ratio2(z);
	
	double za = (GL + ML + SL + VL - a*L*W + (GG + MG + SG + VG - a*G*W)*zeta_ratio(z))/
		(a*(1.0 + a))
		;
	
	double aa = -((GG - M + 2.0*MG + MM + 2.0*MS + 2.0*MV - S + 2.0*SG + SS + 2.0*SV - V + 2.0*VG +
		VV + a*(1.0 + W)*(-2.0*(M + S + V) + a*W) - G*(1.0 + 2.0*a*(1.0 + W)))/
		(sq(a)*sq(1.0 + a)))
		;
	
	covar_type fim(p.size(), p.size());
	fim(pT,pT) = tt; fim(pT,pK) = tk; fim(pT,pR) = tr; fim(pT,pZ) = tz; fim(pT,pA) = ta;
	fim(pK,pT) = tk; fim(pK,pK) = kk; fim(pK,pR) = kr; fim(pK,pZ) = kz; fim(pK,pA) = ka;
	fim(pR,pT) = tr; fim(pR,pK) = kr; fim(pR,pR) = rr; fim(pR,pZ) = rz; fim(pR,pA) = ra;
	fim(pZ,pT) = tz; fim(pZ,pK) = kz; fim(pZ,pR) = rz; fim(pZ,pZ) = zz; fim(pZ,pA) = za;
	fim(pA,pT) = ta; fim(pA,pK) = ka; fim(pA,pR) = ra; fim(pA,pZ) = za; fim(pA,pA) = aa;
		
	invert_matrix(fim, v);
}


void covar_k2p_zeta::add_expectation(ex_type &ex, const ex_type &ad) const
{
	cout << ad << endl;
	
	ex[eMM] += ad[eMM] + 2.0*ex[eM]*ad[eM];
	ex[eSS] += ad[eSS] + 2.0*ex[eS]*ad[eS];
	ex[eVV] += ad[eVV] + 2.0*ex[eV]*ad[eV];
	ex[eGG] += ad[eGG] + 2.0*ex[eG]*ad[eG];
	ex[eLL] += ad[eLL] + 2.0*ex[eL]*ad[eL];
	
	ex[eMS] += ad[eMS] + ex[eM]*ad[eS] + ex[eS]*ad[eM];
	ex[eMV] += ad[eMV] + ex[eM]*ad[eV] + ex[eV]*ad[eM];
	ex[eMG] += ad[eMG] + ex[eM]*ad[eG] + ex[eG]*ad[eM];
	ex[eML] += ad[eML] + ex[eM]*ad[eL] + ex[eL]*ad[eM];
	ex[eSV] += ad[eSV] + ex[eS]*ad[eV] + ex[eV]*ad[eS];
	ex[eSG] += ad[eSG] + ex[eS]*ad[eG] + ex[eG]*ad[eS];
	ex[eSL] += ad[eSL] + ex[eS]*ad[eL] + ex[eL]*ad[eS];
	ex[eVG] += ad[eVG] + ex[eV]*ad[eG] + ex[eG]*ad[eV];
	ex[eVL] += ad[eVL] + ex[eV]*ad[eL] + ex[eL]*ad[eV];
	ex[eGL] += ad[eGL] + ex[eG]*ad[eL] + ex[eL]*ad[eG];

	ex[eM] += ad[eM];
	ex[eS] += ad[eS];
	ex[eV] += ad[eV];
	ex[eG] += ad[eG];
	ex[eL] += ad[eL];
}


std::string covar_k2p_zeta::state() const
{
	ex_type ex;
	ex.clear();
	const exvec_type &ev = get_exvec();
	for(exvec_type::const_iterator cit=ev.begin();cit != ev.end();++cit)
		add_expectation(ex, *cit);	
	
	const covar_type &v = get_covar();
	
	ostringstream ostr;
	ostr << std::setprecision(DBL_DIG);
	ostr << "ex = " << ex << endl;
	ostr << "var = " << v << endl;
	return ostr.str();
}

/***************************************************************************
 * class covar_k2p_geo                                                     *
 ***************************************************************************/

void covar_k2p_geo::preallocate(size_t maxa, size_t maxd)
{
	sz_height = maxa+1;
	sz_width = maxd+1;
		
	size_t sz_max = std::max(sz_height, sz_width);
	
	cache_size.resize(sz_max, 0.0);
	for(size_t u=1;u<sz_max;++u)
		cache_size[u] = static_cast<double>(u);
}

void covar_k2p_geo::expectation_setup(const params_type &params)
{

}

double covar_k2p_geo::expectation(ex_type &ex, const sequence &seq_a, const sequence &seq_d, size_t index) const
{
	size_t sz_anc = seq_a.size();
	size_t sz_dec = seq_d.size();

	const model_type &model = get_model();
	
	enum {xW, xE, xF, xEE, xEF = xEE};

	typedef ublas::cc_vector<double, 4> ex_s;
	ex_s szero = ublas::zero_vector<double>(4);
	table<ex_s> table_s(sz_anc+1,sz_dec+1,szero);
	ex_s row_cache;
	vector<ex_s> col_cache(sz_dec+1,szero);

	table_s(0,0)[xW] = prob_scale;
	double dp = model.p_open;
	int nuc_a, nuc_d;
	// M, MM
	for(size_t d = 1; d <= sz_dec; ++d)
	{
		ex_s &WS = table_s(0,d);
		ex_s &WSK = table_s(0,d-1);
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		noalias(WS) += row_cache*dp;
	}
	row_cache.clear();
	for(size_t a = 1; a <= sz_anc; ++a)
	{
		ex_s &WS = table_s(a,0);
		ex_s &WSK = table_s(a-1,0);
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		noalias(WS) += row_cache*dp;
	}
	for(size_t a=1;a<=sz_anc;++a)
	{
		row_cache.clear();
		nuc_a = seq_a[a-1];
		for(size_t d=1;d<=sz_dec;++d)
		{
			nuc_d = seq_d[d-1];
			dp = model.p_substitution[nuc_a][nuc_d];
				
			ex_s &WS = table_s(a,d);
			ex_s &WSP = table_s(a-1,d-1);
			noalias(WS) = WSP*dp;
			if(nuc_a == nN || nuc_d == nN)
			{
				WS[xE] += WS[xW]*model.p_match;
				WS[xEE] += (WSP[xE]+WSP[xE]+WSP[xW]*model.p_match)*model.p_match*dp;
			}
			else if(nuc_a == nuc_d)
			{
				WS[xE] += WS[xW];
				WS[xEE] += (WSP[xE]+WSP[xE]+WSP[xW])*dp;
			}
			
			dp = model.p_open;
			{
				ex_s &WSK = table_s(a,d-1);
				row_cache *= model.p_extend;
				noalias(row_cache) += WSK;
				noalias(WS) += row_cache*dp;
			}
			
			{
				ex_s &WSK = table_s(a-1,d);
				col_cache[d] *= model.p_extend;
				noalias(col_cache[d]) += WSK;
				noalias(WS) += col_cache[d]*dp;
			}
		}
	}
	ex_s ex1 = table_s(sz_anc,sz_dec);
	double w = ex1[xW];
	ex1 /= w;
	ex[eM] = ex1[xE]; ex[eMM] = ex1[xEE];
	output() << "." << std::flush;

	// S, SS
	row_cache.clear();
	fill(col_cache.begin(), col_cache.end(), szero);
	table_s(0,0)[xW] = prob_scale;
	dp = model.p_open;
	for(size_t d = 1; d <= sz_dec; ++d)
	{
		ex_s &WS = table_s(0,d);
		ex_s &WSK = table_s(0,d-1);
		WS.clear();
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		noalias(WS) += row_cache*dp;
	}
	row_cache.clear();
	for(size_t a = 1; a <= sz_anc; ++a)
	{
		ex_s &WS = table_s(a,0);
		ex_s &WSK = table_s(a-1,0);
		WS.clear();
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		noalias(WS) += row_cache*dp;
	}
	for(size_t a=1;a<=sz_anc;++a)
	{
		row_cache.clear();
		nuc_a = seq_a[a-1];
		for(size_t d=1;d<=sz_dec;++d)
		{
			nuc_d = seq_d[d-1];
			dp = model.p_substitution[nuc_a][nuc_d];
				
			ex_s &WS = table_s(a,d);
			ex_s &WSP = table_s(a-1,d-1);
			noalias(WS) = WSP*dp;
			if(nuc_a == nN || nuc_d == nN)
			{
				WS[xE] += WS[xW]*model.p_ts;
				WS[xEE] += (WSP[xE]+WSP[xE]+WSP[xW]*model.p_ts)*model.p_ts*dp;
			}
			else if(nuc_a != nuc_d && g_pupy[nuc_a] == g_pupy[nuc_d])
			{
				WS[xE] += WS[xW];
				WS[xEE] += (WSP[xE]+WSP[xE]+WSP[xW])*dp;
			}
			
			dp = model.p_open;
			{
				ex_s &WSK = table_s(a,d-1);
				row_cache *= model.p_extend;
				noalias(row_cache) += WSK;
				noalias(WS) += row_cache*dp;
			}
			
			{
				ex_s &WSK = table_s(a-1,d);
				col_cache[d] *= model.p_extend;
				noalias(col_cache[d]) += WSK;
				noalias(WS) += col_cache[d]*dp;
			}
		}
	}
	ex1 = table_s(sz_anc,sz_dec);
	ex1 /= w;
	ex[eS] = ex1[xE]; ex[eSS] = ex1[xEE];
	output() << "." << std::flush;

	// V, VV
	row_cache.clear();
	fill(col_cache.begin(), col_cache.end(), szero);
	table_s(0,0)[xW] = prob_scale;
	dp = model.p_open;
	for(size_t d = 1; d <= sz_dec; ++d)
	{
		ex_s &WS = table_s(0,d);
		ex_s &WSK = table_s(0,d-1);
		WS.clear();
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		noalias(WS) += row_cache*dp;
	}
	row_cache.clear();	
	for(size_t a = 1; a <= sz_anc; ++a)
	{
		ex_s &WS = table_s(a,0);
		ex_s &WSK = table_s(a-1,0);
		WS.clear();
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		noalias(WS) += row_cache*dp;
	}
	for(size_t a=1;a<=sz_anc;++a)
	{
		row_cache.clear();
		nuc_a = seq_a[a-1];
		for(size_t d=1;d<=sz_dec;++d)
		{
			nuc_d = seq_d[d-1];
			dp = model.p_substitution[nuc_a][nuc_d];
				
			ex_s &WS = table_s(a,d);
			ex_s &WSP = table_s(a-1,d-1);
			noalias(WS) = WSP*dp;
			if(nuc_a == nN || nuc_d == nN)
			{
				WS[xE] += WS[xW]*model.p_tv;
				WS[xEE] += (WSP[xE]+WSP[xE]+WSP[xW]*model.p_tv)*model.p_tv*dp;
			}
			else if(nuc_a != nuc_d && g_pupy[nuc_a] != g_pupy[nuc_d])
			{
				WS[xE] += WS[xW];
				WS[xEE] += (WSP[xE]+WSP[xE]+WSP[xW])*dp;
			}
			
			dp = model.p_open;
			{
				ex_s &WSK = table_s(a,d-1);
				row_cache *= model.p_extend;
				noalias(row_cache) += WSK;
				noalias(WS) += row_cache*dp;
			}
			
			{
				ex_s &WSK = table_s(a-1,d);
				col_cache[d] *= model.p_extend;
				noalias(col_cache[d]) += WSK;
				noalias(WS) += col_cache[d]*dp;
			}
		}
	}
	ex1 = table_s(sz_anc,sz_dec);
	ex1 /= w;
	ex[eV] = ex1[xE]; ex[eVV] = ex1[xEE];
	output() << "." << std::flush;
	
	// G, GG
	row_cache.clear();
	fill(col_cache.begin(), col_cache.end(), szero);
	table_s(0,0)[xW] = prob_scale;
	dp = model.p_open;
	for(size_t d = 1; d <= sz_dec; ++d)
	{
		ex_s &WS = table_s(0,d);
		ex_s &WSK = table_s(0,d-1);
		WS.clear();
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		row_cache[xE] += WSK[xW];
		row_cache[xEE] += WSK[xE]+WSK[xE]+WSK[xW];
		noalias(WS) += row_cache*dp;
	}
	row_cache.clear();
	for(size_t a = 1; a <= sz_anc; ++a)
	{
		ex_s &WS = table_s(a,0);
		ex_s &WSK = table_s(a-1,0);
		WS.clear();
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		row_cache[xE] += WSK[xW];
		row_cache[xEE] += WSK[xE]+WSK[xE]+WSK[xW];
		noalias(WS) += row_cache*dp;
	}
	for(size_t a=1;a<=sz_anc;++a)
	{
		row_cache.clear();
		nuc_a = seq_a[a-1];
		for(size_t d=1;d<=sz_dec;++d)
		{
			nuc_d = seq_d[d-1];
			ex_s &WS = table_s(a,d);
			ex_s &WSP = table_s(a-1,d-1);
			dp = model.p_substitution[nuc_a][nuc_d];
			noalias(WS) = WSP*dp;
			dp = model.p_open;
			{
				ex_s &WSK = table_s(a,d-1);
				row_cache *= model.p_extend;
				noalias(row_cache) += WSK;
				row_cache[xE] += WSK[xW];
				row_cache[xEE] += WSK[xE]+WSK[xE]+WSK[xW];
				noalias(WS) += row_cache*dp;
			}
			
			{
				ex_s &WSK = table_s(a-1,d);
				col_cache[d] *= model.p_extend;
				noalias(col_cache[d]) += WSK;
				col_cache[d][xE] += WSK[xW];
				col_cache[d][xEE] += WSK[xE]+WSK[xE]+WSK[xW];
				noalias(WS) += col_cache[d]*dp;
			}
		}
	}
	ex1 = table_s(sz_anc,sz_dec);
	ex1 /= w;
	ex[eG] = ex1[xE]; ex[eGG] = ex1[xEE];
	output() << "." << std::flush;

	// Y, YY
	row_cache.clear();
	fill(col_cache.begin(), col_cache.end(), szero);
	table_s(0,0)[xW] = prob_scale;
	dp = model.p_open;
	for(size_t d = 1; d <= sz_dec; ++d)
	{
		ex_s &WS = table_s(0,d);
		ex_s &WSK = table_s(0,d-1);
		WS.clear();
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		row_cache[xE] += row_cache[xW];
		row_cache[xEE] += WSK[xE]+WSK[xE]+row_cache[xW];
		noalias(WS) += row_cache*dp;
	}
	row_cache.clear();	
	for(size_t a = 1; a <= sz_anc; ++a)
	{
		ex_s &WS = table_s(a,0);
		ex_s &WSK = table_s(a-1,0);
		WS.clear();
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		row_cache[xE] += row_cache[xW];
		row_cache[xEE] += WSK[xE]+WSK[xE]+row_cache[xW];
		noalias(WS) += row_cache*dp;
	}
	for(size_t a=1;a<=sz_anc;++a)
	{
		row_cache.clear();
		nuc_a = seq_a[a-1];
		for(size_t d=1;d<=sz_dec;++d)
		{
			nuc_d = seq_d[d-1];
			ex_s &WS = table_s(a,d);
			ex_s &WSP = table_s(a-1,d-1);
			dp = model.p_substitution[nuc_a][nuc_d];
			noalias(WS) = WSP*dp;
			dp = model.p_open;
			{
				ex_s &WSK = table_s(a,d-1);
				row_cache *= model.p_extend;
				noalias(row_cache) += WSK;
				row_cache[xE] += row_cache[xW];
				row_cache[xEE] += WSK[xE]+WSK[xE]+row_cache[xW];
				noalias(WS) += row_cache*dp;
			}
			
			{
				ex_s &WSK = table_s(a-1,d);
				col_cache[d] *= model.p_extend;
				noalias(col_cache[d]) += WSK;
				col_cache[d][xE] += col_cache[d][xW];
				col_cache[d][xEE] += WSK[xE]+WSK[xE]+col_cache[d][xW];
				noalias(WS) += col_cache[d]*dp;
			}
		}
	}
	ex1 = table_s(sz_anc,sz_dec);
	ex1 /= w;
	ex[eY] = ex1[xE]; ex[eYY] = ex1[xEE];
	output() << "." << std::flush;
	
	// MS
	row_cache.clear();
	fill(col_cache.begin(), col_cache.end(), szero);
	table_s(0,0)[xW] = prob_scale;
	dp = model.p_open;
	for(size_t d = 1; d <= sz_dec; ++d)
	{
		ex_s &WS = table_s(0,d);
		ex_s &WSK = table_s(0,d-1);
		WS.clear();
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		noalias(WS) += row_cache*dp;
	}
	row_cache.clear();
	for(size_t a = 1; a <= sz_anc; ++a)
	{
		ex_s &WS = table_s(a,0);
		ex_s &WSK = table_s(a-1,0);
		WS.clear();
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		noalias(WS) += row_cache*dp;
	}
	for(size_t a=1;a<=sz_anc;++a)
	{
		row_cache.clear();
		nuc_a = seq_a[a-1];
		for(size_t d=1;d<=sz_dec;++d)
		{
			nuc_d = seq_d[d-1];
			dp = model.p_substitution[nuc_a][nuc_d];
				
			ex_s &WS = table_s(a,d);
			ex_s &WSP = table_s(a-1,d-1);
			noalias(WS) = WSP*dp;
			if(nuc_a == nN || nuc_d == nN)
			{
				WS[xE] += WS[xW]*model.p_match;
				WS[xF] += WS[xW]*model.p_ts;
				WS[xEF] += (WSP[xE]*model.p_ts+WSP[xF]*model.p_match+WSP[xW]*model.p_match*model.p_ts)*dp;
			}
			else if(nuc_a == nuc_d)
			{
				WS[xE] += WS[xW];
				WS[xEF] += WSP[xF]*dp;
			}			
			else if(g_pupy[nuc_a] == g_pupy[nuc_d])
			{
				WS[xF] += WS[xW];
				WS[xEF] += WSP[xE]*dp;
			}
			
			dp = model.p_open;
			{
				ex_s &WSK = table_s(a,d-1);
				row_cache *= model.p_extend;
				noalias(row_cache) += WSK;
				noalias(WS) += row_cache*dp;
			}
			
			{
				ex_s &WSK = table_s(a-1,d);
				col_cache[d] *= model.p_extend;
				noalias(col_cache[d]) += WSK;
				noalias(WS) += col_cache[d]*dp;
			}
		}
	}
	ex1 = table_s(sz_anc,sz_dec);
	ex1 /= w;
	ex[eMS] = ex1[xEF];
	output() << "." << std::flush;
	
	// MV
	row_cache.clear();
	fill(col_cache.begin(), col_cache.end(), szero);
	table_s(0,0)[xW] = prob_scale;
	dp = model.p_open;
	for(size_t d = 1; d <= sz_dec; ++d)
	{
		ex_s &WS = table_s(0,d);
		ex_s &WSK = table_s(0,d-1);
		WS.clear();
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		noalias(WS) += row_cache*dp;
	}
	row_cache.clear();
	for(size_t a = 1; a <= sz_anc; ++a)
	{
		ex_s &WS = table_s(a,0);
		ex_s &WSK = table_s(a-1,0);
		WS.clear();
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		noalias(WS) += row_cache*dp;
	}	
	for(size_t a=1;a<=sz_anc;++a)
	{
		row_cache.clear();
		nuc_a = seq_a[a-1];
		for(size_t d=1;d<=sz_dec;++d)
		{
			nuc_d = seq_d[d-1];
			dp = model.p_substitution[nuc_a][nuc_d];
				
			ex_s &WS = table_s(a,d);
			ex_s &WSP = table_s(a-1,d-1);
			noalias(WS) = WSP*dp;
			if(nuc_a == nN || nuc_d == nN)
			{
				WS[xE] += WS[xW]*model.p_match;
				WS[xF] += WS[xW]*model.p_tv;
				WS[xEF] += (WSP[xE]*model.p_tv+WSP[xF]*model.p_match+WSP[xW]*model.p_match*model.p_tv)*dp;
			}
			else if(nuc_a == nuc_d)
			{
				WS[xE] += WS[xW];
				WS[xEF] += WSP[xF]*dp;
			}			
			else if(g_pupy[nuc_a] != g_pupy[nuc_d])
			{
				WS[xF] += WS[xW];
				WS[xEF] += WSP[xE]*dp;
			}
			
			dp = model.p_open;
			{
				ex_s &WSK = table_s(a,d-1);
				row_cache *= model.p_extend;
				noalias(row_cache) += WSK;
				noalias(WS) += row_cache*dp;
			}
			
			{
				ex_s &WSK = table_s(a-1,d);
				col_cache[d] *= model.p_extend;
				noalias(col_cache[d]) += WSK;
				noalias(WS) += col_cache[d]*dp;
			}
		}
	}
	ex1 = table_s(sz_anc,sz_dec);
	ex1 /= w;
	ex[eMV] = ex1[xEF];
	output() << "." << std::flush;

	// SV
	row_cache.clear();
	fill(col_cache.begin(), col_cache.end(), szero);
	table_s(0,0)[xW] = prob_scale;
	dp = model.p_open;
	for(size_t d = 1; d <= sz_dec; ++d)
	{
		ex_s &WS = table_s(0,d);
		ex_s &WSK = table_s(0,d-1);
		WS.clear();
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		noalias(WS) += row_cache*dp;
	}
	row_cache.clear();
	for(size_t a = 1; a <= sz_anc; ++a)
	{
		ex_s &WS = table_s(a,0);
		ex_s &WSK = table_s(a-1,0);
		WS.clear();
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		noalias(WS) += row_cache*dp;
	}
	for(size_t a=1;a<=sz_anc;++a)
	{
		row_cache.clear();
		nuc_a = seq_a[a-1];
		for(size_t d=1;d<=sz_dec;++d)
		{
			nuc_d = seq_d[d-1];
			dp = model.p_substitution[nuc_a][nuc_d];
				
			ex_s &WS = table_s(a,d);
			ex_s &WSP = table_s(a-1,d-1);
			noalias(WS) = WSP*dp;
			if(nuc_a == nN || nuc_d == nN)
			{
				WS[xE] += WS[xW]*model.p_ts;
				WS[xF] += WS[xW]*model.p_tv;
				WS[xEF] += (WSP[xE]*model.p_tv+WSP[xF]*model.p_ts+WSP[xW]*model.p_ts*model.p_tv)*dp;
			}
			else if(nuc_a == nuc_d)
				;
			else if(g_pupy[nuc_a] == g_pupy[nuc_d])
			{
				WS[xE] += WS[xW];
				WS[xEF] += WSP[xF]*dp;
			}
			else
			{
				WS[xF] += WS[xW];
				WS[xEF] += WSP[xE]*dp;
			}
			
			dp = model.p_open;
			{
				ex_s &WSK = table_s(a,d-1);
				row_cache *= model.p_extend;
				noalias(row_cache) += WSK;
				noalias(WS) += row_cache*dp;
			}
			
			{
				ex_s &WSK = table_s(a-1,d);
				col_cache[d] *= model.p_extend;
				noalias(col_cache[d]) += WSK;
				noalias(WS) += col_cache[d]*dp;
			}
		}
	}
	ex1 = table_s(sz_anc,sz_dec);
	ex1 /= w;
	ex[eSV] = ex1[xEF];
	output() << "." << std::flush;
	
	// MG
	row_cache.clear();
	fill(col_cache.begin(), col_cache.end(), szero);
	table_s(0,0)[xW] = prob_scale;
	dp = model.p_open;
	for(size_t d = 1; d <= sz_dec; ++d)
	{
		ex_s &WS = table_s(0,d);
		ex_s &WSK = table_s(0,d-1);
		WS.clear();
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		row_cache[xF] += WSK[xW];
		row_cache[xEF] += WSK[xE];
		noalias(WS) += row_cache*dp;
	}
	row_cache.clear();
	for(size_t a = 1; a <= sz_anc; ++a)
	{
		ex_s &WS = table_s(a,0);
		ex_s &WSK = table_s(a-1,0);
		WS.clear();
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		row_cache[xF] += WSK[xW];
		row_cache[xEF] += WSK[xE];
		noalias(WS) += row_cache*dp;
	}
	for(size_t a=1;a<=sz_anc;++a)
	{
		row_cache.clear();
		nuc_a = seq_a[a-1];
		for(size_t d=1;d<=sz_dec;++d)
		{
			nuc_d = seq_d[d-1];
			dp = model.p_substitution[nuc_a][nuc_d];
				
			ex_s &WS = table_s(a,d);
			ex_s &WSP = table_s(a-1,d-1);
			noalias(WS) = WSP*dp;
			if(nuc_a == nN || nuc_d == nN)
			{
				WS[xE] += WS[xW]*model.p_match;
				WS[xEF] += WSP[xF]*model.p_match*dp;
			}
			else if(nuc_a == nuc_d)
			{
				WS[xE] += WS[xW];
				WS[xEF] += WSP[xF]*dp;
			}
			
			dp = model.p_open;
			{
				ex_s &WSK = table_s(a,d-1);
				row_cache *= model.p_extend;
				noalias(row_cache) += WSK;
				row_cache[xF] += WSK[xW];
				row_cache[xEF] += WSK[xE];
				noalias(WS) += row_cache*dp;
			}
			
			{
				ex_s &WSK = table_s(a-1,d);
				col_cache[d] *= model.p_extend;
				noalias(col_cache[d]) += WSK;
				col_cache[d][xF] += WSK[xW];
				col_cache[d][xEF] += WSK[xE];
				noalias(WS) += col_cache[d]*dp;
			}
		}
	}
	ex1 = table_s(sz_anc,sz_dec);
	ex1 /= w;
	ex[eMG] = ex1[xEF];
	output() << "." << std::flush;
	
	// SG
	row_cache.clear();
	fill(col_cache.begin(), col_cache.end(), szero);
	table_s(0,0)[xW] = prob_scale;
	dp = model.p_open;
	for(size_t d = 1; d <= sz_dec; ++d)
	{
		ex_s &WS = table_s(0,d);
		ex_s &WSK = table_s(0,d-1);
		WS.clear();
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		row_cache[xF] += WSK[xW];
		row_cache[xEF] += WSK[xE];
		noalias(WS) += row_cache*dp;
	}
	row_cache.clear();
	for(size_t a = 1; a <= sz_anc; ++a)
	{
		ex_s &WS = table_s(a,0);
		ex_s &WSK = table_s(a-1,0);
		WS.clear();
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		row_cache[xF] += WSK[xW];
		row_cache[xEF] += WSK[xE];
		noalias(WS) += row_cache*dp;
	}
	for(size_t a=1;a<=sz_anc;++a)
	{
		row_cache.clear();
		nuc_a = seq_a[a-1];
		for(size_t d=1;d<=sz_dec;++d)
		{
			nuc_d = seq_d[d-1];
			dp = model.p_substitution[nuc_a][nuc_d];
				
			ex_s &WS = table_s(a,d);
			ex_s &WSP = table_s(a-1,d-1);
			noalias(WS) = WSP*dp;
			if(nuc_a == nN || nuc_d == nN)
			{
				WS[xE] += WS[xW]*model.p_ts;
				WS[xEF] += WSP[xF]*model.p_ts*dp;
			}
			else if(nuc_a != nuc_d && g_pupy[nuc_a] == g_pupy[nuc_d])
			{
				WS[xE] += WS[xW];
				WS[xEF] += WSP[xF]*dp;
			}
			
			dp = model.p_open;
			{
				ex_s &WSK = table_s(a,d-1);
				row_cache *= model.p_extend;
				noalias(row_cache) += WSK;
				row_cache[xF] += WSK[xW];
				row_cache[xEF] += WSK[xE];
				noalias(WS) += row_cache*dp;
			}
			
			{
				ex_s &WSK = table_s(a-1,d);
				col_cache[d] *= model.p_extend;
				noalias(col_cache[d]) += WSK;
				col_cache[d][xF] += WSK[xW];
				col_cache[d][xEF] += WSK[xE];
				noalias(WS) += col_cache[d]*dp;
			}
		}
	}
	ex1 = table_s(sz_anc,sz_dec);
	ex1 /= w;
	ex[eSG] = ex1[xEF];
	output() << "." << std::flush;
	
	// VG
	row_cache.clear();
	fill(col_cache.begin(), col_cache.end(), szero);
	table_s(0,0)[xW] = prob_scale;
	dp = model.p_open;
	for(size_t d = 1; d <= sz_dec; ++d)
	{
		ex_s &WS = table_s(0,d);
		ex_s &WSK = table_s(0,d-1);
		WS.clear();
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		row_cache[xF] += WSK[xW];
		row_cache[xEF] += WSK[xE];
		noalias(WS) += row_cache*dp;
	}
	row_cache.clear();
	for(size_t a = 1; a <= sz_anc; ++a)
	{
		ex_s &WS = table_s(a,0);
		ex_s &WSK = table_s(a-1,0);
		WS.clear();
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		row_cache[xF] += WSK[xW];
		row_cache[xEF] += WSK[xE];
		noalias(WS) += row_cache*dp;
	}
	for(size_t a=1;a<=sz_anc;++a)
	{
		row_cache.clear();
		nuc_a = seq_a[a-1];
		for(size_t d=1;d<=sz_dec;++d)
		{
			nuc_d = seq_d[d-1];
			dp = model.p_substitution[nuc_a][nuc_d];
				
			ex_s &WS = table_s(a,d);
			ex_s &WSP = table_s(a-1,d-1);
			noalias(WS) = WSP*dp;
			if(nuc_a == nN || nuc_d == nN)
			{
				WS[xE] += WS[xW]*model.p_tv;
				WS[xEF] += WSP[xF]*model.p_tv*dp;
			}
			else if(nuc_a != nuc_d && g_pupy[nuc_a] != g_pupy[nuc_d])
			{
				WS[xE] += WS[xW];
				WS[xEF] += WSP[xF]*dp;
			}
			
			dp = model.p_open;
			{
				ex_s &WSK = table_s(a,d-1);
				row_cache *= model.p_extend;
				noalias(row_cache) += WSK;
				row_cache[xF] += WSK[xW];
				row_cache[xEF] += WSK[xE];
				noalias(WS) += row_cache*dp;
			}
			
			{
				ex_s &WSK = table_s(a-1,d);
				col_cache[d] *= model.p_extend;
				noalias(col_cache[d]) += WSK;
				col_cache[d][xF] += WSK[xW];
				col_cache[d][xEF] += WSK[xE];
				noalias(WS) += col_cache[d]*dp;
			}
		}
	}
	ex1 = table_s(sz_anc,sz_dec);
	ex1 /= w;
	ex[eVG] = ex1[xEF];
	output() << "." << std::flush;

	// MY
	row_cache.clear();
	fill(col_cache.begin(), col_cache.end(), szero);
	table_s(0,0)[xW] = prob_scale;
	dp = model.p_open;
	for(size_t d = 1; d <= sz_dec; ++d)
	{
		ex_s &WS = table_s(0,d);
		ex_s &WSK = table_s(0,d-1);
		WS.clear();
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		row_cache[xF] += row_cache[xW];
		row_cache[xEF] += row_cache[xE];
		noalias(WS) += row_cache*dp;
	}
	row_cache.clear();
	for(size_t a = 1; a <= sz_anc; ++a)
	{
		ex_s &WS = table_s(a,0);
		ex_s &WSK = table_s(a-1,0);
		WS.clear();
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		row_cache[xF] += row_cache[xW];
		row_cache[xEF] += row_cache[xE];
		noalias(WS) += row_cache*dp;
	}
	for(size_t a=1;a<=sz_anc;++a)
	{
		row_cache.clear();
		nuc_a = seq_a[a-1];
		for(size_t d=1;d<=sz_dec;++d)
		{
			nuc_d = seq_d[d-1];
			dp = model.p_substitution[nuc_a][nuc_d];
				
			ex_s &WS = table_s(a,d);
			ex_s &WSP = table_s(a-1,d-1);
			noalias(WS) = WSP*dp;
			if(nuc_a == nN || nuc_d == nN)
			{
				WS[xE] += WS[xW]*model.p_match;
				WS[xEF] += WSP[xF]*model.p_match*dp;
			}
			else if(nuc_a == nuc_d)
			{
				WS[xE] += WS[xW];
				WS[xEF] += WSP[xF]*dp;
			}
			
			dp = model.p_open;
			{
				ex_s &WSK = table_s(a,d-1);
				row_cache *= model.p_extend;
				noalias(row_cache) += WSK;
				row_cache[xF] += row_cache[xW];
				row_cache[xEF] += row_cache[xE];
				noalias(WS) += row_cache*dp;
			}
			{
				ex_s &WSK = table_s(a-1,d);
				col_cache[d] *= model.p_extend;
				noalias(col_cache[d]) += WSK;
				col_cache[d][xF] += col_cache[d][xW];
				col_cache[d][xEF] += col_cache[d][xE];
				noalias(WS) += col_cache[d]*dp;
			}
		}
	}
	ex1 = table_s(sz_anc,sz_dec);
	ex1 /= w;
	ex[eMY] = ex1[xEF];
	output() << "." << std::flush;

	// SY
	row_cache.clear();
	fill(col_cache.begin(), col_cache.end(), szero);
	table_s(0,0)[xW] = prob_scale;
	dp = model.p_open;
	for(size_t d = 1; d <= sz_dec; ++d)
	{
		ex_s &WS = table_s(0,d);
		ex_s &WSK = table_s(0,d-1);
		WS.clear();
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		row_cache[xF] += row_cache[xW];
		row_cache[xEF] += row_cache[xE];
		noalias(WS) += row_cache*dp;
	}
	row_cache.clear();
	for(size_t a = 1; a <= sz_anc; ++a)
	{
		ex_s &WS = table_s(a,0);
		ex_s &WSK = table_s(a-1,0);
		WS.clear();
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		row_cache[xF] += row_cache[xW];
		row_cache[xEF] += row_cache[xE];
		noalias(WS) += row_cache*dp;
	}
	for(size_t a=1;a<=sz_anc;++a)
	{
		row_cache.clear();
		nuc_a = seq_a[a-1];
		for(size_t d=1;d<=sz_dec;++d)
		{
			nuc_d = seq_d[d-1];
			dp = model.p_substitution[nuc_a][nuc_d];
				
			ex_s &WS = table_s(a,d);
			ex_s &WSP = table_s(a-1,d-1);
			noalias(WS) = WSP*dp;
			if(nuc_a == nN || nuc_d == nN)
			{
				WS[xE] += WS[xW]*model.p_ts;
				WS[xEF] += WSP[xF]*model.p_ts*dp;
			}
			else if(nuc_a != nuc_d && g_pupy[nuc_a] == g_pupy[nuc_d])
			{
				WS[xE] += WS[xW];
				WS[xEF] += WSP[xF]*dp;
			}
			
			dp = model.p_open;
			{
				ex_s &WSK = table_s(a,d-1);
				row_cache *= model.p_extend;
				noalias(row_cache) += WSK;
				row_cache[xF] += row_cache[xW];
				row_cache[xEF] += row_cache[xE];
				noalias(WS) += row_cache*dp;
			}
			{
				ex_s &WSK = table_s(a-1,d);
				col_cache[d] *= model.p_extend;
				noalias(col_cache[d]) += WSK;
				col_cache[d][xF] += col_cache[d][xW];
				col_cache[d][xEF] += col_cache[d][xE];
				noalias(WS) += col_cache[d]*dp;
			}
		}
	}
	ex1 = table_s(sz_anc,sz_dec);
	ex1 /= w;
	ex[eSY] = ex1[xEF];
	output() << "." << std::flush;
	
	// VY
	row_cache.clear();
	fill(col_cache.begin(), col_cache.end(), szero);
	table_s(0,0)[xW] = prob_scale;
	dp = model.p_open;
	for(size_t d = 1; d <= sz_dec; ++d)
	{
		ex_s &WS = table_s(0,d);
		ex_s &WSK = table_s(0,d-1);
		WS.clear();
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		row_cache[xF] += row_cache[xW];
		row_cache[xEF] += row_cache[xE];
		noalias(WS) += row_cache*dp;
	}
	row_cache.clear();
	for(size_t a = 1; a <= sz_anc; ++a)
	{
		ex_s &WS = table_s(a,0);
		ex_s &WSK = table_s(a-1,0);
		WS.clear();
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		row_cache[xF] += row_cache[xW];
		row_cache[xEF] += row_cache[xE];
		noalias(WS) += row_cache*dp;
	}
	for(size_t a=1;a<=sz_anc;++a)
	{
		row_cache.clear();
		nuc_a = seq_a[a-1];
		for(size_t d=1;d<=sz_dec;++d)
		{
			nuc_d = seq_d[d-1];
			dp = model.p_substitution[nuc_a][nuc_d];
				
			ex_s &WS = table_s(a,d);
			ex_s &WSP = table_s(a-1,d-1);
			noalias(WS) = WSP*dp;
			if(nuc_a == nN || nuc_d == nN)
			{
				WS[xE] += WS[xW]*model.p_tv;
				WS[xEF] += WSP[xF]*model.p_tv*dp;
			}
			else if(nuc_a != nuc_d && g_pupy[nuc_a] != g_pupy[nuc_d])
			{
				WS[xE] += WS[xW];
				WS[xEF] += WSP[xF]*dp;
			}
			
			dp = model.p_open;
			{
				ex_s &WSK = table_s(a,d-1);
				row_cache *= model.p_extend;
				noalias(row_cache) += WSK;
				row_cache[xF] += row_cache[xW];
				row_cache[xEF] += row_cache[xE];
				noalias(WS) += row_cache*dp;
			}
			
			{
				ex_s &WSK = table_s(a-1,d);
				col_cache[d] *= model.p_extend;
				noalias(col_cache[d]) += WSK;
				col_cache[d][xF] += col_cache[d][xW];
				col_cache[d][xEF] += col_cache[d][xE];
				noalias(WS) += col_cache[d]*dp;
			}
		}
	}
	ex1 = table_s(sz_anc,sz_dec);
	ex1 /= w;
	ex[eVY] = ex1[xEF];
	output() << "." << std::flush;
	
	// GY
	row_cache.clear();
	fill(col_cache.begin(), col_cache.end(), szero);
	table_s(0,0)[xW] = prob_scale;
	dp = model.p_open;
	for(size_t d = 1; d <= sz_dec; ++d)
	{
		ex_s &WS = table_s(0,d);
		ex_s &WSK = table_s(0,d-1);
		WS.clear();
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		row_cache[xEF] += WSK[xF]+row_cache[xE]+row_cache[xW];
		row_cache[xE] += WSK[xW];
		row_cache[xF] += row_cache[xW];
		noalias(WS) += row_cache*dp;
	}
	row_cache.clear();
	for(size_t a = 1; a <= sz_anc; ++a)
	{
		ex_s &WS = table_s(a,0);
		ex_s &WSK = table_s(a-1,0);
		WS.clear();
		row_cache *= model.p_extend;
		noalias(row_cache) += WSK;
		row_cache[xEF] += WSK[xF]+row_cache[xE]+row_cache[xW];
		row_cache[xE] += WSK[xW];
		row_cache[xF] += row_cache[xW];
		noalias(WS) += row_cache*dp;
	}
	for(size_t a=1;a<=sz_anc;++a)
	{
		row_cache.clear();
		nuc_a = seq_a[a-1];
		for(size_t d=1;d<=sz_dec;++d)
		{
			nuc_d = seq_d[d-1];
			dp = model.p_substitution[nuc_a][nuc_d];
				
			ex_s &WS = table_s(a,d);
			ex_s &WSP = table_s(a-1,d-1);
			noalias(WS) = WSP*dp;
			
			dp = model.p_open;
			{
				ex_s &WSK = table_s(a,d-1);
				row_cache *= model.p_extend;
				noalias(row_cache) += WSK;
				row_cache[xEF] += WSK[xF]+row_cache[xE]+row_cache[xW];
				row_cache[xE] += WSK[xW];
				row_cache[xF] += row_cache[xW];
				noalias(WS) += row_cache*dp;
			}
			
			{
				ex_s &WSK = table_s(a-1,d);
				col_cache[d] *= model.p_extend;
				noalias(col_cache[d]) += WSK;
				col_cache[d][xEF] += WSK[xF]+col_cache[d][xE]+col_cache[d][xW];
				col_cache[d][xE] += WSK[xW];
				col_cache[d][xF] += col_cache[d][xW];
				noalias(WS) += col_cache[d]*dp;
			}
		}
	}
	ex1 = table_s(sz_anc,sz_dec);
	ex1 /= w;
	ex[eGY] = ex1[xEF];
	output() << "." << std::flush;

	return log(w) + log(model.p_end) - log(prob_scale)
		+ static_cast<double>(sz_anc+sz_dec)*log(model.nuc_scale)
		+ static_cast<double>(get_seq_info()[index][nN])*log(model.amb_scale);
}

void covar_k2p_geo::covariance(covar_type &v, const params_type &p, const exvec_type &ev)
{
	ex_type ex;
	ex.clear();
	
	for(exvec_type::const_iterator cit=ev.begin();cit != ev.end();++cit)
		add_expectation(ex, *cit);
	
	double W = get_seqs().size();
	double M = ex[eM], S = ex[eS], V = ex[eV], G = ex[eG], Y = ex[eY];
	double MM = ex[eMM], SS = ex[eSS], VV = ex[eVV], GG = ex[eGG], YY = ex[eYY];
	double MS = ex[eMS], MV = ex[eMV], MG = ex[eMG], MY = ex[eMY], SV = ex[eSV];
	double SG = ex[eSG], SY = ex[eSY], VG = ex[eVG], VY = ex[eVY], GY = ex[eGY];
		
	double t = p[pT], k = p[pK], r = p[pR], q = p[pQ], a = p[pA];
	
	double ek = 1.0+k;
	double ek2 = 1.0 + 2.0*k;
	double ea = exp(-t*ek2/ek);
	double eb = exp(-2.0*t/ek);
	double er = exp(-2.0*t*r);
	
	double da = 4.0/(1.0 - 2.0*ea + eb);
	double db = 4.0/(1.0 + 2.0*ea + eb);
	double dc = eb-1.0;
	double de = -2.0*eb/(dc*ek);
	double df = er/(1.0-er);
	
	double tt = ((8.0*de*(V + 2.0*ek*r*(MV + SV - df*VG)))/ek + 4.0*V*sq(de) +
		(2.0*da*(-2.0*eb*(S + 2.0*ek*r*(MS - df*SG + SS)) + 2.0*eb*ek*(de - 2.0*r)*SV +
		ea*ek2*(4.0*ek*MS*r + S + 2.0*k*S -
		2.0*ek*(2.0*df*r*SG + de*SV - 2.0*r*(SS + SV)))))/sq(ek) -
		(2.0*db*(ea*ek2*(M + 2.0*k*M +
		2.0*ek*(-(de*MV) + 2.0*(-(df*MG) + MM + MS + MV)*r)) +
		eb*(2.0*M + da*eb*MS +
		2.0*ek*(-(de*MV) + 2.0*(-(df*MG) + MM + MS + MV)*r)) -
		da*MS*sq(ea + 2.0*ea*k)))/sq(ek) +
		((S - SS)*sq(da)*sq(ea - eb + 2.0*ea*k))/sq(ek) +
		((M - MM)*sq(db)*sq(ea + eb + 2.0*ea*k))/sq(ek) +
		16.0*(-MM - 2.0*MS - 2.0*MV - SS - 2.0*SV +
		df*(G + df*G - df*GG + 2.0*(MG + SG + VG)))*sq(r) -
		(4.0*VV*sq(de)*sq(eb + dc*ek*r))/sq(eb))/4.0
		;
	
	double tk = (-2.0*da*dc*(dc*(ea + eb)*ek*S +
		(-(dc*ea*(2.0*ek*MS*r + S + 2.0*k*S + 2.0*ek*r*(-(df*SG) + SS))) -
		2.0*dc*eb*(S + ek*r*(MS - df*SG + SS)) +
		eb*(-2.0*eb + dc*de*ek - 2.0*dc*ek*r)*SV +
		ea*(dc*de*ek + 2.0*eb*ek2 - 2.0*dc*ek*r)*SV)*t) +
		8.0*eb*(2.0*dc*ek*MV*r*t - (dc*ek + 2.0*t)*V + 2.0*eb*t*VV +
		2.0*dc*ek*r*t*(SV - df*VG + VV)) +
		(ea + eb)*(ea - eb + 2.0*ea*k)*(S - SS)*t*sq(da)*sq(dc) +
		(ea - eb)*(ea + eb + 2.0*ea*k)*(M - MM)*t*sq(db)*sq(dc) +
		2.0*db*dc*(eb*(-(dc*M*(ek - 2.0*t)) +
		(da*dc*eb*MS + 2.0*eb*MV - dc*de*ek*MV +
		2.0*dc*ek*(-(df*MG) + MM + MS + MV)*r)*t) +
		ea*((dc*de*ek*MV + 2.0*eb*ek2*MV -
		2.0*dc*ek*(-(df*MG) + MM + MS + MV)*r)*t + dc*M*(ek - t - 2.0*k*t))
		+ da*dc*ek2*MS*t*sq(ea)))/(4.0*ek*sq(ek)*sq(dc))
		;

	
	double tr = (-2.0*k*M - 2.0*k*S + db*ea*MM*t + 2.0*db*ea*k*MM*t - da*ea*MS*t + db*ea*MS*t -
		2.0*da*ea*k*MS*t + 2.0*db*ea*k*MS*t - 2.0*de*MV*t + db*ea*MV*t -
		2.0*de*k*MV*t + 2.0*db*ea*k*MV*t + 4.0*MM*r*t + 4.0*k*MM*r*t + 8.0*MS*r*t +
		8.0*k*MS*r*t + 8.0*MV*r*t + 8.0*k*MV*r*t - da*ea*SS*t - 2.0*da*ea*k*SS*t +
		4.0*r*SS*t + 4.0*k*r*SS*t - 2.0*de*SV*t - da*ea*SV*t - 2.0*de*k*SV*t -
		2.0*da*ea*k*SV*t + 8.0*r*SV*t + 8.0*k*r*SV*t - 2.0*k*V - 2.0*(M + S + V) +
		dc*df*(2.0*ek*G*(-1.0 + 2.0*r*t) +
		t*(db*(ea + eb + 2.0*ea*k)*MG + 8.0*ek*MG*r -
		(da*(ea - eb + 2.0*ea*k) - 8.0*ek*r)*SG - 2.0*ek*(de - 4.0*r)*VG)) +
		4.0*ek*r*t*VV + eb*(2.0*ek*M + 2.0*S - db*(-1.0 + ea)*MM*t + 2.0*V +
		2.0*k*(S - db*ea*(MM + MS + MV)*t + V +
		t*(de*(MV + SV) + da*ea*(MS + SS + SV) -
		2.0*r*(MM + 2.0*MS + 2.0*MV + SS + 2.0*SV + VV))) +
		t*(-(db*(-1.0 + ea)*(MS + MV)) + da*(1.0 + ea)*(MS + SS + SV) +
		2.0*(de*(MV + SV) - 2.0*VV - 2.0*r*(MM + 2.0*MS + 2.0*MV + SS + 2.0*SV + VV)))
		) + 4.0*dc*ek*(G - GG)*r*t*sq(df) -
		(db*(MM + MS + MV) + da*(MS + SS + SV))*t*sq(eb))/(dc*ek)
		;

	double tq = (db*dc*(ea + eb + 2.0*ea*k)*(MY - MG*q) - 4.0*MY*r + 4.0*eb*MY*r - 4.0*k*MY*r +
		4.0*eb*k*MY*r + 4.0*MG*q*r - 4.0*eb*MG*q*r + 4.0*k*MG*q*r - 4.0*eb*k*MG*q*r -
		4.0*dc*df*ek*(GY - GG*q)*r - da*ea*q*SG + da*eb*q*SG + da*ea*eb*q*SG -
		2.0*da*ea*k*q*SG + 2.0*da*ea*eb*k*q*SG + 4.0*q*r*SG - 4.0*eb*q*r*SG +
		4.0*k*q*r*SG - 4.0*eb*k*q*r*SG + da*ea*SY - da*eb*SY - da*ea*eb*SY +
		2.0*da*ea*k*SY - 2.0*da*ea*eb*k*SY - 4.0*r*SY + 4.0*eb*r*SY - 4.0*k*r*SY +
		4.0*eb*k*r*SY - 2.0*de*q*VG + 2.0*de*eb*q*VG - 2.0*de*k*q*VG + 2.0*de*eb*k*q*VG +
		4.0*q*r*VG - 4.0*eb*q*r*VG + 4.0*k*q*r*VG - 4.0*eb*k*q*r*VG +
		4.0*(eb + dc*ek*r)*VY - da*q*SG*sq(eb) + da*SY*sq(eb))/
		(2.0*dc*ek*(-1.0 + q)*q)
		;
	
	double ta = (db*(ea + eb + 2.0*ea*k)*(MG + MM + MS + MV - a*M*W) -
		da*(ea - eb + 2.0*ea*k)*(MS + SG + SS + SV - a*S*W) +
		2.0*ek*(-(de*(MV + SV + VG + VV - a*V*W)) +
		2.0*r*(MG + MM + 2.0*MS + 2.0*MV + SG + SS + 2.0*SV + VG -
		df*(GG + MG + SG + VG) + VV - a*(-(df*G) + M + S + V)*W)))/
		(2.0*a*(1.0 + a)*ek)
		;

	double kk = (t*(2.0*da*dc*(2.0*dc*(ea + eb)*ek*S +
		(dc*(ea - 2.0*eb)*S - 4.0*eb*(ea + eb)*SV)*t) +
		16.0*eb*((dc*ek + t)*V - eb*t*VV) +
		2.0*db*dc*(-2.0*dc*(ea - eb)*ek*M +
		t*(4.0*(ea - eb)*eb*MV -
		dc*(ea*M + eb*(2.0*M + da*eb*MS) - da*MS*sq(ea)))) +
		(M - MM)*t*sq(db)*sq(dc)*sq(ea - eb) +
		(S - SS)*t*sq(da)*sq(dc)*sq(ea + eb)))/(4.0*sq(sq(ek))*sq(dc))
		;
	
	double kr = ((db*dc*(ea - eb)*(df*MG - MM - MS - MV) +
		da*dc*(ea + eb)*(MS - df*SG + SS + SV) + 4.0*eb*(MV + SV - df*VG + VV))*
		sq(t))/(dc*sq(ek))
		;
	
   double kq = (t*(db*dc*(ea - eb)*(MY - MG*q) + da*dc*(ea + eb)*(q*SG - SY) +
		4.0*eb*(q*VG - VY)))/(2.0*dc*(-1.0 + q)*q*sq(ek))
		;
	
	double ka = (t*(db*dc*(ea - eb)*(MG + MM + MS + MV - a*M*W) -
		da*dc*(ea + eb)*(MS + SG + SS + SV - a*S*W) -
		4.0*eb*(MV + SV + VG + VV - a*V*W)))/(2.0*a*(1.0 + a)*dc*sq(ek))
		;

   double rr = 4.0*(-MM - 2.0*MS - 2.0*MV - SS - 2.0*SV +
		df*(G + df*G - df*GG + 2.0*(MG + SG + VG)) - VV)*sq(t)
		;

   double rq = (2.0*t*(MY + df*(-GY + GG*q) + SY - q*(MG + SG + VG) + VY))/((-1.0 + q)*q)
		;
	
	double ra = (2.0*t*(-(df*(-1.0 + er)*(GG + MG + SG + VG)) +
		(-1.0 + er)*(MG + MM + 2.0*MS + 2.0*MV + SG + SS + 2.0*SV + VG + VV) -
		a*(-M - S - V + er*(G + M + S + V))*W))/(a*(1.0 + a)*(-1.0 + er))
		;


	double qq = -((q*(-2.0*GY + (G + GG)*q - 2.0*Y) + Y + YY)/(sq(-1.0 + q)*sq(q)));
	
	double qa = -((GY + MY + SY + VY - q*(GG + MG + SG + VG - a*G*W) - a*W*Y)/
		(a*(1.0 + a)*(-1.0 + q)*q))
		;
	
	double aa = -((GG - M + 2.0*MG + MM + 2.0*MS + 2.0*MV - S + 2.0*SG + SS + 2.0*SV - V + 2.0*VG +
		VV + a*(1.0 + W)*(-2.0*(M + S + V) + a*W) - G*(1.0 + 2.0*a*(1.0 + W)))/
		(sq(a)*sq(1.0 + a)))
		;
	
	covar_type fim(p.size(), p.size());
	fim(pT,pT) = tt; fim(pT,pK) = tk; fim(pT,pR) = tr; fim(pT,pQ) = tq; fim(pT,pA) = ta;
	fim(pK,pT) = tk; fim(pK,pK) = kk; fim(pK,pR) = kr; fim(pK,pQ) = kq; fim(pK,pA) = ka;
	fim(pR,pT) = tr; fim(pR,pK) = kr; fim(pR,pR) = rr; fim(pR,pQ) = rq; fim(pR,pA) = ra;
	fim(pQ,pT) = tq; fim(pQ,pK) = kq; fim(pQ,pR) = rq; fim(pQ,pQ) = qq; fim(pQ,pA) = qa;
	fim(pA,pT) = ta; fim(pA,pK) = ka; fim(pA,pR) = ra; fim(pA,pQ) = qa; fim(pA,pA) = aa;
		
	invert_matrix(fim, v);
}


void covar_k2p_geo::add_expectation(ex_type &ex, const ex_type &ad) const
{
	cout << ad << endl;
	
	ex[eMM] += ad[eMM] + 2.0*ex[eM]*ad[eM];
	ex[eSS] += ad[eSS] + 2.0*ex[eS]*ad[eS];
	ex[eVV] += ad[eVV] + 2.0*ex[eV]*ad[eV];
	ex[eGG] += ad[eGG] + 2.0*ex[eG]*ad[eG];
	ex[eYY] += ad[eYY] + 2.0*ex[eY]*ad[eY];
	
	ex[eMS] += ad[eMS] + ex[eM]*ad[eS] + ex[eS]*ad[eM];
	ex[eMV] += ad[eMV] + ex[eM]*ad[eV] + ex[eV]*ad[eM];
	ex[eMG] += ad[eMG] + ex[eM]*ad[eG] + ex[eG]*ad[eM];
	ex[eMY] += ad[eMY] + ex[eM]*ad[eY] + ex[eY]*ad[eM];
	ex[eSV] += ad[eSV] + ex[eS]*ad[eV] + ex[eV]*ad[eS];
	ex[eSG] += ad[eSG] + ex[eS]*ad[eG] + ex[eG]*ad[eS];
	ex[eSY] += ad[eSY] + ex[eS]*ad[eY] + ex[eY]*ad[eS];
	ex[eVG] += ad[eVG] + ex[eV]*ad[eG] + ex[eG]*ad[eV];
	ex[eVY] += ad[eVY] + ex[eV]*ad[eY] + ex[eY]*ad[eV];
	ex[eGY] += ad[eGY] + ex[eG]*ad[eY] + ex[eY]*ad[eG];

	ex[eM] += ad[eM];
	ex[eS] += ad[eS];
	ex[eV] += ad[eV];
	ex[eG] += ad[eG];
	ex[eY] += ad[eY];
}


std::string covar_k2p_geo::state() const
{
	ex_type ex;
	ex.clear();
	const exvec_type &ev = get_exvec();
	for(exvec_type::const_iterator cit=ev.begin();cit != ev.end();++cit)
		add_expectation(ex, *cit);
	
	const covar_type &v = get_covar();
	
	ostringstream ostr;
	ostr << std::setprecision(DBL_DIG);
	ostr << "ex = " << ex << endl;
	ostr << "var = " << v << endl;
	return ostr.str();
}

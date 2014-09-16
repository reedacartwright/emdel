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

#include <boost/math/special_functions/zeta.hpp>

#include <algorithm>
#include <cfloat>
#include <iostream>

#include "hist_k2p.h"

#include "series.h"

inline double pow(double p, size_t x){return std::pow(p, (int)x);}
inline double log(size_t x) {return log((double)x);}
inline double zeta(double z) { return boost::math::zeta<double>(z); }
inline double sq(double s) { return s*s; }
inline double csch(double x) { return 1.0/sinh(x); }

using namespace std;

extern const int g_pupy[5];
extern const int g_atgc[5];

/***************************************************************************
 * class hist_k2p_zeta                                                     *
 ***************************************************************************/

void hist_k2p_zeta::preallocate(size_t maxa, size_t maxd)
{
	sz_height = maxa+1;
	sz_width = maxd+1;
	
	size_t sz_max = std::max(sz_height, sz_width);
	p_indel_size.resize(sz_max, 0.0);
	
	ll_nucs.resize(get_seqs().size());
	
	vector<double>::iterator llit = ll_nucs.begin();
	
	for(seq_pair_vector::const_iterator cit = get_seqs().begin();
		cit != get_seqs().end(); ++cit)
	{
		double ncount[nSize] = {0.0,0.0,0.0,0.0,0.0};	
		for(sequence::const_iterator nit = cit->first.begin();
		    nit != cit->first.end(); ++nit)
			++ncount[(int)*nit];
		for(sequence::const_iterator nit = cit->second.begin();
		    nit != cit->second.end(); ++nit)
			++ncount[(int)*nit];
		*llit = log(0.25)*(ncount[nA]+ncount[nC]+ncount[nG]+ncount[nT])
			+ log(1.0)*(ncount[nN]);
		++llit;
	}
}

void hist_k2p_zeta::expectation_setup(const params_type &params)
{
	p_ts = 0.25-0.5*exp(-params[pT]*(2.0*params[pK]+1.0)/(params[pK]+1.0))
			+ 0.25*exp(-2.0*params[pT]/(params[pK]+1.0));
	p_tv = 0.5-0.5*exp(-2.0*params[pT]/(params[pK]+1.0));
	p_match = 1.0-p_ts-p_tv;
	
	p_end = 1.0/(params[pA]+1.0);
	p_2h =  (1.0-p_end)*exp(-2.0*params[pR]*params[pT]);
	p_2g =  (1.0-p_end)*0.5*(1.0-exp(-2.0*params[pR]*params[pT]));
	
	for(size_t i=0;i<nN;++i)
	{
		p_substitution[i][i] = 4.0*p_match*p_2h;
		for(size_t j=i+1;j<nN;++j)
		{
			p_substitution[j][i] = p_substitution[i][j] = p_2h*
				((g_pupy[i] == g_pupy[j]) ? 4.0*p_ts : 2.0*p_tv);
		}
		p_substitution[nN][i] = p_substitution[i][nN] = p_2h;
	}
	p_substitution[nN][nN] = p_2h;
		
	double z = params[pZ];
	double Z = p_2g/zeta(z);
	for(size_t u = 1; u < p_indel_size.size(); ++u)
		p_indel_size[u] = pow((double)u, -z)*Z;
}

void hist_k2p_zeta::expectation(ex_type &ex, ex_type &ex2, const sequence &seq_a, const sequence &seq_d, size_t index) const
{
	for(size_t u = 0; u < ex.size(); ++u)
	{
		double d, d2;
		expectation(d, d2, u, seq_a, seq_d, index);
		ex[u] = d;
		ex2[u] = d2;
	}
}

void hist_k2p_zeta::expectation(double &ex, double &ex2, size_t gap, const sequence &seq_a, const sequence &seq_d, size_t index) const
{
	size_t sz_anc = seq_a.size();
	size_t sz_dec = seq_d.size();
	
	xx_table table_w(sz_anc+1,sz_dec+1,ublas::zero_vector<double>(xSize));
	
	table_w(0,0)[xW] = prob_scale;
	double dt,dp;
	
	for(size_t d = 1; d <= sz_dec; ++d)
	{
		xx_type &W = table_w(0,d);
		for(size_t k = d; k > 0; --k)
		{
			dp = p_indel_size[k];
			xx_type &WK = table_w(0,d-k);
			noalias(W) += WK*dp;
			if(k == gap || gap == 0)
			{
				W[xX] += WK[xW]*dp;
				W[xXX] += 2.0*WK[xX]*dp + WK[xW]*dp;
			}
		}
	}
		
	for(size_t a = 1; a <= sz_anc; ++a)
	{
		xx_type &W = table_w(a,0);
		for(size_t k = a; k > 0; --k)
		{
			dp = p_indel_size[k];
			xx_type &WK = table_w(a-k,0);
			noalias(W) += WK*dp;
			if(k == gap || gap == 0)
			{
				W[xX] += WK[xW]*dp;
				W[xXX] += 2.0*WK[xX]*dp + WK[xW]*dp;
			}
		}
	}

	int nuc_a, nuc_d;
		
	for(size_t a=1;a<=sz_anc;++a)
	{
		for(size_t d=1;d<= sz_dec;++d)
		{
			nuc_a = seq_a[a-1];
			nuc_d = seq_d[d-1];
				
			xx_type &W = table_w(a,d);
			noalias(W) = table_w(a-1,d-1)*p_substitution[nuc_a][nuc_d];
						
			for(size_t k = d; k > 0; --k)
			{
				dp = p_indel_size[k];
				xx_type &WK = table_w(a,d-k);
				noalias(W) += WK*dp;
				if(k == gap || gap == 0)
				{
					W[xX] += WK[xW]*dp;
					W[xXX] += 2.0*WK[xX]*dp + WK[xW]*dp;
				}
			}
			for(size_t k = a; k > 0; --k)
			{		
				dp = p_indel_size[k];
				xx_type &WK = table_w(a-k,d);
				noalias(W) += WK*dp;
				if(k == gap || gap == 0)
				{
					W[xX] += WK[xW]*dp;
					W[xXX] += 2.0*WK[xX]*dp + WK[xW]*dp;
				}
			}
		}
	}
	
	xx_type xx = table_w(sz_anc,sz_dec);
	double w = xx[xW];
	xx /= w;
	ex = xx[xX];
	ex2 = xx[xXX];
}

std::string hist_k2p_zeta::state() const
{
	ex_type mean, var;
	mean.clear();
	var.clear();
	
	for(size_t u = 0; u<exes.size();++u)
	{

		var += ex2s[u] + 2.0*element_prod(mean, exes[u]);
		mean += exes[u];
	}
	noalias(var) -= element_prod(mean, mean);
	
	ostringstream ostr;
	ostr << std::setprecision(DBL_DIG);
	ostr << "Size\tMean\tVariance" << endl;
	for(size_t u = 0; u<mean.size();++u)
	{
		ostr << u << "\t" << mean[u] << "\t" << var[u] << endl;
	}
	return ostr.str();
}

/***************************************************************************
 * class hist_k2p_geo                                                      *
 ***************************************************************************/

void hist_k2p_geo::preallocate(size_t maxa, size_t maxd)
{
	sz_height = maxa+1;
	sz_width = maxd+1;
	
	size_t sz_max = std::max(sz_height, sz_width);
	p_indel_size.resize(sz_max, 0.0);
	cache_size.resize(sz_max, 0.0);
	for(size_t u=1;u<sz_max;++u)
		cache_size[u] = double(u);
	
	ll_nucs.resize(get_seqs().size());
	
	vector<double>::iterator llit = ll_nucs.begin();
	
	for(seq_pair_vector::const_iterator cit = get_seqs().begin();
		cit != get_seqs().end(); ++cit)
	{
		double ncount[nSize] = {0.0,0.0,0.0,0.0,0.0};	
		for(sequence::const_iterator nit = cit->first.begin();
		    nit != cit->first.end(); ++nit)
			++ncount[(int)*nit];
		for(sequence::const_iterator nit = cit->second.begin();
		    nit != cit->second.end(); ++nit)
			++ncount[(int)*nit];
		*llit = log(0.25)*(ncount[nA]+ncount[nC]+ncount[nG]+ncount[nT])
			+ log(1.0)*(ncount[nN]);
		++llit;
	}	
}

void hist_k2p_geo::expectation_setup(const params_type &params)
{
	p_ts = 0.25-0.5*exp(-params[pT]*(2.0*params[pK]+1.0)/(params[pK]+1.0))
			+ 0.25*exp(-2.0*params[pT]/(params[pK]+1.0));
	p_tv = 0.5-0.5*exp(-2.0*params[pT]/(params[pK]+1.0));
	p_match = 1.0-p_ts-p_tv;
	
	p_end = 1.0/(params[pA]+1.0);
	p_2h =  (1.0-p_end)*exp(-2.0*params[pR]*params[pT]);
	p_2g =  (1.0-p_end)*0.5*(1.0-exp(-2.0*params[pR]*params[pT]));
	
	
	for(size_t i=0;i<nN;++i)
	{
		p_substitution[i][i] = 4.0*p_match*p_2h;
		for(size_t j=i+1;j<nN;++j)
		{
			p_substitution[j][i] = p_substitution[i][j] = 4.0*p_2h*
				((g_pupy[i] == g_pupy[j]) ? p_ts : 0.5*p_tv);
		}
		p_substitution[nN][i] = p_substitution[i][nN] = 4.0*p_2h;
	}
	p_substitution[nN][nN] = 16.0*p_2h;
		
	double qq = 1.0/params[pQ];
	double dp = 1.0/(1.0-qq);
	for(size_t u = 1; u < p_indel_size.size(); ++u)
	{
		dp *= 1.0-qq;
		p_indel_size[u] = p_2g*qq*dp;
	}
	p_open = qq*p_2g;
	p_extend = (1.0-qq);
}

void hist_k2p_geo::expectation(ex_type &ex, ex_type &ex2, const sequence &seq_a, const sequence &seq_d, size_t index) const
{
	for(size_t u = 0; u < ex.size(); ++u)
	{
		double d, d2;
		expectation(d, d2, u, seq_a, seq_d, index);
		ex[u] = d;
		ex2[u] = d2;
	}
}

void hist_k2p_geo::expectation(double &ex, double &ex2, size_t gap, const sequence &seq_a, const sequence &seq_d, size_t index) const
{
	size_t sz_anc = seq_a.size();
	size_t sz_dec = seq_d.size();
	
	double p_gap = p_open*pow(p_extend, gap);
	
	xx_table table_w(sz_anc+1,sz_dec+1,ublas::zero_vector<double>(xSize));
	
	table_w(0,0)[xW] = prob_scale;
	
	xx_type row_cache;
	vector<xx_type> col_cache(sz_dec+1, ublas::zero_vector<double>(xSize));
	row_cache.clear();
	for(size_t d = 1; d <= sz_dec; ++d)
	{
		xx_type &W = table_w(0,d);
	
		row_cache *= p_extend;
		noalias(row_cache) += table_w(0,d-1)*p_open;
		if(gap == 0)
		{
			xx_type &WK = table_w(0,d-1);
			row_cache[xX] += WK[xW]*p_open;
			row_cache[xXX] += (WK[xX]+WK[xX]+WK[xW])*p_open;
		}
		else if(d >= gap)
		{
			xx_type &WK = table_w(0,d-gap);
			W[xX] += WK[xW]*p_gap;
			W[xXX] += (WK[xX]+WK[xX]+WK[xW])*p_gap;
		}
		noalias(W) += row_cache;
	}
	
	row_cache.clear();
	for(size_t a = 1; a <= sz_anc; ++a)
	{
		xx_type &W = table_w(a,0);
	
		row_cache *= p_extend;
		noalias(row_cache) += table_w(a-1,0)*p_open;
		if(gap == 0)
		{
			xx_type &WK = table_w(a-1,0);
			row_cache[xX] += WK[xW]*p_open;
			row_cache[xXX] += (WK[xX]+WK[xX]+WK[xW])*p_open;
		}
		else if(a >= gap)
		{
			xx_type &WK = table_w(a-gap,0);
			W[xX] += WK[xW]*p_gap;
			W[xXX] += (WK[xX]+WK[xX]+WK[xW])*p_gap;
		}
		noalias(W) += row_cache;
	}
	
	int nuc_a, nuc_d;
		
	for(size_t a=1;a<=sz_anc;++a)
	{
		row_cache.clear();
		nuc_a = seq_a[a-1];
		for(size_t d=1;d<= sz_dec;++d)
		{
			nuc_d = seq_d[d-1];
				
			xx_type &W = table_w(a,d);
			noalias(W) = table_w(a-1,d-1)*p_substitution[nuc_a][nuc_d];
			
			row_cache *= p_extend;
			noalias(row_cache) += table_w(a,d-1)*p_open;
		
			if(gap == 0)
			{
				xx_type &WK = table_w(a,d-1);
				row_cache[xX] += WK[xW]*p_open;
				row_cache[xXX] += (WK[xX]+WK[xX]+WK[xW])*p_open;
			}
			else if(d >= gap)
			{
				xx_type &WK = table_w(a,d-gap);
				W[xX] += WK[xW]*p_gap;
				W[xXX] += (WK[xX]+WK[xX]+WK[xW])*p_gap;
			}
			noalias(W) += row_cache;
			
			col_cache[d] *= p_extend;
			noalias(col_cache[d]) += table_w(a-1,d)*p_open;
		
			if(gap == 0)
			{
				xx_type &WK = table_w(a-1,d);
				col_cache[d][xX] += WK[xW]*p_open;
				col_cache[d][xXX] += (WK[xX]+WK[xX]+WK[xW])*p_open;
			}
			else if(a >= gap)
			{
				xx_type &WK = table_w(a-gap,d);
				W[xX] += WK[xW]*p_gap;
				W[xXX] += (WK[xX]+WK[xX]+WK[xW])*p_gap;
			}
			noalias(W) += col_cache[d];

		}
	}
	xx_type xx = table_w(sz_anc,sz_dec);
	double w = xx[xW];
	xx /= w;
	ex = xx[xX];
	ex2 = xx[xXX];
	cout << gap << " " << ex << " " << ex2 << endl;
}

std::string hist_k2p_geo::state() const
{
	ex_type mean, var;
	mean.clear();
	var.clear();
	
	for(size_t u = 0; u<exes.size();++u)
	{
		var += ex2s[u] + 2.0*element_prod(mean,exes[u]);
		mean += exes[u];
	}
	noalias(var) -= element_prod(mean,mean);
	
	ostringstream ostr;
	ostr << std::setprecision(DBL_DIG);
	ostr << "Size\tExpected\tVariance" << endl;
	for(size_t u = 0; u<mean.size();++u)
	{
		ostr << u << "\t" << mean[u] << "\t" << var[u] << endl;
	}
	return ostr.str();
}


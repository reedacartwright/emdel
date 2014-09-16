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

#include <boost/math/special_functions/zeta.hpp>

#include "models.h"

inline double zeta(double z) { return boost::math::zeta<double>(z); }

extern const int g_pupy[] = {0,1,0,1,-1};
extern const int g_atgc[] = {0,1,1,0,-1};

const model_k2p_zeta::freq_type model_k2p_zeta::k2p_freqs = {0.25,0.25,0.25,0.25};
const model_k2p_geo::freq_type model_k2p_geo::k2p_freqs = {0.25,0.25,0.25,0.25};

/***************************************************************************
 * class model_k2p_zeta                                                    *
 ***************************************************************************/

void model_k2p_zeta::preallocate(size_t maxa, size_t maxd)
{
	size_t sz_max = std::max(maxa, maxd)+1;
	p_indel_size.resize(sz_max, 0.0);
}

void model_k2p_zeta::expectation_setup(const params_type &params, const freq_type &freq)
{
	p_ts = 0.25-0.5*exp(-params[pT]*(2.0*params[pK]+1.0)/(params[pK]+1.0))
			+ 0.25*exp(-2.0*params[pT]/(params[pK]+1.0));
	p_tv = 0.5-0.5*exp(-2.0*params[pT]/(params[pK]+1.0));
	p_match = 1.0-p_ts-p_tv;
	
	p_end = 1.0/(params[pA]+1.0);
	p_2h =  (1.0-p_end)*exp(-2.0*params[pR]*params[pT]);
	p_2g =  (1.0-p_end)*0.5*(1.0-exp(-2.0*params[pR]*params[pT]));
	
	double nuc_scale2 = 0.25*p_2h
		*pow(p_match, p_match)
		*pow(p_ts, p_ts)
		*pow(0.5*p_tv, p_tv);
	nuc_scale = sqrt(nuc_scale2);
	amb_scale = 4.0;
	
	for(size_t i=0;i<nN;++i)
	{
		p_substitution[i][i] = 0.25*p_2h*p_match/nuc_scale2;
		
		for(size_t j=i+1;j<nN;++j)
		{
			p_substitution[j][i] = p_substitution[i][j] = p_2h*
				((g_pupy[i] == g_pupy[j]) ? 0.25*p_ts : 0.125*p_tv)
				/ nuc_scale2;
		}
		p_substitution[nN][i] = p_substitution[i][nN] = 0.0625*p_2h / nuc_scale2;
	}
	p_substitution[nN][nN] = 0.0625*p_2h / nuc_scale2;
		
	double z = params[pZ];
	double Z = p_2g/zeta(z);
	double s = 0.25/nuc_scale;
	double ss = s;
	for(size_t u = 1; u < p_indel_size.size(); ++u)
	{
		p_indel_size[u] = pow((double)u, -z)*Z*ss;
		ss *= s;
	}
}

/***************************************************************************
 * class model_k2p_geo                                                     *
 ***************************************************************************/

void model_k2p_geo::preallocate(size_t maxa, size_t maxd)
{
	size_t sz_max = std::max(maxa, maxd)+1;
	p_indel_size.resize(sz_max, 0.0);
}

void model_k2p_geo::expectation_setup(const params_type &params, const freq_type &freq)
{
	p_ts = 0.25-0.5*exp(-params[pT]*(2.0*params[pK]+1.0)/(params[pK]+1.0))
			+ 0.25*exp(-2.0*params[pT]/(params[pK]+1.0));
	p_tv = 0.5-0.5*exp(-2.0*params[pT]/(params[pK]+1.0));
	p_match = 1.0-p_ts-p_tv;
	
	p_end = 1.0/(params[pA]+1.0);
	p_2h =  (1.0-p_end)*exp(-2.0*params[pR]*params[pT]);
	p_2g =  (1.0-p_end)*0.5*(1.0-exp(-2.0*params[pR]*params[pT]));
	
	double nuc_scale2 = 0.25*p_2h
		*pow(p_match, p_match)
		*pow(p_ts, p_ts)
		*pow(0.5*p_tv, p_tv);
	nuc_scale = sqrt(nuc_scale2);
	amb_scale = 4.0;
	
	for(size_t i=0;i<nN;++i)
	{
		p_substitution[i][i] = 0.25*p_2h*p_match/nuc_scale2;
		
		for(size_t j=i+1;j<nN;++j)
		{
			p_substitution[j][i] = p_substitution[i][j] = p_2h*
				((g_pupy[i] == g_pupy[j]) ? 0.25*p_ts : 0.125*p_tv)
				/ nuc_scale2;
		}
		p_substitution[nN][i] = p_substitution[i][nN] = 0.0625*p_2h / nuc_scale2;
	}
	p_substitution[nN][nN] = 0.0625*p_2h / nuc_scale2;
		
	double s = 0.25/nuc_scale;
	double qq = 1.0/params[pQ];
	double dp = qq*s;
	for(size_t u = 1; u < p_indel_size.size(); ++u)
	{
		p_indel_size[u] = p_2g*dp;
		dp *= (1.0-qq)*s;
	}
	p_open = qq*p_2g*s;
	p_extend = (1.0-qq)*s;
}


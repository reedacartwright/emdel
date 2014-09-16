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

#define _USE_MATH_DEFINES

#include <iostream>

#include <boost/math/special_functions/zeta.hpp>

#include "sample.h"
#include "sample_k2p.h"
#include "ccvector.h"
#include "series.h"
#include "invert_matrix.h"
#include "table.h"
#include "emdel.h"

using namespace std;

void sample_k2p_zeta::preallocate(size_t maxa, size_t maxd)
{
	p.resize(maxa+1,maxd+1,0.0);
}

void sample_k2p_zeta::presample_step(const params_type &params, const sequence &seq_a, const sequence &seq_d)
{
	size_t sz_max = std::max(sz_height, sz_width);
	size_t sz_anc = seq_a.size();
	size_t sz_dec = seq_d.size();
	sa = seq_a;
	sd = seq_d;
	
	const model_type &model = get_model();
	
	p(0,0) = prob_scale;
	
	for(size_t d = 1; d <= sz_dec; ++d) {
		p(0,d) = 0.0;
		for(size_t k = d; k > 0; --k)
			p(0,d) += p(0,d-k)*model.p_indel_size[k];
	}
		
	for(size_t a = 1; a <= sz_anc; ++a) {
		p(a,0) = 0.0;
		for(size_t k = a; k > 0; --k)
			p(a,0) += p(a-k,0)*model.p_indel_size[k];
	}

	for(size_t a=1;a<=sz_anc;++a) {
		for(size_t d=1;d<= sz_dec;++d) {
			double pt = p(a-1,d-1)*model.p_substitution[sa[a-1]][sd[d-1]];

			for(size_t k = d; k > 0; --k)
				pt += p(a,d-k)*model.p_indel_size[k];
			for(size_t k = a; k > 0; --k)
				pt += p(a-k,d)*model.p_indel_size[k];
			p(a,d) = pt;
		}
	}
}

void sample_k2p_zeta::sample_once(std::string &seq_a, std::string &seq_d)
{
	size_t a = sa.size();
	size_t d = sd.size();
	seq_a.clear();
	seq_d.clear();
	
	const model_type &model = get_model();
	
	while(a != 0 || d != 0) {
		double u = p(a,d)*myrand.uniform01();
		double t = (a > 0 && d > 0) ? p(a-1,d-1)*model.p_substitution[sa[a-1]][sd[d-1]] : 0.0;
		if(u < t) {
			// match
			a = a-1;
			d = d-1;
			seq_a.append(1, ccNuc[sa[a]]);
			seq_d.append(1, ccNuc[sd[d]]);
		} else {
			for(size_t k = 1; k <= a || k <= d; ++k) {
				if(k <= a) {
					t += p(a-k,d)*model.p_indel_size[k];
					if(u < t) {
						// gap in d
						seq_d.append(k, '-');
						while(k--) seq_a.append(1, ccNuc[sa[--a]]);
						break;
					}
				}
				if(k <= d)
				{
					t += p(a,d-k)*model.p_indel_size[k];
					if(u < t) {
						// gap in a
						seq_a.append(k, '-');
						while(k--) seq_d.append(1, ccNuc[sd[--d]]);
						break;
					}
				}
			}
		}
		//seq_a.append(1,',');
		//seq_d.append(1,',');
	}
	reverse(seq_a.begin(), seq_a.end());
	reverse(seq_d.begin(), seq_d.end());
}

void sample_k2p_geo::preallocate(size_t maxa, size_t maxd)
{
	p.resize(maxa+1,maxd+1,0.0);
}

void sample_k2p_geo::presample_step(const params_type &params, const sequence &seq_a, const sequence &seq_d)
{
	size_t sz_max = std::max(sz_height, sz_width);
	size_t sz_anc = seq_a.size();
	size_t sz_dec = seq_d.size();
	sa = seq_a;
	sd = seq_d;
	
	const model_type &model = get_model();
	
	double row_cache = 0.0;
	vector<double> col_cache(sz_dec+1, 0.0);
	
	p(0,0) = prob_scale;
	for(size_t d = 1; d <= sz_dec; ++d) {
		row_cache *= model.p_extend;
		row_cache += p(0,d-1)*model.p_open;
		p(0,d) = row_cache;
	}
	
	row_cache = 0.0;
	for(size_t a = 1; a <= sz_anc; ++a) {
		row_cache *= model.p_extend;
		row_cache += p(a-1,0)*model.p_open;
		p(a,0) = row_cache;
	}

	for(size_t a=1;a<=sz_anc;++a) {
		row_cache = 0.0;
		for(size_t d=1;d<= sz_dec;++d) {
			double pt = p(a-1,d-1)*model.p_substitution[sa[a-1]][sd[d-1]];
			row_cache *= model.p_extend;
			row_cache += p(a,d-1)*model.p_open;
			pt += row_cache;
			col_cache[d] *= model.p_extend;
			col_cache[d] += p(a-1,d)*model.p_open;
			pt += col_cache[d];
			p(a,d) = pt;
		}
	}
}

void sample_k2p_geo::sample_once(std::string &seq_a, std::string &seq_d)
{
	size_t a = sa.size();
	size_t d = sd.size();
	seq_a.clear();
	seq_d.clear();
	
	const model_type &model = get_model();
	
	while(a != 0 || d != 0) {
		double u = p(a,d)*myrand.uniform01();
		double t = (a > 0 && d > 0) ? p(a-1,d-1)*model.p_substitution[sa[a-1]][sd[d-1]] : 0.0;
		if(u < t) {
			// match
			a = a-1;
			d = d-1;
			seq_a.append(1, ccNuc[sa[a]]);
			seq_d.append(1, ccNuc[sd[d]]);
		} else {
			for(size_t k = 1; k <= a || k <= d; ++k) {
				if(k <= a) {
					t += p(a-k,d)*model.p_indel_size[k];
					if(u < t) {
						// gap in d
						seq_d.append(k, '-');
						while(k--) seq_a.append(1, ccNuc[sa[--a]]);
						break;
					}
				}
				if(k <= d)
				{
					t += p(a,d-k)*model.p_indel_size[k];
					if(u < t) {
						// gap in a
						seq_a.append(k, '-');
						while(k--) seq_d.append(1, ccNuc[sd[--d]]);
						break;
					}
				}
			}
		}
		//seq_a.append(1,',');
		//seq_d.append(1,',');
	}
	reverse(seq_a.begin(), seq_a.end());
	reverse(seq_d.begin(), seq_d.end());
}

// Draw from Zipf distribution, with parameter a > 1.0
// Devroye Luc (1986) Non-uniform random variate generation.
//     Springer-Verlag: Berlin. p551
inline unsigned int rand_zipf(double a)
{
	double b = pow(2.0, a-1.0);
	double x,t;
	do {
	 x = floor(pow(myrand.uniform01(), -1.0/(a-1.0)));
	 t = pow(1.0+1.0/x, a-1.0);
	} while( myrand.uniform01()*x*(t-1.0)*b >= t*(b-1.0));
	return (unsigned int)x;
}

const char g_nuc[] = "ACGT";
const char g_nuc2[] = "ACGT" "GTAC" "TGCA" "CATG";

void gen_sample_k2p_zeta::sample_once(std::string &seq_a, std::string &seq_d)
{
	seq_a.clear();
	seq_d.clear();
	
	const model_type &model = get_model();
	const double z = get_params()[model_k2p_zeta::pZ];
	while(1) {
		double p = myrand.uniform01();
		if(p < model.p_2h) {
			// M
			unsigned int x = myrand.uniform(4);
			seq_a.append(1, g_nuc[x]);
			p = myrand.uniform01();
			if(p < model.p_match) {
				seq_d.append(1, g_nuc2[x]);
			} else if(p < model.p_match+model.p_ts) {
				seq_d.append(1, g_nuc2[x+4]);
			} else if(p < 1.0 - 0.5*model.p_ts) {
				seq_d.append(1, g_nuc2[x+8]);
			} else {
				seq_d.append(1, g_nuc2[x+12]);
			}
		} else if(p <  model.p_2h+model.p_2g) {
			// U
			unsigned int x;
			do {
				x = rand_zipf(z);
			} while(x > nmax);
			seq_d.append(x, '-');
			while(x--) {
				seq_a.append(1, g_nuc[myrand.uniform(4)]);
			}
		} else if(p <  1.0-model.p_end) {
			// V
			unsigned int x;
			do {
				x = rand_zipf(z);
			} while(x > nmax);
			seq_a.append(x, '-');
			while(x--) {
				seq_d.append(1, g_nuc[myrand.uniform(4)]);
			}			
		} else {
			// E
			break;
		}
	}
}

void gen_sample_k2p_geo::sample_once(std::string &seq_a, std::string &seq_d)
{
	seq_a.clear();
	seq_d.clear();
	
	const model_type &model = get_model();
	const double q = 1.0/get_params()[model_k2p_geo::pQ];
	while(1) {
		double p = myrand.uniform01();
		if(p < model.p_2h) {
			// M
			unsigned int x = myrand.uniform(4);
			seq_a.append(1, g_nuc[x]);
			p = myrand.uniform01();
			if(p < model.p_match) {
				seq_d.append(1, g_nuc2[x]);
			} else if(p < model.p_match+model.p_ts) {
				seq_d.append(1, g_nuc2[x+4]);
			} else if(p < 1.0 - 0.5*model.p_ts) {
				seq_d.append(1, g_nuc2[x+8]);
			} else {
				seq_d.append(1, g_nuc2[x+12]);
			}
		} else if(p <  model.p_2h+model.p_2g) {
			// U
			unsigned int x = myrand.geometric(q);
			seq_d.append(x, '-');
			while(x--) {
				seq_a.append(1, g_nuc[myrand.uniform(4)]);
			}
		} else if(p <  1.0-model.p_end) {
			// V
			unsigned int x = myrand.geometric(q);
			seq_a.append(x, '-');
			while(x--) {
				seq_d.append(1, g_nuc[myrand.uniform(4)]);
			}
		} else {
			// E
			break;
		}
	}
}



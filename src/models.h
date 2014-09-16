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

#ifndef MODELS_H
#define MODELS_H

#include <vector>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "ccvector.h"

namespace ublas = boost::numeric::ublas;

enum enumN {nA, nC, nG, nT, nN, nSize};

/***************************************************************************
 * class model_k2p_zeta                                                    *
 ***************************************************************************/

class model_k2p_zeta
{
public:
	enum {pT, pK, pR, pZ, pA, pSize};

	typedef ublas::cc_vector<double, pSize> params_type;
	typedef double (freq_type)[4];
	
	static const freq_type k2p_freqs;

	void preallocate(size_t maxa, size_t maxd);
	void expectation_setup(const params_type &params, const freq_type &freq = k2p_freqs);
	
	double p_substitution[nSize][nSize];
	
	double p_match, p_ts, p_tv;
	double p_2h, p_2g;
	double p_end;
	
	double nuc_scale, amb_scale;

	std::vector<double> p_indel_size;
	

};

/***************************************************************************
 * class model_k2p_geo                                                     *
 ***************************************************************************/

class model_k2p_geo
{
public:
	enum {pT, pK, pR, pQ, pA, pSize};

	typedef ublas::cc_vector<double, pSize> params_type;
	typedef double (freq_type)[4];

	static const freq_type k2p_freqs;
	
	void preallocate(size_t maxa, size_t maxd);
	void expectation_setup(const params_type &params, const freq_type &freq = k2p_freqs);

	double p_substitution[nSize][nSize];
		
	double p_match, p_ts, p_tv;
	double p_2h, p_2g;
	double p_end;
	
	double p_open, p_extend;

	double nuc_scale, amb_scale;
	
	std::vector<double> p_indel_size;
};

#endif //MODELS_H


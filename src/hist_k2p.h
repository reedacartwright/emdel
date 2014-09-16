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
#ifndef HIST_K2P_H
#define HIST_K2P_H

#include "hist.h"

/***************************************************************************
 * class hist_k2p_zeta                                                     *
 ***************************************************************************/


class hist_k2p_zeta : public basic_emhist<hist_k2p_zeta, 5, 10>
{
public:
	enum varP {pT, pK, pR, pZ, pA, pSize};
	enum varX {xW, xX, xXX, xSize};
	
	typedef basic_emhist<hist_k2p_zeta, pSize, 10> base_type;
	
	typedef base_type::sequence sequence;
	typedef base_type::params_type params_type;
	typedef base_type::ex_type ex_type;
	typedef base_type::exvec_type exvec_type;
	typedef ublas::cc_vector<double, xSize> xx_type;

	hist_k2p_zeta() {}
	hist_k2p_zeta(const params_type &p) : base_type(p) { }	

	inline void zero_ex(ex_type &e) { e.clear(); }

	void preallocate(size_t maxa, size_t maxd);
	void expectation_setup(const params_type &params);

	void expectation(ex_type &ex, ex_type &ex2, const sequence &seq_a, const sequence &seq_d, size_t index) const;
	void expectation(double &ex, double &ex2, size_t gap, const sequence &seq_a, const sequence &seq_d, size_t index) const;
	
	virtual std::string state() const;

private:
	double p_substitution[base_type::nSize][base_type::nSize];
	
	std::vector<double> ll_nucs;
	
	size_t sz_height, sz_width;

	double p_match, p_ts, p_tv;
	double p_2h, p_2g;
	double p_end;

	std::vector<double> p_indel_size;

	typedef table<xx_type> xx_table;
};

/***************************************************************************
 * class hist_k2p_geo                                                      *
 ***************************************************************************/

class hist_k2p_geo  : public basic_emhist<hist_k2p_geo,5, 10>
{
public:
	enum varP {pT, pK, pR, pQ, pA, pSize};
	enum varX {xW, xX, xXX, xSize};
	
	
	typedef basic_emhist<hist_k2p_geo, pSize, 10> base_type;	
	typedef base_type::sequence sequence;
	typedef base_type::params_type params_type;
	typedef base_type::ex_type ex_type;
	typedef base_type::exvec_type exvec_type;
	typedef ublas::cc_vector<double, xSize> xx_type;	
	

	hist_k2p_geo() {}
	hist_k2p_geo(const params_type &p) : base_type(p) { }	

	inline void zero_ex(ex_type &e) { e.clear(); }

	void preallocate(size_t maxa, size_t maxd);
	void expectation_setup(const params_type &params);

	void expectation(ex_type &ex, ex_type &ex2, const sequence &seq_a, const sequence &seq_d, size_t index) const;
	void expectation(double &ex, double &ex2, size_t gap, const sequence &seq_a, const sequence &seq_d, size_t index) const;
	
	virtual std::string state() const;

private:
	double p_substitution[base_type::nSize][base_type::nSize];

	std::vector<double> ll_nucs;
	
	size_t sz_height, sz_width;

	double p_match, p_ts, p_tv;
	double p_2h, p_2g;
	double p_end;
	
	double p_open, p_extend;

	std::vector<double> p_indel_size;

	std::vector<double> cache_size;

	typedef table<xx_type> xx_table;
};

#endif // EM_HIST_K2P_H

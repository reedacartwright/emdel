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
#ifndef EM_K2P_H
#define EM_K2P_H

#include "em.h"
#include "models.h"

/***************************************************************************
 * class em_k2p_zeta                                                       *
 ***************************************************************************/


class em_k2p_zeta : public basic_em<em_k2p_zeta, model_k2p_zeta, 6>
{
public:
	enum {eW, eM, eS, eV, eG, eL, eSize};
	enum {pT, pK, pR, pZ, pA, pSize};
	
	typedef basic_em<em_k2p_zeta, model_k2p_zeta, eSize> base_type;
		
	typedef base_type::sequence sequence;
	typedef base_type::params_type params_type;
	typedef base_type::model_type model_type;
	typedef base_type::ex_type ex_type;
	typedef base_type::exvec_type exvec_type;
	
	em_k2p_zeta() {}
	em_k2p_zeta(const params_type &p) : base_type(p) { }

	void preallocate(size_t maxa, size_t maxd);
	void expectation_setup(const params_type &params);

	double expectation(ex_type &ex, const sequence &seq_a, const sequence &seq_d, size_t index) const;
	void maximization(params_type &p, const exvec_type &ev);
		
	virtual std::string state() const;

private:

	size_t sz_height, sz_width;

	std::vector<double> cache_log;

	typedef table<ex_type> aln_table;
};

/***************************************************************************
 * class em_k2p_geo                                                        *
 ***************************************************************************/

class em_k2p_geo : public basic_em<em_k2p_geo, model_k2p_geo, 6>
{
public:
	enum {pT, pK, pR, pQ, pA, pSize};
	enum {eW, eM, eS, eV, eG, eL, eSize};

	typedef basic_em<em_k2p_geo, model_k2p_geo, eSize> base_type;
	
	typedef base_type::sequence sequence;
	typedef base_type::params_type params_type;
	typedef base_type::model_type model_type;
	typedef base_type::ex_type ex_type;
	typedef base_type::exvec_type exvec_type;

	em_k2p_geo() {}
	em_k2p_geo(const params_type &p) : base_type(p) { }

	void preallocate(size_t maxa, size_t maxd);
	void expectation_setup(const params_type &params);

	double expectation(ex_type &ex, const sequence &seq_a, const sequence &seq_d, size_t index) const;

	void maximization(params_type &p, const exvec_type &ev);

	virtual std::string state() const;

private:
	
	size_t sz_height, sz_width;

	std::vector<double> cache_size;

	typedef table<ex_type> aln_table;
};

/***************************************************************************
 * class em_k2p_geo_guard                                                  *
 ***************************************************************************/

class em_k2p_geo_guard : public basic2_em<em_k2p_geo_guard, model_k2p_geo, 5>
{
public:
	enum {pT, pK, pR, pQ, pA, pSize};
	enum {eM, eS, eV, eG, eL, eSize};

	typedef basic2_em<em_k2p_geo_guard, model_k2p_geo, eSize> base_type;
	
	typedef base_type::sequence sequence;
	typedef base_type::params_type params_type;
	typedef base_type::model_type model_type;
	typedef base_type::ex_type ex_type;
	typedef base_type::exvec_type exvec_type;

	em_k2p_geo_guard() {}
	em_k2p_geo_guard(const params_type &p) : base_type(p) { }

	void preallocate(size_t maxa, size_t maxd);
	void expectation_setup(const params_type &params);

	double expectation(ex_type &ex, const sequence &seq_a, const sequence &seq_d, size_t index) const;

	void maximization(params_type &p, const exvec_type &ev);

	virtual std::string state() const;

private:
	
	size_t sz_height, sz_width;

	std::vector<double> cache_size;

	typedef table<ex_type> aln_table;
};


#endif // EM_K2P_H

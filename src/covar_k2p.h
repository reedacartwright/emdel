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

#ifndef VAR_K2P_H
#define VAR_K2P_H

#include "covar.h"
#include "models.h"

/***************************************************************************
 * class covar_k2p_zeta                                                    *
 ***************************************************************************/

class covar_k2p_zeta : public basic_covar<covar_k2p_zeta, model_k2p_zeta, 20>
{
public:
	enum {pT, pK, pR, pZ, pA, pSize};
	enum {
		eM,  eS,  eV,  eG,  eL,
		eMM, eSS, eVV, eGG, eLL,
		eMS, eMV, eMG, eML, eSV,
		eSG, eSL, eVG, eVL, eGL,
		eSize};
	enum {e1W, e1M,  e1S,  e1V,  e1G,  e1L};
	enum {e2MM, e2SS, e2VV, e2GG, e2LL};
	enum {e3MS, e3MV, e3MG, e3ML, e3SV};
	enum {e4SG, e4SL, e4VG, e4VL, e4GL};
	
	typedef basic_covar<covar_k2p_zeta, model_k2p_zeta, eSize> base_type;
	
	typedef base_type::sequence sequence;
	typedef base_type::params_type params_type;
	typedef base_type::ex_type ex_type;
	typedef base_type::model_type model_type;
	typedef base_type::exvec_type exvec_type;
	
	covar_k2p_zeta() {}
	covar_k2p_zeta(const params_type &p) : base_type(p) { }

	void preallocate(size_t maxa, size_t maxd);
	void expectation_setup(const params_type &params);

	double expectation(ex_type &ex, const sequence &seq_a, const sequence &seq_d, size_t index) const;
	void covariance(covar_type &v, const params_type &p, const exvec_type &ev);

	virtual std::string state() const;
	
protected:
	void add_expectation(ex_type &ex, const ex_type &ad) const;
		
	size_t sz_height, sz_width;

	std::vector<double> cache_log;

};

/***************************************************************************
 * class covar_k2p_geo                                                     *
 ***************************************************************************/

class covar_k2p_geo : public basic_covar<covar_k2p_geo, model_k2p_geo, 20>
{
public:
	enum {pT, pK, pR, pQ, pA, pSize};
	enum {
		eM,  eS,  eV,  eG,  eY,
		eMM, eSS, eVV, eGG, eYY,
		eMS, eMV, eMG, eMY, eSV,
		eSG, eSY, eVG, eVY, eGY,
		eSize};
	
	typedef basic_covar<covar_k2p_geo, model_k2p_geo, eSize> base_type;
	
	typedef base_type::sequence sequence;
	typedef base_type::params_type params_type;
	typedef base_type::ex_type ex_type;
	typedef base_type::model_type model_type;
	typedef base_type::exvec_type exvec_type;
	
	covar_k2p_geo() {}
	covar_k2p_geo(const params_type &p) : base_type(p) { }

	void preallocate(size_t maxa, size_t maxd);
	void expectation_setup(const params_type &params);

	double expectation(ex_type &ex, const sequence &seq_a, const sequence &seq_d, size_t index) const;
	void covariance(covar_type &v, const params_type &p, const exvec_type &ev);

	virtual std::string state() const;
	
protected:
	void add_expectation(ex_type &ex, const ex_type &ad) const;
		
	size_t sz_height, sz_width;

	std::vector<double> cache_size;

};

#endif

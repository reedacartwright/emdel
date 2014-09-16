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

#ifndef COVAR_H
#define COVAR_H

#include <iomanip>
#include <cfloat>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "task.h"
#include "ccvector.h"

namespace ublas = boost::numeric::ublas;

template<class S, class M, class E>
class covar : public emdel_task
{
public:
	typedef covar<S, M, E> this_type;
	typedef S self_type;
	typedef M model_type;
	typedef E ex_type;
	
	typedef typename model_type::params_type params_type;
	typedef std::vector<E> exvec_type;
	typedef ublas::matrix<double> covar_type;
	
	covar() {}
	covar(const params_type & p) : params(p), cov(p.size(), p.size()) { }
	virtual ~covar() { }
	
	inline void set_params(const params_type& p)
	{
		params = p;
	}
	
	virtual bool run()
	{
		run_setup();
		expectation_step();
		covariance_step();
		log_state();
		report();
		return true;
	}
	
	virtual void report() const
	{
		std::cout << std::endl << "Variance Done" << std::endl << std::endl;
	}
	
	virtual void log_state()
	{
		outlog() << std::setprecision(DBL_DIG)
			<< state() << std::endl;
	}
	
	virtual std::string state() const
	{
		return std::string();
	}	
	
	inline const params_type& get_params() const { return params; }
	inline const exvec_type& get_exvec() const { return exes; }
	inline const covar_type& get_covar() const { return cov; }
	inline const model_type& get_model() const { return model; }
	
	void run_setup()
	{
		// Find Maximum length of an ancestor and descendent
		size_t maxa = 0, maxd = 0;
		for(typename seq_pair_vector::const_iterator cit = get_seqs().begin();
		    cit != get_seqs().end(); ++cit)
		{
			maxa = std::max(cit->first.size(), maxa);
			maxd = std::max(cit->second.size(), maxd);
		}
		
		exes.resize(get_seqs().size());
		
		model.preallocate(maxa, maxd);
		self().preallocate(maxa, maxd);
	}
	
	void expectation_step()
	{
		reset_queue();
		model.expectation_setup(params);
		self().expectation_setup(params);
		expectations_parralel();
	}

	void expectations_parralel()
	{
		ex_type temp;
		seq_pair_vector::const_iterator cit;
		while(next_in_queue(cit))
		{
			boost::timer t;
			size_t pos = (size_t)(cit-sequence_pairs.begin());
			output() << "Pair " << pos << " (" << cit->first.size() << "x" << cit->second.size() << ")..." << std::flush;
			double w = self().expectation(temp, cit->first, cit->second, pos);
			output() << " ll = " << w << " [" << std::setprecision(3) << t.elapsed() << " s]" << std::endl;
			exes[pos] = temp;
		}
	}
		
	void covariance_step()
	{
		self().covariance(cov, params, exes);
	}	

protected:
	self_type & self() { return *static_cast<self_type *>(this); }
	const self_type & self() const { return *static_cast<const self_type *>(this); }

	void reporter()
	{
	}
	
private:
	params_type params;
	exvec_type exes;
	model_type model;
	covar_type cov;
	
};

template<class T, class M, std::size_t E>
class basic_covar : public covar<T, M,
	ublas::cc_vector<double, E> >
{
public:
	typedef covar<T, M, ublas::cc_vector<double, E> > base_type;
	
	typedef typename base_type::params_type params_type;
	typedef typename base_type::ex_type ex_type;
	
	basic_covar() { }
	basic_covar(const params_type &p) : base_type(p) { }
	
	inline void zero_ex(ex_type &e) { e.clear(); }
	
};


#endif

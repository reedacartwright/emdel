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

#ifndef EMHIST_H
#define EMHIST_H

#include <iomanip>
#include <iostream>

#include "task.h"

template<class S, class P, class E>
class emhist : public emdel_task
{
public:
	typedef emhist<S, P, E> this_type;
	typedef S self_type;
	typedef P params_type;
	typedef E ex_type;
	typedef std::vector<E> exvec_type;
	
	emhist() {}
	emhist(const params_type & p) : params(p) { }
	virtual ~emhist() { }
	
	
	inline void set_params(const params_type& p)
	{
		params = p;
	}
	
	virtual bool run()
	{
		run_setup();
		histogram_step();
		log_state();
		report();
		return true;
	}
	
	virtual void report() const
	{
		std::cout << std::endl << "Histogram Done" << std::endl << std::endl;
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
	
	const params_type& get_params() const { return params; }
	
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
		ex2s.resize(get_seqs().size());
		
		self().preallocate(maxa, maxd);
	}
	
	void histogram_step()
	{
		reset_queue();
		self().expectation_setup(params);
		histogram_parralel();
	}
	
	void histogram_parralel()
	{
		ex_type temp;
		ex_type temp2;
		seq_pair_vector::const_iterator cit;
		
		while(next_in_queue(cit))
		{
			size_t pos = (size_t)(cit-sequence_pairs.begin());
			
			self().expectation(temp, temp2, cit->first, cit->second, pos);
			exes[pos] = temp;
			ex2s[pos] = temp2;
		}
	}
	
protected:
	self_type & self() { return *static_cast<self_type *>(this); }
	const self_type & self() const { return *static_cast<const self_type *>(this); }

	void reporter()
	{
	}
	
private:
	
	params_type params;
protected:
	exvec_type exes;
	exvec_type ex2s;
};

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "ccvector.h"
#include "table.h"

namespace ublas = boost::numeric::ublas;

template<class T, std::size_t P, std::size_t E>
class basic_emhist : public emhist<T,
	ublas::cc_vector<double, P>,
	ublas::cc_vector<double, E>
	>
{
public:
	typedef emhist<T, ublas::cc_vector<double, P>,
		ublas::cc_vector<double, E> >	base_type;
	
	typedef typename base_type::params_type params_type;
	typedef typename base_type::ex_type ex_type;
	
	basic_emhist() { }
	basic_emhist(const params_type &p) : base_type(p) { }
};

#endif

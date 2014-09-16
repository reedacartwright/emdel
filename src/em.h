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
#ifndef EM_H
#define EM_H

#include <cmath>
#include <vector>
#include <string>
#include <utility>
#include <iostream>
#include <iomanip>
#include <cfloat>

#include "task.h"

#ifdef _MSC_VER
namespace std {
	inline bool isfinite(double x) { return _finite(x) != 0; }
}
#endif

/***************************************************************************
 * template class em_base                                                  *
 ***************************************************************************/

class em_base : public emdel_task
{
public:
	em_base() : steps(50), threshold(1e-5) { }
	
	inline void set_threshold(double t) { threshold = t; }
	inline void set_steps(int s) { steps = s; }

protected:
	double threshold;
	int steps;	
};

/***************************************************************************
 * template class em                                                       *
 ***************************************************************************/

template<class S, class M, class E>
class em : public em_base
{
public:
	typedef em<S, M, E> this_type;
	typedef S self_type;
	typedef M model_type;
	typedef E ex_type;
	
	typedef typename model_type::params_type params_type;
	typedef std::vector<E> exvec_type;

	em() {}
	em(const params_type & p) : params(p) { }

	inline void set_initial(const params_type& p)
	{
		params = p;
	}
			
	inline bool has_converged() const
	{
		return (loglh-loglh_old < threshold && loglh >= loglh_old);
	}

	virtual bool run()
	{
		run_setup();
		for(int i=1;i<=steps;++i)
		{
			output() << std::endl
				<< "[==== Step " << i << " =======================================================]"
				<< std::endl
				;
			
			expectation_step();
			maximization_step();
			log_state(i);
			report(i);
			if(has_converged())
				return true;
			if(!std::isfinite(loglh))
				return false;
		}
		return false;
	}
	
	virtual void report(int step)
	{
		output() << std::setprecision(DBL_DIG)
			<< "  Log Likelihood: " << loglh
			<< " (change " << ((step > 1) ? loglh-loglh_old : 0.0) << ")"
			<< std::endl
			;
	}
	
	virtual void log_state(int step)
	{
		outlog() << std::setprecision(DBL_DIG)
			<< "[==== Step " << step << " ====]" << std::endl
			<< "  Log-Likelihood: " << loglh
			<< " (change " << ((step > 1) ? loglh-loglh_old : 0.0) << ")"
			<< std::endl << std::endl
			<< state() << std::endl
			;
	}
	
	virtual std::string state() const
	{
		return std::string();
	}

	inline const params_type& get_params() const { return params; }
	inline const exvec_type& get_exvec() const { return exes; }
	inline double get_loglh() const { return loglh; }
	inline const model_type& get_model() const { return model; }
	
protected:
	inline self_type & self() { return *static_cast<self_type *>(this); }
	inline const self_type & self() const { return *static_cast<const self_type *>(this); }
	
	void run_setup()
	{
		loglh = -DBL_MAX;
		loglh_old = loglh;
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
		
		loglh_old = loglh;
		loglh = 0.0;

#ifdef WITH_PARALLEL
		boost::thread report(boost::bind(&this_type::reporter, this));
		boost::thread_group workers;
		for (int i = 0; i < num_threads; ++i)
			workers.create_thread(boost::bind(&this_type::expectations_parralel, this));
		workers.join_all();
		
		stop_reporter.notify_all();
		report.join();
#else
		expectations_parralel();
#endif
	}

	void maximization_step()
	{
		params_old = params;
		self().maximization(params, exes);
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
#ifdef WITH_PARALLEL
			boost::mutex::scoped_lock exlock(ex_mutex);
#endif
			loglh += w;
			exes[pos] = temp;
		}
	}

private:
	
	double loglh, loglh_old;
	
	params_type params, params_old;
	exvec_type exes;
	model_type model;

};

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "ccvector.h"
#include "table.h"
#include "big_prob.h"

namespace ublas = boost::numeric::ublas;

template<class T, class M, std::size_t E>
class basic_em : public em<T, M,
	ublas::cc_vector<double, E>
	>
{
public:
	typedef em<T, M,
		ublas::cc_vector<double, E>
		>	base_type;
	
	typedef typename base_type::params_type params_type;
	typedef typename base_type::ex_type ex_type;
	
	basic_em() { }
	basic_em(const params_type &p) : base_type(p) { }
	
	inline void zero_ex(ex_type &e) { e.clear(); }
	
};

template<std::size_t E>
struct bex_vector {
	typedef bex_vector<E> this_type;
	typedef ublas::cc_vector<float, E> stat_type;
	
	bex_vector() : p(0.0), stat(ublas::zero_vector<float>(E)) { }

	stat_type stat;
	big_prob p;
	
	inline this_type& operator *=(const this_type& r) {
		noalias(stat) += r.stat;
		p *= r.p;
		return *this;
	}
	inline this_type& operator +=(const this_type& r) {
		p += r.p;
		float s = r.p/p;
		//stat *= s;
		//noalias(stat) += (1.0-s)*r.stat;
		//std::cerr << s << std::endl;
		//noalias(stat) = (1.0f-s)*stat + s*r.stat;
		noalias(stat) += s*(r.stat-stat);
		return *this;
	}
	inline this_type operator *(const this_type& r) const {
		this_type ret(*this);
		return (ret *= r);
	}
	inline this_type operator +(const this_type& r) const {
		this_type ret(*this);
		return (ret += r);
	}
	
	inline this_type& operator *=(const big_prob & r) {
		p *= r;
		return *this;
	}
	inline this_type operator *(const big_prob & r) const {
		this_type ret(*this);
		return (ret *= r);
	}
	
	inline this_type& plus2(const this_type& r1, const this_type& r2) {
		big_prob newp = p+r1.p+r2.p;
		float s1 = r1.p/newp;
		float s2 = r2.p/newp;
		//stat *= s;
		//noalias(stat) += (1.0-s)*r.stat;
		noalias(stat) = (1.0f-s1-s2)*stat + s1*r1.stat + s2*r2.stat;
		p = newp;
		return *this;
	}	
	
	inline void clear() {
		p = 0.0;
		stat.clear();
	}
};

template<class T, class M, std::size_t E>
class basic2_em : public em<T, M,
	bex_vector<E>
	>
{
public:
	typedef em<T, M,
		bex_vector<E>
		>	base_type;
	
	typedef typename base_type::params_type params_type;
	typedef typename base_type::ex_type ex_type;
	
	basic2_em() { }
	basic2_em(const params_type &p) : base_type(p) { }
	
	inline void zero_ex(ex_type &e) { e.clear(); }
	
};

#endif

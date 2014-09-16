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

#ifndef SAMPLE_H
#define SAMPLE_H

#include <iomanip>
#include <cfloat>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "task.h"
#include "ccvector.h"

namespace ublas = boost::numeric::ublas;


struct clustal_aln : public std::pair<std::string,std::string> {
	typedef std::string ss;
	typedef std::pair<ss,ss> base_type;
	clustal_aln(const ss &a, const ss &b) : base_type(a,b) { }
	
	template<class OS>
	inline void print(OS &os) const {
		os << "CLUSTAL multiple sequence alignment (Created by " << PACKAGE_STRING;
		//if(msg != NULL)
		//	os << "; " << msg;
		os << ")" << std::endl << std::endl;
	
		size_t sz = first.size();
		// Print interleaved sequences
		size_t a=0, b=0;
		std::string ss;
		std::string nameA("Seq_A"), nameB("Seq_D");
		for(size_t u = 0; u < sz; u+=60)
		{
			// Print a row of each sequence
			ss = first.substr(u, 60);
			a += ss.length() - std::count(ss.begin(), ss.end(), '-');
			os << std::setw(14) << std::setiosflags(std::ios::left) << nameA.substr(0, 14)
				<< " " << ss << " " << a << std::endl;
			
			ss = second.substr(u, 60);
			b += ss.length() - std::count(ss.begin(), ss.end(), '-');
			os << std::setw(14) << std::setiosflags(std::ios::left) << nameB.substr(0, 14)
				<< " " << ss << " " << b << std::endl;
			
			//os << std::setw(14) << std::setiosflags(std::ios::left) << " "   << " " << strC.substr(u, 60) << std::endl;
			os << std::endl;
		}
	}
};

template<class OS>
OS& operator<<(OS &os, const clustal_aln &aln) {
	aln.print(os);
	return os;
}

struct fasta_aln : public std::pair<std::string,std::string> {
	typedef std::string ss;
	typedef std::pair<ss,ss> base_type;
	fasta_aln(const ss &a, const ss &b, unsigned int x) : base_type(a,b), n(x) { }
	
	template<class OS>
	inline void print(OS &os) const {
		size_t sz = first.size();
		// Print sequences
		std::string ss;
		std::string nameA("Seq_A"), nameB("Seq_D");
		os << ">" << nameA << "_" << n << std::endl;
		for(size_t u = 0; u < sz; u+=60) {
			// Print a row of each sequence
			ss = first.substr(u, 60);
			os << ss << std::endl;
		}
		os << std::endl;
		os << ">" << nameB << "_" << n << std::endl;
		for(size_t u = 0; u < sz; u+=60) {
			// Print a row of each sequence
			ss = second.substr(u, 60);
			os << ss << std::endl;
		}
		os << std::endl;
	}
	unsigned int n;
};

template<class OS>
OS& operator<<(OS &os, const fasta_aln &aln) {
	aln.print(os);
	return os;
}

class sample_base : public emdel_task
{
public:
	sample_base() : nsamples(1) { }
	
	inline void set_samples(unsigned int n) { nsamples = n; }
	inline void set_max(unsigned int n) { nmax = n; }
	inline void set_min(unsigned int n) { nmin = n; }
	
protected:
	unsigned int nsamples;
	unsigned int nmax;
	unsigned int nmin;
};

template<class S, class M>
class sample : public sample_base
{
public:
	typedef sample<S, M> this_type;
	typedef S self_type;
	typedef M model_type;
	
	typedef typename model_type::params_type params_type;
	//typedef emdel_task::sequence sequence;
	
	sample() {}
	sample(const params_type & p) : params(p) { }
	virtual ~sample() { }
	
	inline void set_params(const params_type& p) {
		params = p;
	}
		
	virtual bool run()	{
		run_setup();
		do_sample();
		return true;
	}
		
	inline const params_type& get_params() const { return params; }
	inline const model_type& get_model() const { return model; }
	
	void run_setup() {
		// Find length of first ancestor and descendent
		size_t maxa = get_seqs().front().first.size();
		size_t maxd = get_seqs().front().second.size();
		model.preallocate(maxa, maxd);
		self().preallocate(maxa, maxd);
	}
	
	void presample_step(const sequence &seq_a, const sequence &seq_d) {
		model.expectation_setup(params);
		self().presample_step(params, seq_a, seq_d);
	}
	
	void do_sample() {
		std::cerr << "Presampling..." << std::endl;
		presample_step(get_seqs().front().first, get_seqs().front().second);
		std::cerr << "Sampling... " << nsamples << std::endl;
		std::string sa, sd;
		for(unsigned int i=0;i<nsamples;++i) {
			self().sample_once(sa, sd);
			std::cout << sa << std::endl << sd << std::endl << std::endl;
		}
		
	}


protected:
	self_type & self() { return *static_cast<self_type *>(this); }
	const self_type & self() const { return *static_cast<const self_type *>(this); }
	
private:
	params_type params;
	model_type model;
};

template<class S, class M>
class gen_sample : public sample_base
{
public:
	typedef gen_sample<S, M> this_type;
	typedef S self_type;
	typedef M model_type;
	
	typedef typename model_type::params_type params_type;
	//typedef emdel_task::sequence sequence;
	
	gen_sample() {}
	gen_sample(const params_type & p) : params(p) { }
	virtual ~gen_sample() { }
	
	inline void set_params(const params_type& p) {
		params = p;
	}
		
	virtual bool run() {
		run_setup();
		do_sample();
		return true;
	}
	
	inline const params_type& get_params() const { return params; }
	inline const model_type& get_model() const { return model; }
	
	void run_setup() {
		model.preallocate(1, 1);
	}
	
	void presample_step() {
		model.expectation_setup(params);
	}
	
	void do_sample() {
		presample_step();
		std::cerr << "Generating Random Alignments... " << std::endl;
		std::string sa, sd;
		for(unsigned int i=0;i<nsamples;++i) {
			size_t na;
			size_t nd;
			do {
				self().sample_once(sa, sd);
				na = sa.length() - std::count(sa.begin(), sa.end(), '-');
				nd = sd.length() - std::count(sd.begin(), sd.end(), '-');
			} while ( nmin > na || na > nmax || nmin > nd || nd > nmax);
			std::cout << fasta_aln(sa,sd,i) << std::endl;
		}
		
	}


protected:
	self_type & self() { return *static_cast<self_type *>(this); }
	const self_type & self() const { return *static_cast<const self_type *>(this); }
	
private:
	params_type params;
	model_type model;
};


#endif

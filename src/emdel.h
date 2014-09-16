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

#ifndef EMDEL_H
#define EMDEL_H

#include <boost/config.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "em.h"

class emdel_app  {
public:
	emdel_app(int argc, char *argv[]);
	virtual ~emdel_app() { }
	
	virtual int run();
	
	std::string args_comment() const;

	po::options_description desc;	
	
	// use X-Macros to specify argument variables
	struct arg_t {
#	define XCMD(lname, sname, desc, type, def) type V(lname) ;
#	include "emdel.cmds"
#	undef XCMD
	};
protected:
	arg_t arg;

};

#include <ostream>
#include <iterator>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace std {
template<typename _Tp, typename _CharT, typename _Traits>
basic_ostream<_CharT, _Traits>&
operator<<(basic_ostream<_CharT, _Traits>& os, const std::vector<_Tp> &v)
{
	if(v.size() == 1)
	{
		os << v.front();
	}
	else if(v.size() > 1)
	{
		std::copy(v.begin(), v.end()-1, std::ostream_iterator<_Tp, _CharT, _Traits>(os, " "));
		os << v.back();
	} 
	return os;
}
}

#ifdef BOOST_WINDOWS
#	include <process.h>
#	define getpid _getpid
#endif

inline unsigned int rand_seed_start()
{
	unsigned int u = static_cast<unsigned int>(time(NULL));
	unsigned int v = static_cast<unsigned int>(getpid());
	v += (v << 15) + (v >> 3); // Spread 5-decimal PID over 32-bit number
	return u ^ (v + 0x9e3779b9u + (u<<6) + (u>>2));
}

inline unsigned int rand_seed()
{
	static unsigned int u = rand_seed_start();
	return (u = u*1664525u + 1013904223u);
}

class gslrand
{
public:
	typedef unsigned int int_type;
	typedef int_type seed_type;
	typedef double real_type;
	gslrand() {
		//gsl_rng_env_setup();
		T = gsl_rng_mt19937;
		r = gsl_rng_alloc(T);
	}
	virtual ~gslrand() {
		gsl_rng_free(r);
	}
	
	inline void seed(seed_type s) { gsl_rng_set(r, s); }
	inline int_type get() { return gsl_rng_get(r); }
	inline int_type operator()() { return get(); }
	
	inline int_type uniform() { return get(); }
	inline int_type uniform(int_type max) { return gsl_rng_uniform_int(r, max); }
	inline real_type uniform01() { return gsl_rng_uniform(r); }
	inline int_type geometric(real_type p) { return gsl_ran_geometric(r,p); }

private:
	const gsl_rng_type *T;
	gsl_rng *r;
};

extern gslrand myrand;


void set_rand_seed(unsigned int u);


#endif

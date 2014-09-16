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

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

//#define _HAS_ITERATOR_DEBUGGING 0
//#define _SECURE_SCL 0
#define _USE_MATH_DEFINES

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>
#include <fstream>
#include <xmmintrin.h>

#include "emdel.h"
#include "seqdb.h"
#include "em_k2p.h"
#include "covar_k2p.h"
#include "hist_k2p.h"
#include "sample_k2p.h"

using namespace std;

const char emdel_task::ccNuc[] = "ACGTN";

int main(int argc, char *argv[])
{
	_mm_setcsr(_mm_getcsr() | 0x8840);

	int ret = EXIT_FAILURE;
	try {
		emdel_app app(argc, argv);
		ret = app.run();
	} catch(exception &e) {
		cerr << "ERROR: " << e.what() << endl;
	}
	return ret;
}

namespace boost { namespace program_options {
template<>
typed_value<bool>* value() { return bool_switch(); }

template<>
typed_value<bool>* value(bool* v) { return bool_switch(v); }

}}

emdel_app::emdel_app(int argc, char* argv[]) : desc("Allowed Options")
{
	try {
		desc.add_options()
			#define XCMD(lname, sname, desc, type, def) ( \
				S(lname) IFD(sname, "," BOOST_PP_STRINGIZE sname), \
				po::value< type >(&arg.V(lname))->default_value(def), \
				desc )
			#include "emdel.cmds"
			#undef XCMD
			;
		po::variables_map vm;
		po::positional_options_description pdesc;
		pdesc.add("input", -1);
		po::store(po::command_line_parser(argc, argv).options(desc).positional(pdesc).run(), vm);
		po::notify(vm);
		if(!arg.read_args.empty())
		{
			if(arg.read_args == "-")
			{
				po::store(po::parse_config_file(cin, desc), vm);	
			}
			else
			{
				std::ifstream ifs(arg.read_args.c_str());
				if(!ifs.is_open())
				{
					string sse = "unable to open argument file ";
					sse += arg.read_args;
					throw std::runtime_error(sse);
				}
				po::store(po::parse_config_file(ifs, desc), vm);
			}
			po::notify(vm);
		}
	} catch (exception &e) {
			cerr << "ERROR: " << e.what() << endl;
			throw std::runtime_error("unable to process command line");
	}
}

template<typename T>
inline bool do_print_args(const T&) { return true;}

template<>
inline bool do_print_args(const std::string& s) {
	return !s.empty();
}

std::string emdel_app::args_comment() const
{
	ostringstream ostr;
	ostr << std::setprecision(DBL_DIG);
	#define XCMD(lname, sname, desc, type, def) \
	if(do_print_args(arg.V(lname))) \
		ostr << S(lname) << "=" << arg.V(lname) << std::endl;
	#include "emdel.cmds"
	#undef XCMD
		
	return ostr.str();
}

int emdel_app::run()
{
	if(arg.help || arg.version)
	{
		cerr << desc << endl;
		return EXIT_SUCCESS;
	}
	ios_base::openmode om = (arg.append) ? 
		ios_base::out|ios_base::app : 
		ios_base::out|ios_base::trunc ;
	
	ofstream olog(arg.log.c_str(), om);
	if(!olog.is_open())
	{
		cerr << "ERROR: unable to open log file \'" << arg.log << "\'." << endl;
		return EXIT_FAILURE;
	}
	while(arg.seed == 0)
		arg.seed = rand_seed();
	set_rand_seed(arg.seed);
	
	time_t rawtime;
	time(&rawtime);
	string sstime = ctime(&rawtime);
	sstime.erase(24, 1);
	olog << "############################################################" << endl;
	olog << "# EMDEL v1.0                                               #" << endl;
	olog << "# " << sstime         << "                                 #" << endl;
	olog << "############################################################" << endl;
		
	olog << "RUN OPTIONS:"  << endl;
	olog << args_comment() << endl;
	emdel_task *ptask = NULL;
	if(arg.model == "zeta")
	{
		model_k2p_zeta::params_type ini;
		ini[model_k2p_zeta::pT] = arg.branch_length;
		ini[model_k2p_zeta::pK] = arg.ratio;
		ini[model_k2p_zeta::pR] = arg.indel_rate;
		ini[model_k2p_zeta::pZ] = arg.indel_slope;
		ini[model_k2p_zeta::pA] = arg.avgaln;
		
		if(arg.task == "em")
			ptask = new em_k2p_zeta(ini);
		else if(arg.task == "covar")
			ptask = new covar_k2p_zeta(ini);
		else if(arg.task == "hist")
			ptask = new hist_k2p_zeta(ini);
		else if(arg.task == "sample")
			ptask = new sample_k2p_zeta(ini);
		else if(arg.task == "gen")
			ptask = new gen_sample_k2p_zeta(ini);
		else
		{
			cerr << "ERROR: unknown task mode: '" << arg.task << "'.  Options are 'em', 'covar', 'local', 'hist', and 'sample'." << endl;
			return EXIT_FAILURE;
		}
	}
	else if(arg.model == "geo")
	{
		model_k2p_geo::params_type ini;
		ini[model_k2p_geo::pT] = arg.branch_length;
		ini[model_k2p_geo::pK] = arg.ratio;
		ini[model_k2p_geo::pR] = arg.indel_rate;
		ini[model_k2p_geo::pQ] = arg.indel_mean;
		ini[model_k2p_geo::pA] = arg.avgaln;
		
		if(arg.task == "em")
			ptask = new em_k2p_geo(ini);
		else if(arg.task == "covar")
			ptask = new covar_k2p_geo(ini);
		else if(arg.task == "hist")
			ptask = new hist_k2p_geo(ini);
		else if(arg.task == "sample")
			ptask = new sample_k2p_geo(ini);
		else if(arg.task == "gen")
			ptask = new gen_sample_k2p_geo(ini);
		else {
			cerr << "ERROR: unknown task mode: '" << arg.task << "'.  Options are 'em', 'covar', 'local', 'hist', and 'sample'." << endl;
			return EXIT_FAILURE;
		}
	}
	else if(arg.model == "geo-guard")
	{
		model_k2p_geo::params_type ini;
		ini[model_k2p_geo::pT] = arg.branch_length;
		ini[model_k2p_geo::pK] = arg.ratio;
		ini[model_k2p_geo::pR] = arg.indel_rate;
		ini[model_k2p_geo::pQ] = arg.indel_mean;
		ini[model_k2p_geo::pA] = arg.avgaln;
		
		if(arg.task == "em")
			ptask = new em_k2p_geo_guard(ini);
		else {
			cerr << "ERROR: unknown task mode: '" << arg.task << "'.  Options are 'em', 'covar', 'local', 'hist', and 'sample'." << endl;
			return EXIT_FAILURE;
		}
	}
	else {
		cerr << "ERROR: unknown evolutionary model: '" << arg.model << "'.  Options are 'zeta' and 'geo'." << endl;
		return EXIT_FAILURE;
	}
	
	if(arg.task == "em") {
		static_cast<em_base*>(ptask)->set_steps(arg.cycles);
		static_cast<em_base*>(ptask)->set_threshold(arg.threshold);
	} else if(arg.task == "sample" || arg.task == "gen") {
		static_cast<sample_base*>(ptask)->set_samples(arg.sample_size);
		static_cast<sample_base*>(ptask)->set_max(arg.sample_max);
		static_cast<sample_base*>(ptask)->set_min(arg.sample_min);
	}
	
	ptask->set_scale(arg.scale);
	ptask->set_log(&olog);
	if(arg.task != "gen") {
		seq_db mydb;
		for(vector<string>::const_iterator cit = arg.input.begin(); cit != arg.input.end(); ++cit)
		{
			if(!mydb.parse_file(cit->c_str(), true))
			{
				cerr << "ERROR: parsing of \'" << cit->c_str() << "\' failed." << endl;
				return EXIT_FAILURE;
			}
		}
	
		for(size_t u = 0; u+1 < mydb.size(); u+=2) {
			if(arg.seq_max == 0 || ( mydb[u].second.size() < arg.seq_max
					&& mydb[u+1].second.size()  < arg.seq_max ))
				ptask->add_pair( mydb[u].second, mydb[u+1].second );
		}
	}
	ptask->run();
	
	cout << endl << endl;
	
	delete ptask;
		
	return EXIT_SUCCESS;
}

gslrand myrand;

void set_rand_seed(unsigned int u) { myrand.seed(u); }

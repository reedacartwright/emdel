
#ifndef SAMPLE_K2P_H
#define SAMPLE_K2P_H

#include "sample.h"
#include "models.h"
#include "table.h"

/***************************************************************************
 * class sample_k2p_zeta                                                   *
 ***************************************************************************/

class sample_k2p_zeta : public sample<sample_k2p_zeta, model_k2p_zeta>
{
public:
	
	typedef sample<sample_k2p_zeta, model_k2p_zeta> base_type;
	
	typedef base_type::sequence sequence;
	typedef base_type::params_type params_type;
	typedef base_type::model_type model_type;
	
	sample_k2p_zeta() {}
	sample_k2p_zeta(const params_type &p) : base_type(p) { }
	virtual ~sample_k2p_zeta() { }

	void preallocate(size_t maxa, size_t maxd);
	void presample_step(const params_type &params, const sequence &seq_a, const sequence &seq_d);

	void sample_once(std::string &seq_a, std::string &seq_d);
		
protected:
		
	size_t sz_height, sz_width;

	sequence sa, sd;
	
	typedef table<double> aln_table;
	aln_table p;

};

/***************************************************************************
 * class sample_k2p_geo                                                    *
 ***************************************************************************/

class sample_k2p_geo : public sample<sample_k2p_geo, model_k2p_geo>
{
public:
	
	typedef sample<sample_k2p_geo, model_k2p_geo> base_type;
	
	typedef base_type::sequence sequence;
	typedef base_type::params_type params_type;
	typedef base_type::model_type model_type;
	
	sample_k2p_geo() {}
	sample_k2p_geo(const params_type &p) : base_type(p) { }
	virtual ~sample_k2p_geo() { }

	void preallocate(size_t maxa, size_t maxd);
	void presample_step(const params_type &params, const sequence &seq_a, const sequence &seq_d);

	void sample_once(std::string &seq_a, std::string &seq_d);
		
protected:
		
	size_t sz_height, sz_width;

	sequence sa, sd;
	
	typedef table<double> aln_table;
	aln_table p;

};


/***************************************************************************
 * class gen_sample_k2p_zeta                                               *
 ***************************************************************************/

class gen_sample_k2p_zeta : public gen_sample<gen_sample_k2p_zeta, model_k2p_zeta>
{
public:
	
	typedef gen_sample<gen_sample_k2p_zeta, model_k2p_zeta> base_type;
	
	typedef base_type::sequence sequence;
	typedef base_type::params_type params_type;
	typedef base_type::model_type model_type;
	
	gen_sample_k2p_zeta() {}
	gen_sample_k2p_zeta(const params_type &p) : base_type(p) { }
	virtual ~gen_sample_k2p_zeta() { }

	void sample_once(std::string &seq_a, std::string &seq_d);
		
protected:

};

/***************************************************************************
 * class gen_sample_k2p_geo                                                *
 ***************************************************************************/

class gen_sample_k2p_geo : public gen_sample<gen_sample_k2p_geo, model_k2p_geo>
{
public:
	
	typedef gen_sample<gen_sample_k2p_geo, model_k2p_geo> base_type;
	
	typedef base_type::sequence sequence;
	typedef base_type::params_type params_type;
	typedef base_type::model_type model_type;
	
	gen_sample_k2p_geo() {}
	gen_sample_k2p_geo(const params_type &p) : base_type(p) { }
	virtual ~gen_sample_k2p_geo() { }

	void sample_once(std::string &seq_a, std::string &seq_d);
		
protected:

};



#endif

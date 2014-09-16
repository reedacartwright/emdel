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

#ifndef EMDEL_TASK_H
#define EMDEL_TASK_H

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>

#include <boost/timer.hpp>

/***************************************************************************
 * class emdel_task                                                        *
 ***************************************************************************/

class emdel_task
{
public:
	typedef std::vector<char> sequence;
	typedef std::pair<sequence, sequence> seq_pair;
	typedef std::vector<seq_pair> seq_pair_vector;
	typedef std::vector<sequence::size_type> seq_info;
	typedef std::vector<seq_info> seq_pair_info_vector;
	
	enum {nA, nC, nG, nT, nN, nSize};
	static const char ccNuc[];
	
	emdel_task() : num_threads(1), mark_step(60), plog(NULL), total_nucs(0) { }
	virtual ~emdel_task() { }
	
	inline void set_threads(int i)
	{
		num_threads = (i <= 0) ? 1 : i;
	}
	
	inline void set_time_step(int step)
	{
		mark_step = step;
	}
	
	inline void set_scale(double s)
	{
		scale = s;
		prob_scale = pow(2.0, s);
	}
	
 	inline void set_log(std::ostream *p)
 	{
 		plog = p;
 	}
	
	inline std::ostream& outlog() const { return (plog == NULL) ? std::cout : *plog; }
	inline std::ostream& output() const { return std::cout; }
		
	inline const seq_pair_vector& get_seqs() const { return sequence_pairs; }
	inline const seq_pair_info_vector& get_seq_info() const { return sequence_pair_info; }
	
	bool add_pair(const std::string& a, const std::string &d)
	{
		sequence seqA, seqD;
		if(!emdel_task::process_sequence(a, seqA))
			return false;
		if(!emdel_task::process_sequence(d, seqD))
			return false;
		sequence_pairs.push_back(make_pair(seqA, seqD));
		
		seq_info nfo(nSize, 0);
		for(sequence::const_iterator nit = seqA.begin(); nit != seqA.end(); ++nit)
			++nfo[static_cast<sequence::size_type>(*nit)];
		for(sequence::const_iterator nit = seqD.begin(); nit != seqD.end(); ++nit)
			++nfo[static_cast<sequence::size_type>(*nit)];
		sequence_pair_info.push_back(nfo);
		
		total_nucs += seqA.size() + seqD.size();
		
		return true;
	}
	
	virtual bool run() = 0;
	
	static bool process_sequence(const std::string& in, sequence& out)
	{
		static const int nuc_table[] = {
	/*64*/	-1, nA, nN, nC, nN, -1, -1, nG, nN, -1, -1, nN, -1, nN, nN, -1,
	/*80*/	-1, -1, nN, nN, nT, nT, nN,	nN, -1, nN, -1, -1, -1, -1, -1, -1,
	/*96*/	-1, nA, nN, nC, nN, -1, -1, nG,	nN, -1, -1, nN, -1, nN, nN, -1,
	/*112*/	-1, -1, nN, nN, nT, nT, nN,	nN, -1, nN, -1, -1, -1, -1, -1, -1,
		};

		out.clear();
		out.reserve(in.length());
		int temp;
		for(std::string::const_iterator cit = in.begin(); cit != in.end(); ++cit)
		{

			if( *cit >= 64 && (temp = nuc_table[*cit-64]) != -1)
				out.push_back((char)temp);
		}
		return true;
	}
	
	inline void reset_queue()
	{
		qit = sequence_pairs.begin();
	}
	
	inline bool next_in_queue(seq_pair_vector::const_iterator &cit)
	{
#ifdef WITH_PARALLEL
		boost::mutex::scoped_lock qlock(q_mutex);
#endif
		if(qit == sequence_pairs.end())
			return false;
		cit = qit++;
		return true;
	}	
	
	
protected:
	seq_pair_vector sequence_pairs;
	seq_pair_info_vector sequence_pair_info;
	sequence::size_type total_nucs;
	
	int num_threads;
	int mark_step;
	
	double prob_scale;
	double scale;
	
private:
	seq_pair_vector::const_iterator qit;
	std::ostream *plog;
};

template<typename T>
struct sum_elements
{
	T sum;
	sum_elements(T t) : sum(t) { }
	sum_elements() : sum(static_cast<T>(0)) { }
	void operator()(const T& e) {sum += e;}
};

template<typename T>
struct prod_elements {
	T prod;
	prod_elements(T t) : prod(t) { }
	prod_elements() : prod(static_cast<T>(1)) { }
	void operator()(const T& e) {prod *= e;}
};

#endif

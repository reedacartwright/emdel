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
 
#include <boost/spirit/include/classic_core.hpp>
#include <boost/spirit/include/classic.hpp>
#include <iostream>
#include <stack>

#include "seqdb.h"

bool parse_file(const char* cs, seq_db &rdb);

struct push_string
{
	push_string(std::stack<std::string>& st) : my_stack(st) { }
	
	template<typename It> void operator() ( It first, It last) const
	{
		my_stack.push(std::string(first, last));
		//std::cout << "Pushed " << std::string(first, last) << std::endl;
	}
		
	std::stack<std::string>& my_stack;
};

inline void trim_string(std::string &ss, const char *ch = " \n\t\r\v\f")
{
	std::string::size_type sz = ss.find_last_not_of(ch);
	if(sz != std::string::npos) ss.erase(sz+1); 
	sz = ss.find_first_not_of(ch);
	if(sz != std::string::npos) ss.erase(0, sz);
}

inline void sanitize_string(std::string &ss, const char *ch = " \n\t\r\v\f")
{
	std::string::size_type sz = ss.find_first_of(ch);
	while(sz != std::string::npos)
	{
		ss.erase(sz, 1);
		sz = ss.find_first_of(ch, sz);
	}
}

struct pop_sequence
{
	pop_sequence(std::stack<std::string>& st, seq_db& db) : my_stack(st), rdb(db) { }
	
	template<typename It> void operator() ( It first, It last) const
	{
		std::string seq(my_stack.top());
		my_stack.pop();
		std::string name(my_stack.top());
		my_stack.pop();
		trim_string(name);
		sanitize_string(seq);
		
		rdb.add(name, seq);
	}
		
	std::stack<std::string>& my_stack;
	seq_db &rdb;
};


using namespace boost::spirit::classic;

struct seq_grammar : public grammar<seq_grammar>
{
	seq_grammar(std::stack<std::string>& st, seq_db& db ) : string_stack(st), rdb(db) {}
	
	template <typename ScannerT> struct definition
	{
		definition(seq_grammar const& self)
		{
			file_format =
					(*space_p) >>
					fasta_format
//				|	clustal
//				|	phylip
				;
			fasta_format = 
					*fasta_seq
				;
			fasta_seq =
				(	fasta_seq_head
				>>	fasta_seq_body
				)	[pop_sequence(self.string_stack, self.rdb)]
					;
			fasta_seq_head = 
					ch_p('>')
				>>	(+(graph_p|blank_p))[push_string(self.string_stack)]
				>>	eol_p
				;
			fasta_seq_body =
					(+(~ch_p('>')))[push_string(self.string_stack)]
				;
			
		}

		rule<ScannerT> file_format;
		rule<ScannerT> fasta_format;
		rule<ScannerT> fasta_seq;
		rule<ScannerT> fasta_seq_head;
		rule<ScannerT> fasta_seq_body;
		rule<ScannerT> const& start() const { return file_format; }
	};
	
	std::stack<std::string> &string_stack;
	seq_db &rdb;
};

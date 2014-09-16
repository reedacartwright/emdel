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

#ifndef SEQDB_H
#define SEQDB_H

#include <string>
#include <vector>
#include <map>

//bool parse_file(const char* csFile, StringVec &vNames, SeqVec &vSeqs);

class seq_db
{
public:
	typedef std::string sequence;
	typedef std::string name;
	typedef std::pair<name, sequence> value_type; 
	typedef std::vector<value_type>::size_type size_type;
	typedef std::map<name, size_type> name_map;
	
	inline bool add(const name& id, const sequence& s)
	{
		std::pair<name_map::iterator, bool> p = data_map.insert(make_pair(id, size()));
		if(p.second)
			data_vec.push_back(make_pair(id, s));
		else
			data_vec[p.first->second].second.append(s);
		return p.second;
	}
	inline void append(size_type idx, const sequence& s)
	{
		data_vec[idx].second.append(s);
	}
	
	inline size_type size() const { return data_vec.size(); }
	
	inline const value_type& operator[](size_type sz) const
	{
		return data_vec[sz];
	}
			
	inline void clear()
	{
		data_vec.clear();
		data_map.clear();
	}
	
	bool parse_file(const char *csfile, bool bappend=false);
	
protected:
	std::vector<value_type> data_vec;
	name_map data_map;
};

#endif

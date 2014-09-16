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
 
#ifndef EMDEL_SERIES_H
#define EMDEL_SERIES_H

#include <vector>

template<typename T>
class taylor_series
{
public:
	typedef T value_type;
	typedef std::vector<value_type> data_type;
	
	taylor_series(const value_type &p, const value_type* d, size_t sz) : point(p), data(d, d+sz) { }
	
	template<typename V, std::size_t N>
	taylor_series(const V &p, const V (&d)[N]) : point(p), data(d,d+N) { } 
	
	inline value_type operator()(value_type r) const
	{
		r -= point;
		if(data.empty())
			return value_type();
		value_type ret = data.back();
		for(typename data_type::const_reverse_iterator cit = data.rbegin()-1;
			cit != data.rend(); ++cit )
		{
			ret *= r;
			ret += *cit;
		}
		return ret;
	}

private:
	value_type point;
	data_type data;
};


#endif

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
 
#ifndef EMDEL_TABLE_H
#define EMDEL_TABLE_H

#include <vector>

template<typename T, typename A = std::allocator<T> >
class table : public std::vector<T, A>
{
	public:
		typedef std::vector<T, A> base_type;
		typedef typename base_type::size_type size_type;
		typedef typename base_type::value_type value_type;
		
		table(size_type a, size_type b, const value_type &v=value_type()) : base_type(a*b,v), sz_a(a), sz_b(b)  { }
		table() : base_type(), sz_a(0), sz_b(0) { }
		
		inline void resize(size_type a, size_type b, const value_type &v=value_type() )
		{
			sz_a = a;
			sz_b = b;
			base_type::resize(a*b,v);
		}
		
		inline void assign(size_type a, size_type b, const value_type &v=value_type() )
		{
			sz_a = a;
			sz_b = b;
			base_type::assign(a*b,v);
		}
		inline void assign(const value_type &v)
		{
			assign(size1(),size2(), v);
		}
		
		inline size_type size1() const { return sz_a; }
		inline size_type size2() const { return sz_b; }
		
		inline const value_type& operator()(size_type a, size_type b) const
		{
			return this->at((a)*sz_b+(b));
		}
		
		inline value_type& operator()(size_type a, size_type b)
		{
			return this->at((a)*sz_b+(b));
		}
		
	private:
		size_type sz_a;
		size_type sz_b;
};

template<typename T, typename A = std::allocator<T> >
		class torus : public std::vector<T, A>
{
	public:
		typedef std::vector<T, A> base_type;
		typedef typename base_type::size_type size_type;
		typedef typename base_type::value_type value_type;
		
		torus(size_type a, size_type b, const value_type &v=value_type()) : base_type(a*b,v), sz_a(a), sz_b(b)  { }
		torus() : base_type(), sz_a(0), sz_b(0) { }
		
		inline void resize(size_type a, size_type b, const value_type &v=value_type() )
		{
			sz_a = a;
			sz_b = b;
			base_type::resize(a*b,v);
		}
		
		inline void assign(size_type a, size_type b, const value_type &v=value_type() )
		{
			sz_a = a;
			sz_b = b;
			base_type::assign(a*b,v);
		}
		inline void assign(const value_type &v)
		{
			assign(size1(),size2(), v);
		}
		
		inline size_type size1() const { return sz_a; }
		inline size_type size2() const { return sz_b; }
		
		inline const value_type& operator()(size_type a, size_type b) const
		{
			return at((a%sz_a)*sz_b+(b%sz_b));
		}
		
		inline value_type& operator()(size_type a, size_type b)
		{
			return at((a%sz_a)*sz_b+(b%sz_b));
		}
		
	private:
		size_type sz_a;
		size_type sz_b;
};

#endif

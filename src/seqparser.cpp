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

#include "seqparser.h"

using namespace std;

bool parse_file(const char* cs, seq_db &rdb)
{
	file_iterator<char> file_first(cs);
	if(!file_first)
	{
		cerr << "ERROR: unable to open \'" << cs << "\'" << endl;
		return false;			
	}
	file_iterator<char>  file_last = file_first.make_end();
		
	stack<string> my_stack;
	seq_grammar my_grammar(my_stack, rdb);
	parse_info< file_iterator<char> > info = parse(file_first, file_last, my_grammar);
	if (!info.full)
	{
		cerr << "ERROR: unable to parse \'" << cs << "\'" << endl;
		return false;
	}
	return true;
}


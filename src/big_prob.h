#ifndef BIG_PROB_H
#define BIG_PROB_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <cfloat>
#include <iostream>

struct big_prob {
	static const int M = 512;
	double q;
	int n;
	
	big_prob(double a=0.0, int b=0) : q(a), n(b) { fix(); }
	
	operator double () const {
		return ldexp(q, M*n);
	}
	
	big_prob& operator = (double d) {
		q = d;
		n = 0;
		return fix();
	}
	
	inline big_prob& fix() {
		return *this;
		/*
		double f = fabs(q);
		if(7.458340731e-155 < f || f < 1.340780793e+154) {
			return *this;
		}
		
		int np;
		double qp = frexp(q,&np);
		if(np > M) {
			q = ldexp(qp, np-M);
			++n;
		} else {
			q = ldexp(qp, np+M);
			--n;
		}
		std::cerr << "x";
		return *this;
		*/
	}
	
	inline big_prob operator * (const big_prob& r) const {
		return big_prob(q*r.q, n+r.n);
	}
	inline big_prob operator / (const big_prob& r) const {
		return big_prob(q/r.q, n-r.n);
	}
	inline big_prob& operator *= (const big_prob& r) {
		q *= r.q;
		n += r.n;
		return fix();
	}
	inline big_prob& operator /= (const big_prob& r) {
		q /= r.q;
		n -= r.n;
		return fix();
	}
	
	inline big_prob operator + (const big_prob& r) const {
		big_prob ret;
		if(n == r.n) {
			ret.q = q+r.q;
			ret.n = n;
		} else if(n > r.n) {
			ret.q = q+ldexp(r.q, (r.n-n)*M);
			ret.n = n;
		} else {
			ret.q = r.q+ldexp(q, (n-r.n)*M);
			ret.n = r.n;
		}
		return ret.fix();
	}
	inline big_prob operator - (const big_prob& r) const {
		big_prob ret;
		if(n == r.n) {
			ret.q = q-r.q;
			ret.n = n;
		} else if(n > r.n) {
			ret.q = q-ldexp(r.q, (r.n-n)*M);
			ret.n = n;
		} else {
			ret.q = r.q-ldexp(q, (n-r.n)*M);
			ret.n = r.n;
		}
		return ret.fix();
	}
	inline big_prob& operator += (const big_prob& r) {
		if(n == r.n) {
			q += r.q;
		} else if(n > r.n) {
			q += ldexp(r.q, (r.n-n)*M);
		} else {
			q = r.q+ldexp(q, (n-r.n)*M);
			n = r.n;
		}
		return fix();
	}
	inline big_prob& operator -= (const big_prob& r) {
		if(n == r.n) {
			q -= r.q;
		} else if(n > r.n) {
			q -= ldexp(r.q, (r.n-n)*M);
		} else {
			q = r.q-ldexp(q, (n-r.n)*M);
			n = r.n;
		}
		return fix();
	}
	inline big_prob operator -() const {
		big_prob ret(*this);
		ret.q = -ret.q;
		return ret;
	}
	inline double log() const {
		return ::log(q)+(M*n)*M_LN2;
	}
	
};

template<class OS>
OS& operator<<(OS& os, const big_prob& r) {
	os << r.q << "x" << "2^" << big_prob::M*r.n;
	return os;
}

#endif

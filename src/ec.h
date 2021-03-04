/*
Copyright 2021, Yves Gallot

f12ecm is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include "res.h"

// Montgomery curve: B y^2 = x^3 + A x^2 + x
template<typename VComplex>
class EC
{
public:
	class Point
	{
	private:
		Res<VComplex> _x, _z;

	public:
		Point() {}	// init(thread_index) must be called
		Point(const size_t thread_index) : _x(thread_index), _z(thread_index) {}
		Point(const Point &) = delete;
		Point & operator=(const Point &) = delete;

		void init(const size_t thread_index) { _x.init(thread_index); _z.init(thread_index); }
		
		const Res<VComplex> & x() const { return _x; }
		const Res<VComplex> & z() const { return _z; }

		void set(const size_t i, const mpz_class & x, const mpz_class & z) { _x.set_z(i, x); _z.set_z(i, z); }
		void set(const Point & P) { _x.set(P._x); _z.set(P._z); }
		void swap(Point & P) { _x.swap(P._x); _z.swap(P._z); }

		void addsub() { Res<VComplex>::addsub(_x, _z); }
		void addsub(const Point & P) { _x.add(P._x, P._z); _z.sub(P._x, P._z); }

		void sqr() { _x.sqr(); _z.sqr(); }
		void to_multiplier() { _x.to_multiplier(); _z.to_multiplier(); }
		void cross(const Point & P) { _x.mul_m(P._z); _z.mul_m(P._x); }
		void mul_mm(const Point & P) { _x.mul_mm(P._x); _z.mul_mm(P._z); }

		void vmul(Point & P) { _z.to_multiplier(); _x.mul_m(_z); _z.set(P._x); P._z.to_multiplier(); _z.mul_m(P._z); }
		void mulc_xz(const Point & P, const Res<VComplex> & C, Res<VComplex> & t) { t.sub(P._x, P._z); _x.set(t); t.mul_m(C); _z.add(P._z, t); }
	};

private:
	Res<VComplex> _C;	// C = (A + 2) / 4
	Point _P1, _P2, _T;	// buffers
	Res<VComplex> _t;	// buffer

	// P + P
	void dbl(Point & P)
	{
		P.addsub();
		P.sqr();
		_T.mulc_xz(P, _C, _t);
		P.vmul(_T);
	}

	// P1 + P2, where Pm = P1 - P2
	void sum(Point & P1, const Point & P2, const Point & Pm)
	{
		P1.addsub();
		_T.addsub(P2);
		_T.to_multiplier();
		P1.cross(_T);
		P1.addsub();
		P1.sqr();
		_T.set(Pm);
		_T.to_multiplier();
		P1.cross(_T);
	}

	void sum_dbl(Point & P1, Point & P2, const Point & Pm)	// 28 transforms ~ 14 squares
	{
		P1.addsub();
		P2.addsub();
		_T.set(P2);
		_T.to_multiplier();		// 2 transforms
		P1.cross(_T);			// 4 transforms
		P1.addsub();
		P1.sqr();				// 4 transforms
		_T.set(Pm);
		_T.to_multiplier();		// 2 transforms
		P1.cross(_T);			// 4 transforms
		P2.sqr();				// 4 transforms
		_T.mulc_xz(P2, _C, _t);	// 2 transforms
		P2.vmul(_T);			// 6 transforms
	}

public:
	EC(const size_t thread_index) : _C(thread_index), _P1(thread_index), _P2(thread_index), _T(thread_index), _t(thread_index) {}

	void set(const size_t i, const mpz_class & A, const mpz_class & n)
	{
		mpz_class t; mpz_invert(t.get_mpz_t(), mpz_class(4).get_mpz_t(), n.get_mpz_t());
		_C.set_z(i, ((A + 2) * t) % n);
		_C.to_multiplier();
	}

	void dbl(Point & Po, const Point & P) { Po.set(P); dbl(Po); }
	void sum(Point & Po, const Point & P1, const Point & P2, const Point & Pm) { _P1.set(P1); sum(_P1, P2, Pm); _P1.swap(Po); }

	// m * P
	void mul(Point & Po, const Point & P, const uint64_t m)	// m >= 3
	{
		_P1.set(P); dbl(_P2, P);
		const int log2_m = 63 - __builtin_clzll(m);
		for (uint64_t b = uint64_t(1) << (log2_m - 1); b > 1; b >>= 1)
		{
			const bool isSet = ((m & b) != 0);
			Point & Q1 = isSet ? _P1 : _P2;
			Point & Q2 = isSet ? _P2 : _P1;
			sum_dbl(Q1, Q2, P);
		}
		if (m % 2 != 0) sum(_P1, _P2, P); else dbl(_P1);
		Po.swap(_P1);
	}
};

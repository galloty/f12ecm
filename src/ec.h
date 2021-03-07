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
		void addsub(const Point & P) { if (this == &P) addsub(); else { _x.add(P._x, P._z); _z.sub(P._x, P._z); } }

		void sqr() { _x.sqr(); _z.sqr(); }
		void to_multiplier() { _x.to_multiplier(); _z.to_multiplier(); }
		void cross(const Point & P) { _x.mul_m(P._z); _z.mul_m(P._x); }
		void mul_mm(const Point & P) { _x.mul_mm(P._x); _z.mul_mm(P._z); }

		void vmul(Point & P) { _z.to_multiplier(); _x.mul_m(_z); _z.set(P._x); P._z.to_multiplier(); _z.mul_m(P._z); }
		void mulc_xz(const Point & P, const Res<VComplex> & C, Res<VComplex> & t) { t.sub(P._x, P._z); _x.set(t); t.mul_m(C); _z.add(P._z, t); }
	};

private:
	Res<VComplex> _A2_4, _t;		// _A2_4 = (A + 2) / 4
	Point _A, _B, _C, _T1, _T2, _T, _Tm;

public:
	// Pr = P + P
	void dbl(Point & Pr, const Point & P)	// 2 mul + 2 sqr + 1 mul_const
	{
		Pr.addsub(P);				// P'.x = P.x + P.z, P'.z = P.x - P.z
		Pr.sqr();					// P'.x = (P.x + P.z)^2, P'.z = (P.x - P.z)^2
		_T.mulc_xz(Pr, _A2_4, _t);	// T.x = P'.x - P'.z, T.z = C * T.x + P'.z
		Pr.vmul(_T);				// P'.x =  P'.x * P'.z, P'.z = T.x * T.z
	}

	// Pr = P1 + P2, where P1 != P2 and Pm = P1 - P2 or Pm = P2 - P1
	void add(Point & Pr, const Point & P1, const Point & P2, const Point & Pm)	// 4 mul + 2 sqr
	{
		// Pr.set(Pm); return;
		_T.addsub(P2);				// T.x = P2.x + P2.z, T.z = P2.x - P2.z
		_T.to_multiplier();
		_Tm.set(Pm);
		_Tm.to_multiplier();
		Pr.addsub(P1);				// P'.x = P1.x + P1.z, P'.z = P1.x - P1.z
		Pr.cross(_T);				// P'.x = P'.x * T.z, P'.z = P'.z * T.x
		Pr.addsub();				// P".x = P'.x + P'.z, P".z = P'.x - P'.z
		Pr.sqr();					// P".x = (P'.x + P'.z)^2, P".z = (P'.x - P'.z)^2
		Pr.cross(_Tm);				// P".x = P".x * Pm.z, P1".z = P".z * Pm.x
	}

public:
	EC(const size_t thread_index) : _A2_4(thread_index), _t(thread_index),
		_A(thread_index), _B(thread_index), _C(thread_index), _T1(thread_index), _T2(thread_index), _T(thread_index), _Tm(thread_index) {}

	void set(const size_t i, const mpz_class & A, const mpz_class & n)
	{
		mpz_class t; mpz_invert(t.get_mpz_t(), mpz_class(4).get_mpz_t(), n.get_mpz_t());
		_A2_4.set_z(i, ((A + 2) * t) % n);
		_A2_4.to_multiplier();
	}

	// Pr = n * P
	// PRAC algorithm: Montgomery, Peter L. (1983). "Evaluating recurrences of form X_{m+n} = f(X_m, X_n, X_{m-n}) via Lucas Chains", unpublished.
	void mul(Point & Pr, const Point & P, const uint64_t n, const uint64_t p = 0)	// n >= 3, p | n
	{
		_A.set(P);
		uint64_t d = n;

		while (d != 1)
		{
			const uint64_t r = (p == 0) ? std::llrint(d * ((std::sqrt(5) - 1) / 2)) : (d / p) * std::llrint(p * ((std::sqrt(5) - 1) / 2));
			uint64_t e = 2 * r - d; d -= r;
			_B.set(_A); _C.set(_A); dbl(_A, _A);

			while (d != e)
			{
				if (d < e)
				{
					_A.swap(_B); // C = -C;
					std::swap(d, e);
				}
				// e <= d

				// #1: A' = 2A + B, B' = 2B + A
				if ((4 * d <= 5 * e) && ((d + e) % 3 == 0))
				{
					const uint64_t t = (2 * e - d) / 3; d = (2 * d - e) / 3; e = t;
					add(_T1, _A, _B, _C); add(_T2, _T1, _A, _B);
					add(_B, _T1, _B, _A); _A.swap(_T2);
				}
				// #2
				else if ((4 * d <= 5 * e) && ((d - e) % 6 == 0))
				{
					d = (d - e) / 2;
					add(_B, _A, _B, _C); dbl(_A, _A);
				}
				// #3
				else if (d < 4 * e)
				{
					d -= e;
					add(_C, _A, _B, _C); _B.swap(_C); // C = -C;
				}
				// #4
				else if ((d - e) % 2 == 0)
				{
					d = (d - e) / 2;
					add(_B, _A, _B, _C); dbl(_A, _A);
				}
				// #5
				else if (d % 2 == 0)
				{
					d /= 2;
					add(_C, _A, _C, _B); dbl(_A, _A);
				}
				// #6
				else if (d % 3 == 0)
				{
					d = d / 3 - e;
					dbl(_T1, _A); add(_T2, _A, _B, _C);
					add(_A, _T1, _A, _A); add(_C, _T1, _T2, _C); _B.swap(_C); // C = -C;
				}
				// #7
				else if ((d + e) % 3 == 0)
				{
					d = (d - 2 * e) / 3;
					add(_T1, _A, _B, _C); dbl(_T2, _A);
					add(_B, _T1, _A, _B); add(_A, _T2, _A, _A);
				}
				// #8
				else if ((d - e) % 3 == 0)
				{
					d = (d - e) / 3;
					add(_T1, _A, _B, _C); add(_C, _A, _C, _B);
					_B.swap(_T1);
					dbl(_T1, _A); add(_A, _T1, _A, _A);
				}
				// #9
				else if (e % 2 == 0)
				{
					e = e / 2;
					add(_C, _C, _B, _A); dbl(_B, _B);	// C' = C + (-B)
				}
				else
				{
					throw std::runtime_error("PRAC failed");
				}
			}

			add(_C, _A, _B, _C); _A.swap(_C);
		}

		Pr.swap(_A);
	}
};

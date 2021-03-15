/*
Copyright 2021, Yves Gallot

f12ecm is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include "res.h"
#include "prac.h"

// Montgomery curve: B y^2 = x^3 + A x^2 + x
template<typename VComplex>
class EC_m
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
	static mpz_class cub_mod(const mpz_class & x, const mpz_class & n) { mpz_class t = (x * x * x) % n; if (t < 0) t += n; return t; }

public:
	// A such that the group order is divisible by 12 (H. Suyama)
	static mpz_class A_12(const size_t i, const uint64_t sigma, const mpz_class & n, Point & P0)
	{
		const mpz_class s = mpz(sigma);
		const mpz_class u = s * s - 5, v = 4 * s;
		const mpz_class x = cub_mod(u, n), z = cub_mod(v, n);
		P0.set(i, x, z);
		mpz_class t = 4 * x * v;
		mpz_invert(t.get_mpz_t(), t.get_mpz_t(), n.get_mpz_t());
		t = (t * cub_mod(v - u, n) * (3 * u + v)) % n;
		t -= 2; if (t < 0) t += n;
		return t;
	}

private:
	Res<VComplex> _A2_4, _t;		// _A2_4 = (A + 2) / 4
	Point _A, _B, _C, _T1, _T2, _T, _Tm;

public:
	// Pr = P + P
	void dbl(Point & Pr, const Point & P)	// 2 mul + 2 sqr + 1 mul_const: 12 transforms
	{
		Pr.addsub(P);				// P'.x = P.x + P.z, P'.z = P.x - P.z
		Pr.sqr();					// P'.x = (P.x + P.z)^2, P'.z = (P.x - P.z)^2
		_T.mulc_xz(Pr, _A2_4, _t);	// T.x = P'.x - P'.z, T.z = C * T.x + P'.z
		Pr.vmul(_T);				// P'.x =  P'.x * P'.z, P'.z = T.x * T.z
	}

	// Pr = P1 + P2, where P1 != P2 and Pm = P1 - P2 or Pm = P2 - P1
	void add(Point & Pr, const Point & P1, const Point & P2, const Point & Pm)	// 4 mul + 2 sqr: 16 transforms
	{
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
	EC_m(const size_t thread_index) : _A2_4(thread_index), _t(thread_index),
		_A(thread_index), _B(thread_index), _C(thread_index), _T1(thread_index), _T2(thread_index), _T(thread_index), _Tm(thread_index) {}

	void set(const size_t i, const mpz_class & A, const mpz_class & n)
	{
		mpz_class t; mpz_invert(t.get_mpz_t(), mpz_class(4).get_mpz_t(), n.get_mpz_t());
		_A2_4.set_z(i, ((A + 2) * t) % n);
	}

	void init()
	{
		_A2_4.to_multiplier();
	}

	// Pr = n * P
	// PRAC algorithm: Montgomery, Peter L. (1983). "Evaluating recurrences of form X_{m+n} = f(X_m, X_n, X_{m-n}) via Lucas Chains", unpublished.
	void mul(Point & Pr, const Point & P, const uint64_t n, const uint64_t p = 0)	// n >= 3, p | n
	{
		const double phi = 1.6180339887498948482045868343656381177;
		const double alpha[10] = { 1 / phi, 5 / (phi + 7), 1429 / (phi + 2311), 3739 / (6051 - phi), 79 / (129 - phi),
								   31 / (phi + 49), 209 / (phi + 337), 11 / (19 - phi), 545 / (883 - phi), 1 / (3 - phi) };

		size_t c_min = size_t(-1), best_i = 0;
		for (size_t i = 0; i < 10; ++i)
		{
			size_t sqr_count, mul_count;
			const size_t c = PRAC::cost(n, p, alpha[i], sqr_count, mul_count);
			if (c < c_min)
			{
				c_min = c; best_i = i;
			}
		}

		_A.set(P);
		uint64_t d = n;

		while (d != 1)
		{
			const uint64_t r = (p == 0) ? std::llrint(d * alpha[best_i]) : (d / p) * std::llrint(p * alpha[best_i]);
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

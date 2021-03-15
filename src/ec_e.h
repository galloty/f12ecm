/*
Copyright 2021, Yves Gallot

f12ecm is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include "res.h"

// Edwards curve: x^2 + y^2 = 1 + d x^2 y^2
template<typename VComplex>
class EC_e
{
private:
	// Montgomery construction: Daniel J. Bernstein; Tanja Lange; Peter Birkner; Christiane Peters, ECM using Edwards curves, § 7.7
	class EC_mc
	{
	public:
		struct Point
		{
			mpz_class x, y;
		};

	private:
		static Point add(const Point & p1, const Point & p2, const mpz_class & m)
		{
			mpz_class lambda, den;
			if (p1.x != p2.x)
			{
				lambda = p1.y - p2.y;
				den = p1.x - p2.x;
			}
			else
			{
				lambda = 3 * p1.x * p1.x - 12;
				den = p1.y + p2.y;
			}
			mpz_invert(den.get_mpz_t(), den.get_mpz_t(), m.get_mpz_t());
			lambda = (lambda * den) % m;

			Point p3;
			p3.x = (lambda * lambda - p1.x - p2.x) % m;
			p3.y = (m - (lambda * (p3.x - p1.x) + p1.y)) % m;
			if (p3.x < 0) p3.x += m;
			if (p3.y < 0) p3.y += m;
			return p3;
		}

	public:
		// (s, t) = n * [-2, 4] on y^2 = x^3 − 12*x (modulo m)
		static Point get(const uint64_t n, const mpz_class & m)
		{
			Point p0; p0.x = m - 2; p0.y = 4;
			Point p = p0;

			for (int b = 62 - __builtin_clzll(n); b >= 0; --b)
			{
				p = add(p, p, m);
				if ((n & (uint64_t(1) << b)) != 0) p = add(p, p0, m);
			}

			return p;
		}
	};

public:
	class Point
	{
	private:
		Res<VComplex> _x, _y, _z;

	public:
		Point() {}	// init(thread_index) must be called
		Point(const size_t thread_index) : _x(thread_index), _y(thread_index), _z(thread_index) {}
		Point(const Point &) = delete;
		Point & operator=(const Point &) = delete;

		void init(const size_t thread_index) { _x.init(thread_index); _y.init(thread_index); _z.init(thread_index); }
		
		const Res<VComplex> & x() const { return _x; }
		const Res<VComplex> & y() const { return _y; }
		const Res<VComplex> & z() const { return _z; }

		Res<VComplex> & x() { return _x; }
		Res<VComplex> & y() { return _y; }
		Res<VComplex> & z() { return _z; }

		void set(const size_t i, const mpz_class & x, const mpz_class & y) { _x.set_z(i, x); _y.set_z(i, y); _z.set_z(i, mpz_class(1));}
		void set(const Point & P) { _x.set(P._x); _y.set(P._y); _z.set(P._z); }
		void swap(Point & P) { _x.swap(P._x); _y.swap(P._y); _z.swap(P._z); }
	};

private:
	static mpz_class sqr_mod(const mpz_class & x, const mpz_class & n) { return (x * x) % n; }
	static mpz_class cub_mod(const mpz_class & x, const mpz_class & n) { mpz_class t = (x * x * x) % n; if (t < 0) t += n; return t; }

public:
	// Montgomery construction of Edwards curve: the group order is divisible by 12
	static mpz_class d_12(const size_t i, const uint64_t sigma, const mpz_class & n, Point & P0)
	{
		// (s, t) = sigma * [-2, 4] on y^2 = x^3 − 12*x (modulo n)
		const typename EC_mc::Point P = EC_mc::get(sigma, n);
		const mpz_class s = P.x, t = P.y;

		const mpz_class s2 = sqr_mod(s, n), t2 = sqr_mod(t, n);
		const mpz_class a = s2 - 12 * s - 12, b = ((s - 2) * (s + 6)) % n;
		mpz_class den;

		den = (1024 * s2 * t2) % n; mpz_invert(den.get_mpz_t(), den.get_mpz_t(), n.get_mpz_t());
		mpz_class d = (-a * cub_mod(b, n) * den) % n;
		if (d < 0) d += n;

		den = (b * (s2 + 12 * s - 12)) % n; mpz_invert(den.get_mpz_t(), den.get_mpz_t(), n.get_mpz_t());
		mpz_class x = (8 * t * (s2 + 12) * den) % n;
		if (x < 0) x += n;

		den = (b * (s2 - 12)) % n; mpz_invert(den.get_mpz_t(), den.get_mpz_t(), n.get_mpz_t());
		mpz_class y = (-4 * s * a * den) % n;
		if (y < 0) y += n;

		P0.set(i, x, y);
		return d;
	}

private:
	Res<VComplex> _d, _t;
	Res<VComplex> _A, _B, _C, _D, _E;
	Point _P0;

public:
	// Pr = P + P
	void dbl(Point & Pr, const Point & P)	// 4 T + 4 MM + 1 S + 2 M: 14 transforms
	{
		_A.set(P.x()); _A.to_multiplier();
		_B.set(P.y()); _B.to_multiplier();

		_C.set(_A);
		_C.mul_mm(_B);							// C = xy
		_C.add(_C, _C);							// C = 2xy
		_A.mul_mm(_A);							// A = x^2, C = 2xy
		_B.mul_mm(_B);							// A = x^2, B = y^2, C = 2xy

		Res<VComplex>::addsub(_A, _B);			// A = x^2 + y^2, B = x^2 - y^2, C = 2xy

		_D.set(P.z());
		_D.sqr();								// A = x^2 + y^2, B = x^2 - y^2, C = 2xy, D = z^2

		_D.add(_D, _D);
		_D.sub(_A, _D);							// A = x^2 + y^2, B = x^2 - y^2, C = 2xy, D = x^2 + y^2 - 2z^2

		_D.to_multiplier();
		_C.mul_m(_D);
		Pr.x().swap(_C);

		_A.to_multiplier();
		_B.mul_m(_A);
		Pr.y().swap(_B);

		_A.mul_mm(_D);
		Pr.z().swap(_A);
	}

	// Pr = P1 + P2
	void add(Point & Pr, const Point & P1, const Point & P2)	// 11 T + 8 + 5*2: 29 transforms
	{
		_A.set(P1.x()); _A.to_multiplier();
		_B.set(P2.x()); _B.to_multiplier();
		_C.set(P1.y()); _C.to_multiplier();
		_D.set(P2.y()); _D.to_multiplier();

		_E.set(_A);
		_E.mul_mm(_B);							// E = x12

		_A.mul_mm(_D);							// A = x1y2, E = x12
		_D.mul_mm(_C);							// A = x1y2, D = y12, E = x12
		_C.mul_mm(_B);							// A = x1y2, C = x2y1, D = y12, E = x12
		_D.sub(_D, _E);							// A = x1y2, C = x2y1, D = y12 - x12

		_B.set(_A);
		_t.set(_C); _t.to_multiplier();
		_B.mul_m(_t);
		_B.mul_m(_d);							// A = x1y2, B = d * x1x2y1y2, C = x2y1, D = y12 - x12
		_C.add(_C, _A);							// B = d * x1x2y1y2, C = x1y2 + x2y1, D = y12 - x12

		_A.set(P1.z());
		_t.set(P2.z()); _t.to_multiplier();
		_A.mul_m(_t);							// A = z12, B = d * x1x2y1y2, C = x1y2 + x2y1, D = y12 - x12

		_A.to_multiplier();
		_C.mul_m(_A);
		_D.mul_m(_A);
		_A.mul_mm(_A);							// A = z12^2, B = d * x1x2y1y2, C = z12 * (x1y2 + x2y1), D = z12 * (y12 - x12)

		Res<VComplex>::addsub(_A, _B);			// A = z12^2 + d * x1x2y1y2, B = z12^2 - d * x1x2y1y2, C = z12 * (x1y2 + x2y1), D = z12 * (y12 - x12)

		_B.to_multiplier();
		_C.mul_m(_B);
		Pr.x().swap(_C);

		_A.to_multiplier();
		_D.mul_m(_A);
		Pr.y().swap(_D);

		_A.mul_mm(_B);
		Pr.z().swap(_A);
	}

public:
	EC_e(const size_t thread_index) : _d(thread_index), _t(thread_index),
		_A(thread_index), _B(thread_index), _C(thread_index), _D(thread_index), _E(thread_index),
		_P0(thread_index) {}

	void set(const size_t i, const mpz_class & d)
	{
		_d.set_z(i, d);
	}

	void init()
	{
		_d.to_multiplier();
	}

	// Pr = n * P
	void mul(Point & Pr, const Point & P, const uint64_t n)
	{
		_P0.set(P);
		if (&Pr != &P) Pr.set(P);

		for (int b = 62 - __builtin_clzll(n); b >= 0; --b)
		{
			dbl(Pr, Pr);
			if ((n & (uint64_t(1) << b)) != 0) add(Pr, Pr, _P0);
		}
	}
};

/*
Copyright 2021, Yves Gallot

f12ecm is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <vector>

#include "res.h"

// Edwards curve: x^2 + y^2 = 1 + d x^2 y^2
template<typename VComplex>
class EC_e
{
private:
	template<int a>
	class EC_modular
	{
	private:
		struct Point
		{
			mpz_class x, y;
		};

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
				lambda = 3 * p1.x * p1.x + a;
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
		// (s, t) = n * [s1, t1] on T^2 = S^3 + a*S + b (modulo m) [b is not needed]
		static void get(const uint64_t n, const int s1, const int t1, const mpz_class & m,  mpz_class & s,  mpz_class & t)
		{
			Point p0;
			p0.x = s1; if (s1 < 0) p0.x += m;
			p0.y = t1; if (t1 < 0) p0.y += m;
			Point p = p0;

			for (int b = 62 - __builtin_clzll(n); b >= 0; --b)
			{
				p = add(p, p, m);
				if ((n & (uint64_t(1) << b)) != 0) p = add(p, p0, m);
			}

			s = p.x; t = p.y;
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

		void set(const size_t i, const mpz_class & x, const mpz_class & y) { _x.set_z(i, x); _y.set_z(i, y); _z.set1(i);}
		void set(const Point & P) { _x.set(P._x); _y.set(P._y); _z.set(P._z); }
		void swap(Point & P) { _x.swap(P._x); _y.swap(P._y); _z.swap(P._z); }
	};

private:
	static mpz_class sqr_mod(const mpz_class & x, const mpz_class & n) { return (x * x) % n; }
	static mpz_class cub_mod(const mpz_class & x, const mpz_class & n) { mpz_class t = (x * x * x) % n; if (t < 0) t += n; return t; }

public:
	// Daniel J. Bernstein; Tanja Lange; Peter Birkner; Christiane Peters, ECM using Edwards curves, § 7.7

	// Montgomery construction: the group order of Edwards curve is divisible by 12
	static mpz_class d_12m(const size_t i, const uint64_t m, const mpz_class & n, Point & P0)
	{
		// (s, t) = m * [-2, 4] on t^2 = s^3 − 12*s (modulo n)
		mpz_class s, t;
		EC_modular<-12>::get(m, -2, 4, n, s, t);

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

	// Atkin-Morain construction: the group order of Edwards curve is divisible by 16
	static mpz_class d_16am(const size_t i, const uint64_t m, const mpz_class & n, Point & P0)
	{
		// (s, t) = m * [12, 40] on t^2 = s^3 − 8*s - 32 (modulo n)
		mpz_class s, t;
		EC_modular<-8>::get(m, 12, 40, n, s, t);

		mpz_class den;

		den = s - 9; if (den < 0) den += n; mpz_invert(den.get_mpz_t(), den.get_mpz_t(), n.get_mpz_t());
		den = ((t + 25) * den + 1) % n;
		mpz_class alpha; mpz_invert(alpha.get_mpz_t(), den.get_mpz_t(), n.get_mpz_t());

		den = 8 * sqr_mod(alpha, n) - 1; if (den < 0) den += n; mpz_invert(den.get_mpz_t(), den.get_mpz_t(), n.get_mpz_t());
		const mpz_class beta = (2 * alpha * (4 * alpha + 1) * den) % n;

		mpz_class t1 = 2 * beta - 1; if (t1 < 0) t1 += n;
		const mpz_class t2 = (t1 * t1) % n;

		mpz_class d = 2 * t2 - 1; if (d < 0) d += n;
		den = (t2 * t2) % n; mpz_invert(den.get_mpz_t(), den.get_mpz_t(), n.get_mpz_t());
		d = (d * den) % n;

		den = 6 * beta - 5; if (den < 0) den += n; mpz_invert(den.get_mpz_t(), den.get_mpz_t(), n.get_mpz_t());
		mpz_class x = (t1 * (4 * beta - 3) * den) % n;
		if (x < 0) x += n;

		den = ((t + 3 * s - 2) * (t + s + 16)) % n; mpz_invert(den.get_mpz_t(), den.get_mpz_t(), n.get_mpz_t());
		mpz_class y = (t1 * (sqr_mod(t, n) + 50 * t - 2 * cub_mod(s, n) + 27 * sqr_mod(s, n) - 104)) % n;
		y = (y * den) % n;
		if (y < 0) y += n;

		P0.set(i, x, y);
		return d;
	}

	// Gallot construction: the group order of Edwards curve is divisible by 16
	static mpz_class d_16g(const size_t i, const uint64_t m, const mpz_class & n, Point & P0)
	{
		// (s, t) = m * [4, 8] on t^2 = s^3 + 4*s - 16 (modulo n)
		mpz_class s, t;
		EC_modular<4>::get(m, 4, 8, n, s, t);

		mpz_class den;

		den = s - 4; if (den < 0) den += n; mpz_invert(den.get_mpz_t(), den.get_mpz_t(), n.get_mpz_t());
		const mpz_class alpha = ((t + 8) * den) % n;
		const mpz_class alpha2 = sqr_mod(alpha, n); 

		den = 8 - alpha2; if (den < 0) den += n; mpz_invert(den.get_mpz_t(), den.get_mpz_t(), n.get_mpz_t());
		const mpz_class r = ((8 + 2 * alpha) * den) % n;
		const mpz_class t1 = sqr_mod(2 * r - 1, n);

		mpz_class d = (8 * r * (r - 1) + 1) % n; if (d < 0) d += n;
		den = sqr_mod(t1, n); mpz_invert(den.get_mpz_t(), den.get_mpz_t(), n.get_mpz_t());
		d = (d * den) % n;

		mpz_class x = ((8 - alpha2) * (2 * r * r - 1)) % n;
		den = 2 * s - alpha2 + 4; mpz_invert(den.get_mpz_t(), den.get_mpz_t(), n.get_mpz_t());
		x = (x * den) % n; if (x < 0) x += n;

		den = 4 * r - 3; if (den < 0) den += n; mpz_invert(den.get_mpz_t(), den.get_mpz_t(), n.get_mpz_t());
		const mpz_class y = (t1 * den) % n;

		P0.set(i, x, y);
		return d;
	}

public:
	static const size_t W = size_t(1) << 8;
	typedef short naf_type;

private:
	Res<VComplex> _d, _D, _E, _F;
	Point _Pi[W / 4];
	size_t _dbl_count, _add_count;

public:
	// P = P + P
	void dbl(Point & P)	// 5 T + 5 MM + 1 S(2) + 1 M(2): 14 transforms
	{
		++_dbl_count;

		Res<VComplex> & A = P.x();
		Res<VComplex> & B = P.y();
		Res<VComplex> & C = P.z();

		A.to_m();
		B.to_m();

		_D.set(A);
		_D.mul_mm_add(B);					// A = Tx, B = Ty, C = z, D = 2xy
		A.mul_mm(A);
		B.mul_mm(B);						// A = x^2, B = y^2, C = z, D = 2xy

		Res<VComplex>::addsub(A, B);		// A = x^2 + y^2, B = x^2 - y^2, C = z, D = 2xy

		C.sqr_add();
		C.sub(A, C);						// A = x^2 + y^2, B = x^2 - y^2, C = x^2 + y^2 - 2z^2, D = 2xy

		B.mul(A);							// A = T(x^2 + y^2), B = y', C = x^2 + y^2 - 2z^2, D = 2xy
		A.mul_t(C);							// A = z', B = y', C = T(x^2 + y^2 - 2z^2), D = 2xy
		C.mul_t(_D);						// A = z', B = y', C = x'
		A.swap(C);
	}

	// if m = 1 then P1 = P1 + P2, if m = -1 then P1 = P1 - P2
	void add(Point & P1, const Point & P2, const double m = 1)	// 9 T + 6 MM + 7 M(2): 29 transforms
	{
		++_add_count;

		Res<VComplex> & A = P1.x();
		Res<VComplex> & B = P1.y();
		Res<VComplex> & C = P1.z();

		_D.set(P2.x());
		_E.set(P2.y());						// A = x1, B = y1, C = z1, D = x2, E = y2

		A.to_m();
		_F.set(A);
		_F.mul_t(_D);						// A = Tx1, B = y1, C = z1, D = Tx2, E = y2, F = x1x2
		A.mul_t(_E);						// A = x1y2, B = y1, C = z1, D = Tx2, E = Ty2, F = x1x2
		_E.mul_t(B);						// A = x1y2, B = Ty1, C = z1, D = Tx2, E = y1y2, F = x1x2
		B.mul_mm(_D);						// A = x1y2, B = x2y1, C = z1, E = y1y2, F = x1x2
		_E.sub(_E, _F, m);					// A = x1y2, B = x2y1, C = z1, E = y1y2 - x1x2

		_D.set(P2.z());
		C.mul(_D);							// A = x1y2, B = x2y1, C = z1z2, E = y1y2 - x1x2

		_D.set(A);
		A.add(A, B, m);						// A = x1y2 + x2y1, C = z1z2, D = x1y2, E = y1y2 - x1x2
		_D.mul(B);
		_D.mul_m(_d);						// A = x1y2, C = z1z2, D = d * x1x2y1y2, E = y1y2 - x1x2

		A.mul(C);
		_E.mul_m(C);
		C.mul_mm(C);						// A = z1z2 * (x1y2 + x2y1), C = z1z2^2, D = d * x1x2y1y2, E = z1z2 * (y1y2 - x1x2)

		Res<VComplex>::addsub(C, _D, m);	// A = z1z2 * (x1y2 + x2y1), C = z1z2^2 + d * x1x2y1y2, D = z1z2^2 - d * x1x2y1y2, E = z1z2 * (y1y2 - x1x2)

		A.mul(_D);							// A = x', C = z1z2^2 + d * x1x2y1y2, D = T(z1z2^2 - d * x1x2y1y2), E = z1z2 * (y1y2 - x1x2)
		_E.mul(C);							// A = x', C = T(z1z2^2 + d * x1x2y1y2), D = T(z1z2^2 - d * x1x2y1y2), E = y'
		C.mul_mm(_D);						// A = x', C = z', E = y'
		B.swap(_E);
	}

public:
	EC_e(const size_t thread_index) : _d(thread_index), _D(thread_index), _E(thread_index), _F(thread_index)
	{
		for (size_t i = 0; i < W / 4; ++i) _Pi[i].init(thread_index);
	}

	void set(const size_t i, const mpz_class & d)
	{
		_d.set_z(i, d);
	}

	void init()
	{
		_d.to_m();
		_dbl_count = 0; _add_count = 0;
	}

	void getCounters(size_t & dbl_count, size_t & add_count, size_t & cost) const
	{
		dbl_count = _dbl_count;
		add_count = _add_count;
		cost = 14 * _dbl_count + 29 * add_count;
	}

	// P = n * P
	void mul(Point & P, const mpz_class & n)
	{
		mpz_class ec = n;
		mpz_ptr e = ec.get_mpz_t();

		std::vector<naf_type> naf_vec; naf_vec.reserve(mpz_sizeinbase(e, 2));

		for (; mpz_size(e) != 0; mpz_fdiv_q_2exp(e, e, 1))
		{
			int ei = 0;
			if (mpz_odd_p(e))
			{
				ei = mpz_getlimbn(e, 0) & (W - 1);
				if (ei >= int(W/2)) ei -= W;
				if (ei > 0) mpz_sub_ui(e, e, ei); else mpz_add_ui(e, e, -ei);
			}
			naf_vec.push_back(naf_type(ei));
		}

		_Pi[0].set(P);
		_Pi[W / 4 - 1].set(P); dbl(_Pi[W / 4 - 1]);
		for (size_t i = 1; i < W / 4 - 1; ++i)
		{
			_Pi[i].set(_Pi[i - 1]);
			add(_Pi[i], _Pi[W / 4 - 1]);
		}
		add(_Pi[W / 4 - 1], _Pi[W / 4 - 2]);
	
		const size_t len = naf_vec.size() - 1;
		const naf_type * const naf = naf_vec.data();

		P.set(_Pi[(naf[len] - 1) / 2]);

		for (size_t i = 0; i < len; ++i)
		{
			dbl(P);
			const int ni = naf[len - 1 - i];
			if (ni != 0) add(P, _Pi[(std::abs(ni) - 1) / 2], (ni > 0) ? 1.0 : -1.0);
		}
	}
};

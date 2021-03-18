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

		void set(const size_t i, const mpz_class & x, const mpz_class & y) { _x.set_z(i, x); _y.set_z(i, y); _z.set1(i);}
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

		A.to_multiplier();
		B.to_multiplier();

		_D.set(A);
		_D.mul_mm(B, 2);					// A = Tx, B = Ty, C = z, D = 2xy
		A.mul_mm(A);
		B.mul_mm(B);						// A = x^2, B = y^2, C = z, D = 2xy

		Res<VComplex>::addsub(A, B);		// A = x^2 + y^2, B = x^2 - y^2, C = z, D = 2xy

		C.sqr(2);
		C.sub(A, C);						// A = x^2 + y^2, B = x^2 - y^2, C = x^2 + y^2 - 2z^2, D = 2xy

		B.mul(A);							// A = T(x^2 + y^2), B = y', C = x^2 + y^2 - 2z^2, D = 2xy
		A.mul_t(C);						// A = z', B = y', C = T(x^2 + y^2 - 2z^2), D = 2xy
		C.mul_t(_D);						// A = z', B = y', C = x'
		A.swap(C);
	}

	// P1 = P1 + P2
	void add(Point & P1, const Point & P2)	// 9 T + 6 MM + 7 M(2): 29 transforms
	{
		++_add_count;

		Res<VComplex> & A = P1.x();
		Res<VComplex> & B = P1.y();
		Res<VComplex> & C = P1.z();

		_D.set(P2.x());
		_E.set(P2.y());						// A = x1, B = y1, C = z1, D = x2, E = y2

		A.to_multiplier();
		_F.set(A);
		_F.mul_t(_D);						// A = Tx1, B = y1, C = z1, D = Tx2, E = y2, F = x1x2
		A.mul_t(_E);						// A = x1y2, B = y1, C = z1, D = Tx2, E = Ty2, F = x1x2
		_E.mul_t(B);						// A = x1y2, B = Ty1, C = z1, D = Tx2, E = y1y2, F = x1x2
		B.mul_mm(_D);						// A = x1y2, B = x2y1, C = z1, E = y1y2, F = x1x2
		_E.sub(_E, _F);						// A = x1y2, B = x2y1, C = z1, E = y1y2 - x1x2

		_D.set(P2.z());
		C.mul(_D);							// A = x1y2, B = x2y1, C = z1z2, E = y1y2 - x1x2

		_D.set(A);
		A.add(A, B);						// A = x1y2 + x2y1, C = z1z2, D = x1y2, E = y1y2 - x1x2
		_D.mul(B);
		_D.mul_m(_d);						// A = x1y2, C = z1z2, D = d * x1x2y1y2, E = y1y2 - x1x2

		A.mul(C);
		_E.mul_m(C);
		C.mul_mm(C);						// A = z1z2 * (x1y2 + x2y1), C = z1z2^2, D = d * x1x2y1y2, E = z1z2 * (y1y2 - x1x2)

		Res<VComplex>::addsub(C, _D);		// A = z1z2 * (x1y2 + x2y1), C = z1z2^2 + d * x1x2y1y2, D = z1z2^2 - d * x1x2y1y2, E = z1z2 * (y1y2 - x1x2)

		A.mul(_D);							// A = x', C = z1z2^2 + d * x1x2y1y2, D = T(z1z2^2 - d * x1x2y1y2), E = z1z2 * (y1y2 - x1x2)
		_E.mul(C);							// A = x', C = T(z1z2^2 + d * x1x2y1y2), D = T(z1z2^2 - d * x1x2y1y2), E = y'
		C.mul_mm(_D);						// A = x', C = z', E = y'
		B.swap(_E);
	}

	// P1 = P1 - P2
	void sub(Point & P1, const Point & P2)
	{
		++_add_count;

		Res<VComplex> & A = P1.x();
		Res<VComplex> & B = P1.y();
		Res<VComplex> & C = P1.z();

		_D.set(P2.x());
		_E.set(P2.y());

		A.to_multiplier();
		_F.set(A);
		_F.mul_t(_D);
		A.mul_t(_E);
		_E.mul_t(B);
		B.mul_mm(_D);
		_E.add(_E, _F);						// E = y1y2 + x1x2

		_D.set(P2.z());
		C.mul(_D);

		_D.set(A);
		A.sub(A, B);						// A = x1y2 - x2y1, E = y1y2 + x1x2
		_D.mul(B);
		_D.mul_m(_d);

		A.mul(C);
		_E.mul_m(C);
		C.mul_mm(C);

		Res<VComplex>::addsub(C, _D);		// A = z1z2 * (x1y2 - x2y1), C = z1z2^2 - d * x1x2y1y2, D = z1z2^2 + d * x1x2y1y2, E = z1z2 * (y1y2 + x1x2)

		A.mul(C);
		_E.mul(_D);
		C.mul_mm(_D);
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
		_d.to_multiplier();
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
			if (ni > 0) add(P, _Pi[(ni - 1) / 2]);
			else if (ni < 0) sub(P, _Pi[(-ni - 1) / 2]);
		}
	}
};

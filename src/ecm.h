/*
Copyright 2021, Yves Gallot

f12ecm is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>
#include <atomic>
#include <thread>
#include <mutex>
#include <chrono>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "pool.h"
#include "prm.h"
#include "ec.h"

inline mpz_class mpz(const uint64_t n)
{
	mpz_class z;
	mp_limb_t * const p_limb = mpz_limbs_write(z.get_mpz_t(), 1);
	p_limb[0] = n;
	mpz_limbs_finish(z.get_mpz_t(), 1);
	return z;
}

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
			den = 2 * p1.y;
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
		bool s = false;
		Point p0; p0.x = m - 2; p0.y = 4;
		Point p = p0;
		for (int b = 63; b >= 0; --b)
		{
			if (s) p = add(p, p, m);

			if ((n & (uint64_t(1) << b)) != 0)
			{
				if (s) p = add(p, p0, m);
				s = true;
			}
		}
		return p;
	}
};

template<typename VComplex>
class ECM
{
private:
	struct deleter { void operator()(const ECM * const p) { delete p; } };

public:
	ECM() : F_12((mpz_class(1) << (1 << 12)) + 1) {}
	virtual ~ECM() {}

	static ECM & getInstance()
	{
		static std::unique_ptr<ECM, deleter> pInstance(new ECM());
		return *pInstance;
	}

public:
	void quit() { _quit = true; }

protected:
	volatile bool _quit = false;

private:
	uint64_t _B1, _B2, _sigma_0;
	size_t _thread_count;
	std::atomic<size_t> _running_threads;
	std::mutex _output_mutex;

	const mpz_class F_12;

public:
	static const size_t D = 128;

private:
	static mpz_class gcd(const size_t i, const Res<VComplex> & x, const mpz_class & n)
	{
		mpz_class t = x.get_z(i); if (t < 0) t += n;
		return ::gcd(t, n);
	}

	static mpz_class sqr_mod(const mpz_class & x, const mpz_class & n) { return (x * x) % n; }
	static mpz_class cub_mod(const mpz_class & x, const mpz_class & n) { mpz_class t = (x * x * x) % n; if (t < 0) t += n; return t; }

	// A such that the group order is divisible by 12 (H. Suyama)
	static mpz_class A_12(const size_t i, const uint64_t sigma, const mpz_class & n, typename EC<VComplex>::Point & P0)
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

	// Montgomery construction of Edwards curve: the group order is divisible by 12
	static mpz_class A_Ed(const size_t i, const uint64_t sigma, const mpz_class & n, typename EC<VComplex>::Point & P0)
	{
		// (s, t) = sigma * [-2, 4] on y^2 = x^3 − 12*x (modulo n)
		const EC_mc::Point P = EC_mc::get(sigma, n);
		const mpz_class s = P.x, t = P.y;

		// Edwards curve
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
		mpz_class y = (-4 * s * a * den);
		if (y < 0) y += n;

		// Montgomery curve
		den = 1 - d; mpz_invert(den.get_mpz_t(), den.get_mpz_t(), n.get_mpz_t());
		const mpz_class c = 2 * den;
		mpz_class A = ((1 + d) * c) % n;
		if (A < 0) A += n;
		mpz_class B = (2 * c) % n;
		if (B < 0) B += n;
		den = 1 - y; mpz_invert(den.get_mpz_t(), den.get_mpz_t(), n.get_mpz_t());
		mpz_class X = ((1 + y) * den) % n;
		if (X < 0) X += n;
		den = x; mpz_invert(den.get_mpz_t(), den.get_mpz_t(), n.get_mpz_t());
		mpz_class Y = X * den;
		if (Y < 0) Y += n;

		P0.set(i, X, mpz_class(1));
		return A;
	}

	void output(const int stage, const std::vector<std::pair<uint64_t, mpz_class>> & sol, const long time)
	{
		static const mpz_class p[6] = {
			mpz_class("114689"), mpz_class("26017793"), mpz_class("63766529"), mpz_class("190274191361"),
			mpz_class("1256132134125569"), mpz_class("568630647535356955169033410940867804839360742060818433") };

		std::stringstream ss;
		ss << "Stage " << stage << ", " << time << " seconds:" << std::endl;
		for (const auto & pair : sol)
		{
			ss << "sigma = " << pair.first << ":";
			mpz_class r = pair.second;
			for (size_t i = 0; i < sizeof(p)/sizeof(mpz_class); ++i)
			{
				if (mpz_divisible_p(r.get_mpz_t(), p[i].get_mpz_t()))
				{
					ss << " " << p[i].get_str();
					r /= p[i];
				}
			}
			if (r != 1) ss << " " << r.get_str() << " !!!";
			ss << std::endl;
		}

		std::lock_guard<std::mutex> guard(_output_mutex);
		std::cout << "                                                  \r" << ss.str();
		std::ofstream logFile("f12.log", std::ios::app);
		if (logFile.is_open())
		{
			logFile << ss.str();
			logFile.flush(); logFile.close();
		}
	}

	void test(const size_t thread_index)
	{
		const uint64_t B1 = _B1, B2 = _B2;

		EC<VComplex> ec(thread_index);
		typename EC<VComplex>::Point P(thread_index);

		typename EC<VComplex>::Point S[D]; for (size_t i = 0; i < D; ++i) S[i].init(thread_index);

		typename EC<VComplex>::Point Se(thread_index), T(thread_index), R(thread_index), Rm(thread_index);
		Res<VComplex> g(thread_index), t1(thread_index), t2(thread_index);

		PseudoPrmGen prmGen;

		const size_t v_size = sizeof(VComplex) / sizeof(Complex);

		for (uint64_t sigma = _sigma_0 + v_size * thread_index; sigma < _sigma_0 + v_size * _thread_count; sigma += v_size * _thread_count)
		{
			const auto clock0 = std::chrono::steady_clock::now();

			for (size_t i = 0; i < v_size; ++i)
			{
				const mpz_class A = A_12(i, sigma + i, F_12, P);
				ec.set(i, A, F_12);
			}
			ec.init();

			for (size_t i = 0; i < 64; ++i) ec.dbl(P, P);

			uint64_t p = prmGen.first(), disp = 100000;
			for (p = prmGen.next(); p <= B1; p = prmGen.next())
			{
				uint64_t m = p; while (m * double(p) <= B1) m *= p;
				ec.mul(P, P, m, p);
				if (_quit) break;
				if ((thread_index == 0) && (p > disp))
				{
					disp += 100000;
					const std::chrono::duration<double> dt = std::chrono::steady_clock::now() - clock0;
					std::cout << "p = " << p << " (B1:" << (100 * p) / B1 << "%), " << std::lrint(dt.count()) << " seconds\r";
				}
			}

			const auto clock1 = std::chrono::steady_clock::now();
			const std::chrono::duration<double> dt1 = clock1 - clock0;

			std::vector<std::pair<uint64_t, mpz_class>> sol1;
			for (size_t i = 0; i < v_size; ++i)
			{
				const mpz_class g1 = gcd(i, P.z(), F_12);
				if (g1 != 1) sol1.push_back(std::make_pair(sigma + i, g1));
			}
			if (_quit) break;

			if (!sol1.empty()) output(1, sol1, std::lrint(dt1.count()));

			ec.dbl(S[0], P); ec.dbl(S[1], S[0]);
			for (size_t d = 2; d < D; ++d) ec.add(S[d], S[d - 1], S[0], S[d - 2]);
			Se.set(S[D - 1]);
			for (size_t d = 0; d < D; ++d) S[d].to_multiplier();

			const uint64_t r_min = p - 2;
			ec.mul(T, P, r_min - 2 * D), ec.mul(R, P, r_min);

			for (size_t i = 0; i < v_size; ++i) g.set_z(i, mpz_class(1));	// TODO set_1

			for (uint64_t r = r_min; r < B2; r += 2 * D)
			{
				Rm.set(R); Rm.to_multiplier();

				for (; p <= r + 2 * D; p = prmGen.next())
				{
					const size_t delta = (p - r) / 2 - 1;
					t1.set(Rm.x()); t1.mul_mm(S[delta].z());
					t2.set(Rm.z()); t2.mul_mm(S[delta].x());
					t1.sub(t1, t2); t1.to_multiplier();
					g.mul_m(t1);
				}

				ec.add(T, R, Se, T);
				T.swap(R);
				if (_quit) break;

				if ((thread_index == 0) && (p > disp))
				{
					disp += 10000000;
					const std::chrono::duration<double> dt = std::chrono::steady_clock::now() - clock0;
					std::cout << "p = " << p << " (B2:" << (100 * p) / B2 << "%), " << std::lrint(dt.count()) << " seconds\r";
				}
			}

			const auto clock2 = std::chrono::steady_clock::now();
			const std::chrono::duration<double> dt2 = clock2 - clock1;

			std::vector<std::pair<uint64_t, mpz_class>> sol2;
			for (size_t i = 0; i < v_size; ++i)
			{
				const mpz_class g2 = gcd(i, g, F_12);
				if (g2 != 1) sol2.push_back(std::make_pair(sigma + i, g2));
			}

			if (!sol2.empty()) output(2, sol2, std::lrint(dt2.count()));
		};

		mainPool.free(thread_index);

		_running_threads--;
	}

public:
	void run(const uint64_t B1, const uint64_t B2, const uint64_t sigma_0, const size_t thread_count, const std::string & ext)
	{
		_B1 = B1; _B2 = B2; _sigma_0 = sigma_0;
		_thread_count = thread_count; _running_threads = 0;

		const size_t v_size = sizeof(VComplex) / sizeof(Complex);

		std::cout << "Testing " << v_size * thread_count << " curves (" << ext << ", " << thread_count << " thread(s)), B1 = "
				<< B1 << ", B2 = " << B2 << "." << std::endl;

		mainPool.init(ECM::D, v_size * sizeof(Complex), thread_count);

		for (size_t i = 0; i < thread_count; ++i)
		{
			_running_threads++;
			std::thread t_test([=] { test(i); });
			t_test.detach();
		}

		while (_running_threads != 0) std::this_thread::sleep_for(std::chrono::seconds(1));

		mainPool.release();
	}
};

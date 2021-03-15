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
#include "ec_m.h"
#include "ec_e.h"

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
	bool _isEdwards;
	size_t _thread_count;
	std::atomic<size_t> _running_threads;
	std::mutex _output_mutex;

	const mpz_class F_12;

public:
	static const size_t D = 128;

private:
	static mpz_class gcd(const size_t i, const Res<VComplex> & x, const mpz_class & n)
	{
		mpz_class t = x.get_z(i) % n; if (t < 0) t += n;
		return ::gcd(t, n);
	}

	static void Edwards2Montgomery(const typename EC_e<VComplex>::Point & P_e, typename EC_m<VComplex>::Point & P_m, const mpz_class & n)
	{
		for (size_t i = 0, v_size = sizeof(VComplex) / sizeof(Complex); i < v_size; ++i)
		{
			mpz_class y = P_e.y().get_z(i) % n; if (y < 0) y += n;
			mpz_class z = P_e.z().get_z(i) % n; if (z < 0) z += n;
			mpz_class X = z + y; if (X >= n) X -= n;
			mpz_class Z = z - y; if (Z < 0) Z += n;
			P_m.set(i, X, Z);
		}
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

		EC_e<VComplex> ec_e(thread_index);
		typename EC_e<VComplex>::Point P_e(thread_index);

		EC_m<VComplex> ec_m(thread_index);
		typename EC_m<VComplex>::Point P_m(thread_index);

		typename EC_m<VComplex>::Point S[D]; for (size_t i = 0; i < D; ++i) S[i].init(thread_index);

		typename EC_m<VComplex>::Point Se(thread_index), T(thread_index), R(thread_index), Rm(thread_index);
		Res<VComplex> g(thread_index), t1(thread_index), t2(thread_index);

		PseudoPrmGen prmGen;

		const size_t v_size = sizeof(VComplex) / sizeof(Complex);

		for (uint64_t sigma = _sigma_0 + v_size * thread_index; sigma < _sigma_0 + v_size * _thread_count; sigma += v_size * _thread_count)
		{
			const auto clock0 = std::chrono::steady_clock::now();

			if (_isEdwards)
			{
				for (size_t i = 0; i < v_size; ++i)
				{
					const mpz_class d = EC_e<VComplex>::d_12(i, sigma + i, F_12, P_e);
					ec_e.set(i, d);

					mpz_class den = 1 - d; mpz_invert(den.get_mpz_t(), den.get_mpz_t(), F_12.get_mpz_t());
					mpz_class A = (2 * (1 + d) * den) % F_12; if (A < 0) A += F_12;
					ec_m.set(i, A, F_12);
				}
				ec_e.init();
				ec_m.init();

				for (size_t i = 0; i < 64; ++i) ec_e.dbl(P_e);
			}
			else
			{
				for (size_t i = 0; i < v_size; ++i)
				{
					const mpz_class A = EC_m<VComplex>::A_12(i, sigma + i, F_12, P_m);
					ec_m.set(i, A, F_12);
				}
				ec_m.init();

				for (size_t i = 0; i < 64; ++i) ec_m.dbl(P_m, P_m);
			}

			uint64_t p = prmGen.first(), disp = 100000;
			for (p = prmGen.next(); p <= B1; p = prmGen.next())
			{
				uint64_t m = p; while (m * double(p) <= B1) m *= p;
				if (_isEdwards) ec_e.mul(P_e, m); else ec_m.mul(P_m, P_m, m, p);
				if (_quit) break;
				if ((thread_index == 0) && (p > disp))
				{
					disp += 100000;
					const std::chrono::duration<double> dt = std::chrono::steady_clock::now() - clock0;
					std::cout << "p = " << p << " (B1:" << (100 * p) / B1 << "%), " << std::lrint(dt.count()) << " seconds\r";
				}
			}

			if (_isEdwards) Edwards2Montgomery(P_e, P_m, F_12);

			const auto clock1 = std::chrono::steady_clock::now();
			const std::chrono::duration<double> dt1 = clock1 - clock0;

			std::vector<std::pair<uint64_t, mpz_class>> sol1;
			for (size_t i = 0; i < v_size; ++i)
			{
				const mpz_class g1 = gcd(i, P_m.z(), F_12);
				if (g1 != 1) sol1.push_back(std::make_pair(sigma + i, g1));
			}
			if (_quit) break;

			if (!sol1.empty()) output(1, sol1, std::lrint(dt1.count()));

			ec_m.dbl(S[0], P_m); ec_m.dbl(S[1], S[0]);
			for (size_t d = 2; d < D; ++d) ec_m.add(S[d], S[d - 1], S[0], S[d - 2]);
			Se.set(S[D - 1]);
			for (size_t d = 0; d < D; ++d) S[d].to_multiplier();

			const uint64_t r_min = p - 2;
			ec_m.mul(T, P_m, r_min - 2 * D), ec_m.mul(R, P_m, r_min);

			g.set1();

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

				ec_m.add(T, R, Se, T);
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
	void run(const uint64_t B1, const uint64_t B2, const uint64_t sigma_0, const bool isEdwards, const size_t thread_count, const std::string & ext)
	{
		_B1 = B1; _B2 = B2; _sigma_0 = sigma_0; _isEdwards = isEdwards;
		_thread_count = thread_count; _running_threads = 0;

		const size_t v_size = sizeof(VComplex) / sizeof(Complex);

		const char * const ctype = isEdwards ? " Edwards" : " Montgomery";
		std::cout << "Testing " << v_size * thread_count << ctype << " curves (" << ext << ", " << thread_count << " thread(s)), B1 = "
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

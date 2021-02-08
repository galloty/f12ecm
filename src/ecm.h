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
		std::cout << ss.str();
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
		Res<VComplex> beta[D]; for (size_t i = 0; i < D; ++i) beta[i].init(thread_index);

		typename EC<VComplex>::Point T(thread_index), R(thread_index);;
		Res<VComplex> g(thread_index), alpha(thread_index), t1(thread_index), t2(thread_index);

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

			uint64_t p = prmGen.first();
			for (; p <= B1; p = prmGen.next())
			{
				uint64_t m = p; while (m * p <= B1) m *= p;
				ec.mul(P, P, m);
				if (_quit) break;
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
			for (size_t d = 2; d < D; ++d) ec.sum(S[d], S[d - 1], S[0], S[d - 2]);

			for (size_t d = 0; d < D; ++d) beta[d].mul(S[d].x(), S[d].z());

			const uint64_t r_min = p - 2;
			ec.mul(T, P, r_min - 2 * D), ec.mul(R, P, r_min);

			for (size_t i = 0; i < v_size; ++i) g.set_z(i, mpz_class(1));	// TODO set_1

			for (uint64_t r = r_min; r < B2; r += 2 * D)
			{
				alpha.mul(R.x(), R.z());

				for (; p <= r + 2 * D; p = prmGen.next())
				{
					const size_t delta = (p - r) / 2 - 1;
					t1.sub(R.x(), S[delta].x()); t2.add(R.z(), S[delta].z());
					t1.mul(t2); t1.sub(t1, alpha); t1.add(t1, beta[delta]);
					g.mul(t1);
				}

				ec.sum(T, R, S[D - 1], T);
				T.swap(R);
				if (_quit) break;
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

		std::cout << "Testing " << v_size * thread_count << " curves with " << ext << " extension..." << std::endl;

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

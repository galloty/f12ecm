/*
Copyright 2021, Yves Gallot

f12ecm is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include <cstdint>
#include <cstdlib>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <thread>
#include <mutex>
#include <atomic>
#include <chrono>

#include "pool.h"
#include "prm.h"
#include "ec.h"

#include <gmpxx.h>

inline mpz_class mpz(const uint64_t n)
{
	mpz_class z;
	mp_limb_t * const p_limb = mpz_limbs_write(z.get_mpz_t(), 1);
	p_limb[0] = n;
	mpz_limbs_finish(z.get_mpz_t(), 1);
	return z;
}

static mpz_class gcd(const Res & x, const mpz_class & n)
{
	mpz_class t = x.get_z(); if (t < 0) t += n;
	return ::gcd(t, n);
}

class ECM
{
private:
	const uint64_t _B1, _B2, _sigma_0;
	const size_t _thread_count;
	std::atomic<size_t> _running_threads;
	std::mutex _output_mutex;

	static const mpz_class F_12;

public:
	static const size_t D = 128;

private:
	static mpz_class cub_mod(const mpz_class & x, const mpz_class & n) { mpz_class t = (x * x * x) % n; if (t < 0) t += n; return t; }

	// A such that the group order is divisible by 12 (H. Suyama)
	static mpz_class A_12(const uint64_t sigma, const mpz_class & n, EC::Point & P0)
	{
		const mpz_class s = mpz(sigma);
		const mpz_class u = s * s - 5, v = 4 * s;
		const mpz_class x = cub_mod(u, n), z = cub_mod(v, n);
		P0.set(x, z);
		mpz_class t = 4 * x * v;
		mpz_invert(t.get_mpz_t(), t.get_mpz_t(), n.get_mpz_t());
		t = (t * cub_mod(v - u, n) * (3 * u + v)) % n;
		t -= 2; if (t < 0) t += n;
		return t;
	}

	void output(const uint64_t sigma, const int level, const uint64_t B, const mpz_class & g)
	{
		static const mpz_class p[6] = {
			mpz_class("114689"), mpz_class("26017793"), mpz_class("63766529"), mpz_class("190274191361"),
			mpz_class("1256132134125569"), mpz_class("568630647535356955169033410940867804839360742060818433") };

		mpz_class r = g;

		std::stringstream ss;
		ss << "sigma = " << sigma << ", B" << level << " = " << B << ":";
		for (size_t i = 0; i < sizeof(p)/sizeof(mpz_class); ++i)
		{
			if (mpz_divisible_p(r.get_mpz_t(), p[i].get_mpz_t()))
			{
				ss << " " << p[i].get_str();
				r /= p[i];
			}
		}
		if (r != 1) ss << " " << r.get_str() << " !!!";

		std::lock_guard<std::mutex> guard(_output_mutex);
		std::cout << ss.str() << std::endl;
		std::ofstream logFile("f12.log", std::ios::app);
		if (logFile.is_open())
		{
			logFile << ss.str() << std::endl;
			logFile.flush(); logFile.close();
		}
	}

	void test(const size_t thread_index)
	{
		const uint64_t B1 = _B1, B2 = _B2;

		EC ec(thread_index);
		EC::Point P(thread_index);

		EC::Point S[D]; for (size_t i = 0; i < D; ++i) S[i].init(thread_index);
		Res beta[D]; for (size_t i = 0; i < D; ++i) beta[i].init(thread_index);

		EC::Point T(thread_index), R(thread_index);;
		Res g(thread_index), alpha(thread_index), t1(thread_index), t2(thread_index);

		PseudoPrmGen prmGen;

		for (uint64_t sigma = _sigma_0 + thread_index; sigma < _sigma_0 + 3 * _thread_count; sigma += _thread_count)
		{
			const mpz_class A = A_12(sigma, F_12, P);
			ec.set(A, F_12);

			uint64_t p = prmGen.first();
			for (; p <= B1; p = prmGen.next())
			{
				uint64_t m = p; while (m * p <= B1) m *= p;
				ec.mul(P, P, m);
			}

			const mpz_class g1 = gcd(P.z(), F_12);
			if (g1 != 1) output(sigma, 1, B1, g1);

			ec.dbl(S[0], P); ec.dbl(S[1], S[0]);
			for (size_t d = 2; d < D; ++d) ec.sum(S[d], S[d - 1], S[0], S[d - 2]);

			for (size_t d = 0; d < D; ++d) beta[d].mul(S[d].x(), S[d].z());

			const uint64_t r_min = p - 2;
			ec.mul(T, P, r_min - 2 * D), ec.mul(R, P, r_min);

			g.set_z(mpz_class(1));

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
			}

			const mpz_class g2 = gcd(g, F_12);
			if (g2 != g1) output(sigma, 2, B2, g2 / g1);
		};

		mainPool.free(thread_index);

		_running_threads--;
	}

public:
	ECM(const uint64_t B1, const uint64_t B2, const uint64_t sigma_0, const size_t thread_count) : _B1(B1), _B2(B2), _sigma_0(sigma_0), _thread_count(thread_count)
	{
		_running_threads = 0;

		for (size_t i = 0; i < thread_count; ++i)
		{
			_running_threads++;
			std::thread t_test([=] { test(i); });
			t_test.detach();
		}

		while (_running_threads != 0)
		{
			std::this_thread::sleep_for(std::chrono::seconds(1));
		}
	}
};

const mpz_class ECM::F_12 = (mpz_class(1) << (1 << 12)) + 1;

int main(int argc, char * argv[])
{
	std::cerr << "f12ecm: factorize the 12th Fermat number with the Elliptic Curve Method" << std::endl;
	std::cerr << " Copyright (c) 2021, Yves Gallot" << std::endl;
	std::cerr << " f12ecm is free source code, under the MIT license." << std::endl << std::endl;
	std::cerr << " Usage: f12ecm <n_threads> <B1> <B2> <sigma_0>" << std::endl;
	std::cerr << "   n_threads: the number of threads (default 3)" << std::endl;
	std::cerr << "   B1: the bound of stage 1 (default 10000)" << std::endl;
	std::cerr << "   B2: the bound of stage 2 (default 100*B1)" << std::endl;
	std::cerr << "   sigma_0: the number of threads (default random)" << std::endl << std::endl;

	std::random_device rd;
	std::uniform_int_distribution<uint64_t> dist(6, uint64_t(-1));

	const size_t n_threads = (argc > 1) ? std::atoi(argv[1]) : 3;
	const uint64_t B1 = (argc > 2) ? std::max(std::atoll(argv[2]), 100ll) : 10000;
	const uint64_t B2 = (argc > 3) ? std::max(std::atoll(argv[3]), 10000ll) : 100 * B1;
	const uint64_t sigma_0 = (argc > 4) ? std::max(std::atoll(argv[4]), 6ll) : dist(rd);

	mainPool.init(ECM::D, n_threads);

	ECM(B1, B2, sigma_0, n_threads);

	return EXIT_SUCCESS;
}
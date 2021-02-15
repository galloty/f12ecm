/*
Copyright 2021, Yves Gallot

f12ecm is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include <cstdint>
#include <cmath>
#include <vector>

class Prob
{
public:
	struct Rho
	{
		double u;
		double k1;	// Dickman function
		double k2;	// Knuth and Trabb Pardo generalization, k = 2
		Rho(const double u, const double k1, const double k2) : u(u), k1(k1), k2(k2) {}
	};

private:
	std::vector<Rho> _rho;

public:
	Prob()
	{
		const size_t n = 65536, u_max = 120, step = 4096;
		const double dt = 1 / double(n);
		double * const f1 = new double[2 * n];
		double * const f2 = new double[2 * n];

		_rho.reserve(u_max * n / step);

		for (size_t i = 0; i < n; ++i)
		{
			const double u = 1 + i * dt;
			f1[i] = 1 - std::log(u); f2[i] = 1;
			if (i % 4096 == 0) _rho.push_back(Rho(u, f1[i], f2[i]));
		}

		double F1 = 2 * n - 1;

		for (size_t j = 2; j <= u_max; ++j)
		{
			double S1 = 0; for (size_t k = 0; k < n - 1; ++k) S1 += f1[k];
			double S2 = 0; for (size_t k = 0; k < n - 1; ++k) S2 += f2[k];

			for (size_t i = 0; i < n; ++i)
			{
				const double u = j + i * dt;

				S1 += f1[i + n - 1] - f1[i];
				f1[i + n] = (f1[i] + 2 * S1) / (2 * n * u - 1);

				S2 += f2[i + n - 1] - f2[i];
				f2[i + n] = (f2[i] + 2 * S2 + f1[i] + F1) / (2 * n * u - 1);
				F1 += 2 * f1[i];

				if (i % 4096 == 0) _rho.push_back(Rho(u, f1[i + n], f2[i + n]));
			}
			for (size_t i = 0; i < n; ++i) f1[i] = f1[i + n];
			for (size_t i = 0; i < n; ++i) f2[i] = f2[i + n];
		}

		delete[] f1;
		delete[] f2;
	}

private:
	double rho(const double u, const int k) const
	{
		if (u <= k) return 1;
		for (size_t i = 1, n = _rho.size(); i < n; ++i)
		{
			if (u < _rho[i].u)
			{
				const double u1 = _rho[i - 1].u, u2 = _rho[i].u;
				const double a = (u1 - u) / (u1 - u2);
				const double r1 = (k == 1) ? _rho[i - 1].k1 : _rho[i - 1].k2;
				const double r2 = (k == 1) ? _rho[i].k1 : _rho[i].k2;
				return r1 * (1 - a) + r2 * a;
			}
		}
		return 0;
	}

public:
	double rho1(const double u) const { return rho(u, 1); }
	double rho2(const double u) const { return rho(u, 2); }
	double F1(const double alpha) const { return rho(1 / alpha, 1); }
	double F2(const double alpha) const { return rho(1 / alpha, 2); }
};

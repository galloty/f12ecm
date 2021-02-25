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
		const size_t n = 262144, u_max = 120, step = 256;
		const double dt = 1 / double(n);
		double * const f1 = new double[2 * n];
		double * const f2 = new double[2 * n];

		_rho.reserve(u_max * n / step);

		for (size_t i = 0; i < n; ++i)
		{
			const double u = 1 + i * dt;
			f1[i] = 1 - std::log(u); f2[i] = 1;
			if (i % step == 0) _rho.push_back(Rho(u, f1[i], f2[i]));
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

				if (i % step == 0) _rho.push_back(Rho(u, f1[i + n], f2[i + n]));
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

	double F_inv(const double p, const double err, const int k) const
	{
		double alpha_min = 0, alpha_max = 1;

		while (true)
		{
			const double alpha = 0.5 * (alpha_min + alpha_max);
			const double p_est = (k == 1) ? F1(alpha) : F2(alpha);
			if (alpha_max - alpha_min < err * p) return alpha;
			if (p_est < p) alpha_min = alpha; else alpha_max = alpha;
		}
	}

public:
	double rho1(const double u) const { return rho(u, 1); }
	double rho2(const double u) const { return rho(u, 2); }
	double F1(const double alpha) const { return rho(1 / alpha, 1); }
	double F2(const double alpha) const { return rho(1 / alpha, 2); }
	double F1_inv(const double p, const double err = 1e-5) const { return F_inv(p, err, 1); }
	double F2_inv(const double p, const double err = 1e-5) const { return F_inv(p, err, 2); }

	double sigma(const double u, const double v) const
	{
		double s = rho1(v);
		const double a = v + 1 - v / u, b = v;
		for (int k = 0, n = std::max(1024, int((b - a) * 1024)); k <= n; ++k)
		{
			const double dt = (b - a) / double(n), t = a + k * dt;
			double f = rho1(t - 1) * dt / (v + 1 - t);
			if ((k == 0) || (k == n)) f *= 0.5;
			s += f;
		}
		return s;
	}

	double G(const double alpha, const double beta) const { return (beta < alpha) ? sigma(1 / alpha, 1 / beta) : F1(alpha); }

	double B2(const int digits) const
	{
		const double logP = (digits - 1) * log(10.0) + log(5), eps = 1e-9;
		double B2_min = 100, B2_max = 1e15;
		while (true)
		{
			const double B2 = std::sqrt(B2_min * B2_max);
			const double B2_m = B2 * (1 - eps), logB2_m = log(B2_m);
			const double f_m = F1(logB2_m / logP) * logB2_m / B2_m;
			const double B2_p = B2 * (1 + eps), logB2_p = log(B2_p);
			const double f_p = F1(logB2_p / logP) * logB2_p / B2_p;
			if (f_p > f_m) B2_min = B2; else B2_max = B2;
			if (fabs(B2_max - B2_min) / B2_max < eps) break;
		}
		return 0.5 * (B2_min + B2_max);
	}

	double B1(const int digits, const double B2) const
	{
		const double logP = (digits - 1) * log(10.0) + log(5), eps = 1e-9;
		const double logB2 = log(B2), alpha = logB2 / logP;
		double B1_min = B2 / 1000, B1_max = B2;
		while (true)
		{
			const double B1 = std::sqrt(B1_min * B1_max);
			const double B1_m = B1 * (1 - eps), logB1_m = log(B1_m);
			const double f_m = G(alpha, logB1_m / logP) / (B1_m + B2 / logB2);
			const double B1_p = B1 * (1 + eps), logB1_p = log(B1_p);
			const double f_p = G(alpha, logB1_p / logP) / (B1_p + B2 / logB2);
			if (f_p > f_m) B1_min = B1; else B1_max = B1;
			if (fabs(B1_max - B1_min) / B1_max < eps) break;
		}
		return 0.5 * (B1_min + B1_max);
	}

	double B2(const int digits, const double B1, const double B2_0) const
	{
		const double logP = (digits - 1) * log(10.0) + log(5), eps = 1e-9;
		const double beta = log(B1) / logP;
		double B2_min = 0.5 * B2_0, B2_max = 2 * B2_0;
		while (true)
		{
			const double B2 = 0.5 * (B2_min + B2_max);
			const double B2_m = B2 * (1 - eps), logB2_m = log(B2_m);
			const double f_m = G(logB2_m / logP, beta) / (B1 + B2_m / logB2_m);
			const double B2_p = B2 * (1 + eps), logB2_p = log(B2_p);
			const double f_p = G(logB2_p / logP, beta) / (B1 + B2_p / logB2_p);
			if (f_p > f_m) B2_min = B2; else B2_max = B2;
			if (fabs(B2_max - B2_min) / B2_max < eps) break;
		}
		return 0.5 * (B2_min + B2_max);
	}

	void B12(const int digits, double & rB1, double & rB2, int & n) const
	{
		double bB2 = B2(digits), bB1 = B1(digits, bB2);
		double err = 1e100;
		while (true)
		{
			const double nB2 = B2(digits, bB1, bB2);
			const double nB1 = B1(digits, nB2);
			const double e = (fabs(nB2 - bB2) / bB2) * (fabs(nB1 - bB1) / bB1);
			if (e > err) break;
			err = e; bB2 = nB2; bB1 = nB1;
		}
		rB2 = bB2; rB1 = bB1;
		const double logP = (digits - 1) * log(10.0) + log(5);
		const double alpha = log(rB2) / logP, beta = log(rB1) / logP;
		const double g = G(alpha, beta);
		n = std::lrint(1 / g);
	}
};

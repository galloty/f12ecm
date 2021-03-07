/*
Copyright 2021, Yves Gallot

f12ecm is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "prm.h"

// PRAC algorithm: Montgomery, Peter L. (1983). "Evaluating recurrences of form X_{m+n} = f(X_m, X_n, X_{m-n}) via Lucas Chains", unpublished.

class PRAC
{
private:
	size_t _rule[10] = { 0 };
	size_t _sqr_count = 0, _mul_count = 0;;
	int64_t _A, _B, _C, _T1, _T2;

	const double phi = 1.6180339887498948482045868343656381177;
	const double _alpha[10] = { 1 / phi, 5 / (phi + 7), 1429 / (phi + 2311), 3739 / (6051 - phi), 79 / (129 - phi),
								31 / (phi + 49), 209 / (phi + 337), 11 / (19 - phi), 545 / (883 - phi), 1 / (3 - phi) };

private:
	static void add(int64_t & R, const int64_t & A, const int64_t & B, const int64_t & C)
	{
		const int64_t T = A + B;
		if (T + C != 2 * A) std::cout << "error add: " << A << ", " << B << ", " << C << std::endl;
		R = T;
	}

	static void dbl(int64_t & R, const int64_t & A) { R = A + A; }

	static void check(const int64_t A, const int64_t B, const int64_t C, const uint64_t a, const uint64_t b, const int64_t P, const int i)
	{
		if (A != int64_t(a) * P) std::cout << "error A (" << i << "): " << A << " != " << a << " * " << P << std::endl;
		if (B != int64_t(b) * P) std::cout << "error B (" << i << "): " << B << " != " << b << " * " << P << std::endl;
		if (C != A - B) std::cout << "error C (" << i << "): " << A << ", " << B << ", " << C << std::endl;
		if (C == 0) std::cout << "warning (" << i << ")" << std::endl;
	}

public:
	static size_t cost(const uint64_t n, const uint64_t p, const double alpha, size_t & sqr_count, size_t & mul_count)
	{
		size_t dbl_count = 0, add_count = 0;

		uint64_t d = n;

		while (d != 1)
		{
			const uint64_t r = (p == 0) ? std::llrint(d * alpha) : (d / p) * std::llrint(p * alpha);

			uint64_t e = 2 * r - d;
			d -= r;
			dbl_count += 1;

			while (d != e)
			{
				if (d < e) std::swap(d, e);

				if ((4 * d <= 5 * e) && ((d + e) % 3 == 0))	// #1
				{
					const uint64_t t = (2 * e - d) / 3; d = (2 * d - e) / 3; e = t;
					add_count += 3;
				}
				else if ((4 * d <= 5 * e) && ((d - e) % 6 == 0))	// #2
				{
					d = (d - e) / 2;
					dbl_count += 1; add_count += 1;
				}
				else if (d < 4 * e)	// #3
				{
					d -= e;
					add_count += 1;
				}
				else if ((d - e) % 2 == 0)	// #4
				{
					d = (d - e) / 2;
					dbl_count += 1; add_count += 1;
				}
				else if (d % 2 == 0)	// #5
				{
					d /= 2;
					dbl_count += 1; add_count += 1;
				}
				else if (d % 3 == 0)	// #6
				{
					d = d / 3 - e;
					dbl_count += 1; add_count += 3;
				}
				else if ((d + e) % 3 == 0)	// #7
				{
					d = (d - 2*e) / 3;
					dbl_count += 1; add_count += 3;
				}
				else if ((d - e) % 3 == 0)	// #8
				{
					d = (d - e) / 3;
					dbl_count += 1; add_count += 3;
				}
				else if (e % 2 == 0)	// #9
				{
					e = e / 2;
					dbl_count += 1; add_count += 1;
				}
				else
				{
					throw std::runtime_error("PRAC");
				}
			}

			add_count += 1;
		}

		sqr_count = 3 * dbl_count + 2 * add_count;
		mul_count = 2 * dbl_count + 4 * add_count;
		return 12 * dbl_count + 16 * add_count;
	}

public:
	void valid(const uint64_t n, const uint64_t p, const int64_t P)
	{
		if (n == 0)
		{
			size_t s = 0;
			for (size_t i = 1; i <= 9; ++i)
			{
				std::cout << " R" << i << ": " << _rule[i];
				s += _rule[i];
			}
			std::cout << std::endl << s;
			for (size_t i = 1; i <= 9; ++i)
			{
				std::cout << " R" << i << ": " << _rule[i] * 100.0 / s;
			}
			const size_t transform_count = 2 * _sqr_count + 3 * _mul_count;
			double f = std::log(2) / p;	// p is B if n = 0. We have log B# ~ B => log_2 B# ~ B / log(2)
			std::cout << std::endl << _sqr_count << " squares and " << _mul_count << " muls: " << transform_count << " transforms" << std::endl;
			std::cout << "Per bit: " << _sqr_count * f << " squares and " << _mul_count * f << " muls: " << transform_count * f << " transforms" << std::endl;
			return;
		}

		size_t c_min = size_t(-1), best_i = 0, sqr_min = 0, mul_min = 0;
		for (size_t i = 0; i < 10; ++i)
		{
			size_t sqr_count, mul_count;
			const size_t c = cost(n, p, _alpha[i], sqr_count, mul_count);
			if (c < c_min)
			{
				c_min = c; best_i = i;
				sqr_min = sqr_count; mul_min = mul_count;
			}
		}
		_sqr_count += sqr_min; _mul_count += mul_min;

		_A = P;
		uint64_t a = 1, d = n;

		size_t cnt = 0, cnt_loop = 0;

		while (d != 1)
		{
			const uint64_t r = (d / p) * std::llrint(p * _alpha[best_i]);

			uint64_t b = a, e = 2 * r - d;
			a *= 2; d -= r;
			_B = _A; _C =  _A; dbl(_A, _A);

			while (d != e)
			{
				if (d < e)
				{
					std::swap(_A, _B); _C = -_C;
					std::swap(a, b); std::swap(d, e);
				}

				int i = 0;

				// e <= d
				if ((4 * d <= 5 * e) && ((d + e) % 3 == 0))	// #1
				{
					i = 1;
					const uint64_t f = (2 * e - d) / 3; d = (2 * d - e) / 3; e = f;
					add(_T1, _A, _B, _C); add(_T2, _T1, _A, _B); add(_B, _T1, _B, _A); std::swap(_A, _T2);
					const uint64_t t = a + b; a += t; b += t;
				}
				else if ((4 * d <= 5 * e) && ((d - e) % 6 == 0))	// #2
				{
					i = 2;
					d = (d - e) / 2;
					add(_B, _A, _B, _C); dbl(_A, _A);
					b += a; a *= 2;
				}
				else if (d < 4 * e)	// #3
				{
					i = 3;
					d -= e;
					add(_C, _A, _B, _C); std::swap(_B, _C); _C = -_C;
					b += a;
				}
				else if ((d - e) % 2 == 0)	// #4
				{
					i = 4;
					d = (d - e) / 2;
					add(_B, _A, _B, _C); dbl(_A, _A);
					b += a; a *= 2;
				}
				else if (d % 2 == 0)	// #5
				{
					i = 5;
					d /= 2;
					add(_C, _A, _C, _B); dbl(_A, _A);
					a *= 2;
				}
				else if (d % 3 == 0)	// #6
				{
					i = 6;
					d = d / 3 - e;
					dbl(_T1, _A); add(_T2, _A, _B, _C);
					add(_A, _T1, _A, _A); add(_C, _T1, _T2, _C); std::swap(_B, _C); _C = -_C;
					a *= 3; b += a;
				}
				else if ((d + e) % 3 == 0)	// #7
				{
					i = 7;
					d = (d - 2*e) / 3;
					add(_T1, _A, _B, _C); dbl(_T2, _A);
					add(_B, _T1, _A, _B); add(_A, _T2, _A, _A);
					b += 2 * a; a *= 3;
				}
				else if ((d - e) % 3 == 0)	// #8
				{
					i = 8;
					d = (d - e) / 3;
					add(_T1, _A, _B, _C); add(_C, _A, _C, _B);
					std::swap(_B, _T1);
					dbl(_T1, _A); add(_A, _T1, _A, _A);
					b += a; a *= 3;
				}
				else if (e % 2 == 0)	// #9
				{
					i = 9;
					e = e / 2;
					add(_C, _C, -_B, _A); dbl(_B, _B);
					b *= 2;
				}
				else
				{
					throw std::runtime_error("PRAC");
				}

				check(_A, _B, _C, a, b, P, i);
				if ((_A < 0) | (_B < 0)) std::cout << "error: " << _A << ", " << _B << ", " << _C << std::endl;

				// std::cout << cnt << ", " << rule << ", " << _A << ", " << _B << ", " << _C << ": " << a << ", " << b << ", " << d << ", " << e << std::endl;
				cnt++;
				_rule[i]++;
			}

			add(_C, _A, _B, _C); std::swap(_A, _C);
			a += b;
			check(_A, _B, _C, a, b, P, 0);

			// if ((d != 1) || (cnt_loop != 0))
			// {
			// 	std::cout << cnt_loop << ", " << cnt << ", " << n << ", " << _A << ", " << _B << ", " << _C << ": " << a << ", " << b << ", " << d << ", " << e << std::endl;
			// }
			cnt++;
			cnt_loop++;
		}

		if (_A != int64_t(n) * P) std::cout << "error" << _A << " != " << n << " * " << P << std::endl;
	}

	void check()
	{
		PseudoPrmGen prmGen;

		uint64_t B = 100000000;
		uint64_t p = prmGen.first();
		for (p = prmGen.next(); p <= B; p = prmGen.next())
		{
			uint64_t m = p; while (m * p <= B) m *= p;
			valid(m, p, 10);
		}
		valid(0, B, 0);
	}
};

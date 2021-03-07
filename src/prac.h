/*
Copyright 2021, Yves Gallot

f12ecm is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>
#include <cmath>
#include <iostream>

#include "prm.h"

// PRAC algorithm: Montgomery, Peter L. (1983). "Evaluating recurrences of form X_{m+n} = f(X_m, X_n, X_{m-n}) via Lucas Chains", unpublished.

class PRAC
{
private:
	size_t _rule[10] = { 0 };
	int64_t _A, _B, _C, _T1, _T2;

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
			std::cout << std::endl;
			return;
		}

		_A = P;
		uint64_t a = 1, d = n;

		size_t cnt = 0, cnt_loop = 0;

		while (d != 1)
		{
			const uint64_t r = (d / p) * std::llrint(p * ((std::sqrt(5) - 1) / 2));

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

		uint64_t B = 1000000000;
		uint64_t p = prmGen.first();
		for (p = prmGen.next(); p <= B; p = prmGen.next())
		{
			uint64_t m = p; while (m * p <= B) m *= p;
			valid(m, p, 10);
		}
		valid(0, 0, 0);
	}
};

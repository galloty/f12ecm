/*
Copyright 2021, Yves Gallot

f12ecm is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>

#include "complex.h"

struct Transform
{
	static void forward4(const size_t m, Complex * const z, const Complex * const w)
	{
		for (size_t k = 0; k < m; ++k)
		{
			Complex & z0 = z[k + 0 * m]; Complex & z1 = z[k + 1 * m]; Complex & z2 = z[k + 2 * m]; Complex & z3 = z[k + 3 * m];
			const Complex u0 = z0, u2 = z2.mulW(w[0]), u1 = z1.mulW(w[1]), u3 = z3.mulW(w[2]);
			const Complex v0 = u0 + u2, v2 = u0 - u2, v1 = u1 + u3, v3 = u1 - u3;
			z0 = v0 + v1; z1 = v0 - v1; z2 = v2.addi(v3); z3 = v2.subi(v3);
		}
	}

	static void backward4(const size_t m, Complex * const z, const Complex * const w)
	{
		for (size_t k = 0; k < m; ++k)
		{
			Complex & z0 = z[k + 0 * m]; Complex & z1 = z[k + 1 * m]; Complex & z2 = z[k + 2 * m]; Complex & z3 = z[k + 3 * m];
			const Complex v0 = z0, v1 = z1, v2 = z2, v3 = z3;
			const Complex u0 = v0 + v1, u1 = v0 - v1, u2 = v2 + v3, u3 = v2 - v3;
			z0 = u0 + u2; z2 = Complex(u0 - u2).mulWconj(w[0]); z1 = u1.subi(u3).mulWconj(w[1]); z3 = u1.addi(u3).mulWconj(w[2]);
		}
	}

	static void forward4_0(const size_t m, Complex * const z)
	{
		const double csqrt2_2 = 0.707106781186547524400844362104849039284835937688;
		const Complex cs2pi_16 = Complex(0.92387953251128675612818318939678828682, 0.41421356237309504880168872420969807857);

		for (size_t k = 0; k < m; ++k)
		{
			Complex & z0 = z[k + 0 * m]; Complex & z1 = z[k + 1 * m]; Complex & z2 = z[k + 2 * m]; Complex & z3 = z[k + 3 * m];
			const Complex u0 = z0, u1 = z1, t2 = z2, t3 = z3;
			const Complex u2 = Complex(t2.real() - t2.imag(), t2.real() + t2.imag()) * csqrt2_2;
			const Complex u3 = Complex(t3.real() - t3.imag(), t3.real() + t3.imag()) * csqrt2_2;
			const Complex v0 = u0 + u2, v2 = u0 - u2, v1 = Complex(u1 + u3).mulW(cs2pi_16), v3 = Complex(u1 - u3).mulW(cs2pi_16);
			z0 = v0 + v1; z1 = v0 - v1; z2 = v2.addi(v3); z3 = v2.subi(v3);
		}
	}

	static void backward4_0(const size_t m, Complex * const z)
	{
		const double csqrt2_2 = 0.707106781186547524400844362104849039284835937688;
		const Complex cs2pi_16 = Complex(0.92387953251128675612818318939678828682, 0.41421356237309504880168872420969807857);

		for (size_t k = 0; k < m; ++k)
		{
			Complex & z0 = z[k + 0 * m]; Complex & z1 = z[k + 1 * m]; Complex & z2 = z[k + 2 * m]; Complex & z3 = z[k + 3 * m];
			const Complex v0 = z0, v1 = z1, v2 = z2, v3 = z3;
			const Complex u0 = v0 + v1, u1 = Complex(v0 - v1).mulWconj(cs2pi_16), u2 = v2 + v3, u3 = Complex(v2 - v3).mulWconj(cs2pi_16);
			z0 = u0 + u2; const Complex t2 = (u0 - u2) * csqrt2_2; z2 = Complex(t2.imag() + t2.real(), t2.imag() - t2.real());
			z1 = u1.subi(u3); const Complex t3 = u1.addi(u3) * csqrt2_2; z3 = Complex(t3.imag() + t3.real(), t3.imag() - t3.real());
		}
	}

	static void square2(Complex * const z, const Complex & w)
	{
		const Complex u0 = z[0], u1 = z[1];
		z[1] = (u0 + u0) * u1;
		z[0] = u0.sqr() + u1.sqr().mulW(w);

		const Complex u2 = z[2], u3 = z[3];
		z[3] = (u2 + u2) * u3;
		z[2] = u2.sqr() - u3.sqr().mulW(w);

		const Complex u4 = z[4], u5 = z[5];
		z[5] = (u4 + u4) * u5;
		z[4] = u4.sqr().addi(u5.sqr().mulW(w));

		const Complex u6 = z[6], u7 = z[7];
		z[7] = (u6 + u6) * u7;
		z[6] = u6.sqr().subi(u7.sqr().mulW(w));
	}

	static void mul2(Complex * const z, const Complex * const zr, const Complex & w)
	{
		const Complex u0 = z[0], u1 = z[1], v0 = zr[0], v1 = zr[1];
		z[1] = u0 * v1 + u1 * v0;
		z[0] = u0 * v0 + Complex(u1 * v1).mulW(w);

		const Complex u2 = z[2], u3 = z[3], v2 = zr[2], v3 = zr[3];
		z[3] = u2 * v3 + u3 * v2;
		z[2] = u2 * v2 - Complex(u3 * v3).mulW(w);

		const Complex u4 = z[4], u5 = z[5], v4 = zr[4], v5 = zr[5];
		z[5] = u4 * v5 + u5 * v4;
		z[4] = Complex(u4 * v4).addi(Complex(u5 * v5).mulW(w));

		const Complex u6 = z[6], u7 = z[7], v6 = zr[6], v7 = zr[7];
		z[7] = u6 * v7 + u7 * v6;
		z[6] = Complex(u6 * v6).subi(Complex(u7 * v7).mulW(w));
	}

	// -65536 < z[i] < 65536
	static void norm(Complex * const z, const size_t m)
	{
		// static double max_err = 0;

		Complex f = Complex(0.0, 0.0);	//, e = Complex(0.0, 0.0);

		for (size_t k = 0; k < m; ++k)
		{
			const Complex r = z[k] * (2.0 / m), r_i = r.round();
			// e = e.max_abs(r - r_i);
			const Complex fi = f + r_i;
			f = Complex(fi * (1 / 65536.0)).round();
			z[k] = fi - f * 65536.0;
		}

		bool is_zero = f.is_zero();
		// const double err = ::fmax(e.real(), e.imag());
		// if (err > max_err) { max_err = err; std::cout << err << std::endl; }

		while (!is_zero)
		{
			f = f.rotate();

			for (size_t k = 0; k < m; ++k)
			{
				const Complex fi = f + z[k];
				f = Complex(fi * (1 / 65536.0)).round();
				is_zero = f.is_zero();
				z[k] = fi - f * 65536.0;
				if (is_zero) break;
			}
		}
	}
};

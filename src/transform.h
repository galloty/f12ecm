/*
Copyright 2021, Yves Gallot

f12ecm is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>

#include "complex.h"

template<typename VComplex>
struct Transform
{
	finline static void forward4(const size_t m, VComplex * __restrict__ const z, const Complex * __restrict__ const w)
	{
		for (size_t k = 0; k < m; ++k)
		{
			VComplex & z0 = z[k + 0 * m]; VComplex & z1 = z[k + 1 * m]; VComplex & z2 = z[k + 2 * m]; VComplex & z3 = z[k + 3 * m];
			const VComplex u0 = z0, u2 = z2.mulW(w[0]), u1 = z1.mulW(w[1]), u3 = z3.mulW(w[2]);
			const VComplex v0 = u0 + u2, v2 = u0 - u2, v1 = u1 + u3, v3 = u1 - u3;
			z0 = v0 + v1; z1 = v0 - v1; z2 = v2.addi(v3); z3 = v2.subi(v3);
		}
	}

	finline static void backward4(const size_t m, VComplex * __restrict__ const z, const Complex * __restrict__ const w)
	{
		for (size_t k = 0; k < m; ++k)
		{
			VComplex & z0 = z[k + 0 * m]; VComplex & z1 = z[k + 1 * m]; VComplex & z2 = z[k + 2 * m]; VComplex & z3 = z[k + 3 * m];
			const VComplex v0 = z0, v1 = z1, v2 = z2, v3 = z3;
			const VComplex u0 = v0 + v1, u1 = v0 - v1, u2 = v2 + v3, u3 = v2 - v3;
			z0 = u0 + u2; z2 = VComplex(u0 - u2).mulWconj(w[0]); z1 = u1.subi(u3).mulWconj(w[1]); z3 = u1.addi(u3).mulWconj(w[2]);
		}
	}

	finline static void forward4_0(const size_t m, VComplex * __restrict__ const z)
	{
		const double csqrt2_2 = 0.707106781186547524400844362104849039284835937688;
		const Complex cs2pi_16 = Complex(0.92387953251128675612818318939678828682, 0.41421356237309504880168872420969807857);

		for (size_t k = 0; k < m; ++k)
		{
			VComplex & z0 = z[k + 0 * m]; VComplex & z1 = z[k + 1 * m]; VComplex & z2 = z[k + 2 * m]; VComplex & z3 = z[k + 3 * m];
			const VComplex u0 = z0, u1 = z1, t2 = z2, t3 = z3;
			const VComplex u2 = VComplex(t2.real() - t2.imag(), t2.real() + t2.imag()) * csqrt2_2;
			const VComplex u3 = VComplex(t3.real() - t3.imag(), t3.real() + t3.imag()) * csqrt2_2;
			const VComplex v0 = u0 + u2, v2 = u0 - u2, v1 = VComplex(u1 + u3).mulW(cs2pi_16), v3 = VComplex(u1 - u3).mulW(cs2pi_16);
			z0 = v0 + v1; z1 = v0 - v1; z2 = v2.addi(v3); z3 = v2.subi(v3);
		}
	}

	finline static void backward4_0(const size_t m, VComplex * __restrict__ const z)
	{
		const double csqrt2_2 = 0.707106781186547524400844362104849039284835937688;
		const Complex cs2pi_16 = Complex(0.92387953251128675612818318939678828682, 0.41421356237309504880168872420969807857);

		for (size_t k = 0; k < m; ++k)
		{
			VComplex & z0 = z[k + 0 * m]; VComplex & z1 = z[k + 1 * m]; VComplex & z2 = z[k + 2 * m]; VComplex & z3 = z[k + 3 * m];
			const VComplex v0 = z0, v1 = z1, v2 = z2, v3 = z3;
			const VComplex u0 = v0 + v1, u1 = VComplex(v0 - v1).mulWconj(cs2pi_16), u2 = v2 + v3, u3 = VComplex(v2 - v3).mulWconj(cs2pi_16);
			z0 = u0 + u2; const VComplex t2 = (u0 - u2) * csqrt2_2; z2 = VComplex(t2.imag() + t2.real(), t2.imag() - t2.real());
			z1 = u1.subi(u3); const VComplex t3 = u1.addi(u3) * csqrt2_2; z3 = VComplex(t3.imag() + t3.real(), t3.imag() - t3.real());
		}
	}

	finline static void square2(VComplex * __restrict__ const z, const Complex & w)
	{
		const VComplex u0 = z[0], u1 = z[1];
		z[1] = (u0 + u0) * u1;
		z[0] = u0.sqr() + u1.sqr().mulW(w);

		const VComplex u2 = z[2], u3 = z[3];
		z[3] = (u2 + u2) * u3;
		z[2] = u2.sqr() - u3.sqr().mulW(w);

		const VComplex u4 = z[4], u5 = z[5];
		z[5] = (u4 + u4) * u5;
		z[4] = u4.sqr().addi(u5.sqr().mulW(w));

		const VComplex u6 = z[6], u7 = z[7];
		z[7] = (u6 + u6) * u7;
		z[6] = u6.sqr().subi(u7.sqr().mulW(w));
	}

	finline static void mul2(VComplex * __restrict__ const z, const VComplex * __restrict__ const zr, const Complex & w)
	{
		const VComplex u0 = z[0], u1 = z[1], v0 = zr[0], v1 = zr[1];
		z[1] = u0 * v1 + u1 * v0;
		z[0] = u0 * v0 + VComplex(u1 * v1).mulW(w);

		const VComplex u2 = z[2], u3 = z[3], v2 = zr[2], v3 = zr[3];
		z[3] = u2 * v3 + u3 * v2;
		z[2] = u2 * v2 - VComplex(u3 * v3).mulW(w);

		const VComplex u4 = z[4], u5 = z[5], v4 = zr[4], v5 = zr[5];
		z[5] = u4 * v5 + u5 * v4;
		z[4] = VComplex(u4 * v4).addi(VComplex(u5 * v5).mulW(w));

		const VComplex u6 = z[6], u7 = z[7], v6 = zr[6], v7 = zr[7];
		z[7] = u6 * v7 + u7 * v6;
		z[6] = VComplex(u6 * v6).subi(VComplex(u7 * v7).mulW(w));
	}

	// -65536 < z[i] < 65536
	finline static void norm(VComplex * __restrict__ const z, const size_t m, const int multiplier)
	{
		const double base = 65536, base_inv = 1 / 65536.0;

		// static double max_err = 0;

		VComplex * const z0 = &z[0 * m / 4];
		VComplex * const z1 = &z[1 * m / 4];
		VComplex * const z2 = &z[2 * m / 4];
		VComplex * const z3 = &z[3 * m / 4];

		VComplex f0 = VComplex(0.0, 0.0), f1 = VComplex(0.0, 0.0), f2 = VComplex(0.0, 0.0), f3 = VComplex(0.0, 0.0);
		// VComplex e = VComplex(0.0, 0.0);

		const double c = multiplier * (2.0 / m);
		for (size_t k = 0; k < m / 4; ++k)
		{
			const VComplex r0 = z0[k] * c, ri0 = r0.round(), fi0 = f0 + ri0;
			f0 = VComplex(fi0 * base_inv).round();
			z0[k] = fi0 - f0 * base;

			const VComplex r1 = z1[k] * c, ri1 = r1.round(), fi1 = f1 + ri1;
			f1 = VComplex(fi1 * base_inv).round();
			z1[k] = fi1 - f1 * base;

			const VComplex r2 = z2[k] * c, ri2 = r2.round(), fi2 = f2 + ri2;
			f2 = VComplex(fi2 * base_inv).round();
			z2[k] = fi2 - f2 * base;

			const VComplex r3 = z3[k] * c, ri3 = r3.round(), fi3 = f3 + ri3;
			f3 = VComplex(fi3 * base_inv).round();
			z3[k] = fi3 - f3 * base;

			// e = e.max(VComplex(r0 - ri0).abs()); e = e.max(VComplex(r1 - ri1).abs());
			// e = e.max(VComplex(r2 - ri2).abs()); e = e.max(VComplex(r3 - ri3).abs());
		}

		bool is_zero = f0.is_zero() & f1.is_zero() & f2.is_zero() & f3.is_zero();
		// const double err = e.max();	// SSE4 only
		// if (err > max_err) { max_err = err; std::cout << err << std::endl; }

		while (!is_zero)
		{
			const VComplex t = f3; f3 = f2; f2 = f1; f1 = f0; f0 = t.rotate();

			for (size_t k = 0; k < m / 4; ++k)
			{
				const VComplex fi0 = f0 + z0[k], fi1 = f1 + z1[k], fi2 = f2 + z2[k], fi3 = f3 + z3[k];
				f0 = VComplex(fi0 * base_inv).round(); f1 = VComplex(fi1 * base_inv).round();
				f2 = VComplex(fi2 * base_inv).round(); f3 = VComplex(fi3 * base_inv).round();
				is_zero = f0.is_zero() & f1.is_zero() & f2.is_zero() & f3.is_zero();
				z0[k] = fi0 - f0 * base; z1[k] = fi1 - f1 * base; z2[k] = fi2 - f2 * base; z3[k] = fi3 - f3 * base;
				if (is_zero) break;
			}
		}
	}
};

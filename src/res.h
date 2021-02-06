/*
Copyright 2021, Yves Gallot

f12ecm is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include "pool.h"
#include "transform.h"

#include <gmpxx.h>

class Res : public Transform
{
private:
	static const size_t _n = 256 / 2;	// F_12 = 65536^256 + 1
	Complex * _z;

public:
	Res() : _z(nullptr) {}	// init(thread_index) must be called
	Res(const size_t thread_index) { init(thread_index); }
	Res(const Res &) = delete;
	Res & operator=(const Res &) = delete;

	void init(const size_t thread_index)
	{
		_z = static_cast<Complex *>(mainPool.alloc(thread_index, _n * sizeof(Complex)));
	}

	mpz_class get_z() const
	{
		const Complex * const z = _z;

		mpz_class r = 0;
		for (size_t k = 0; k < _n; ++k)
		{
			const int h = int(z[_n - 1 - k].imag());
			r = (r << 16) + h;
		}
		for (size_t k = 0; k < _n; ++k)
		{
			const int l = int(z[_n - 1 - k].real());
			r = (r << 16) + l;
		}
		return r;
	}

	void set_z(const mpz_class & r)
	{
		const size_t r_size = mpz_size(r.get_mpz_t());
		// if ((r < 0) ||(r_size * sizeof(uint64_t) > 2 * _n * sizeof(uint16_t)))
		// {
		// 	std::cout << mpz_class((mpz_class(1) << (1 << 12)) + 1).get_str() << std::endl << r.get_str() << std::endl;
		// 	throw;
		// }

		uint16_t s[2 * _n];
		const mp_limb_t * r_ptr = mpz_limbs_read(r.get_mpz_t());
		for (size_t i = 0; i < r_size; ++i)
		{
			const uint64_t r_i = r_ptr[i];
			for (size_t j = 0; j < 4; ++j) s[4 * i + j] = uint16_t(r_i >> (16 * j));
		}
		for (size_t i = 4 * r_size; i < 256; ++i) s[i] = 0;

		Complex * const z = _z;
		for (size_t k = 0; k < _n; ++k)
		{
			z[k] = Complex(double(s[k]), double(s[k + _n]));
		}
	}

	void set(const Res & x)
	{
		const Complex * const zr = x._z;
		Complex * const z = _z;
		for (size_t k = 0; k < _n; ++k) z[k] = zr[k];
	}

	void swap(Res & x)
	{
		std::swap(_z, x._z);
	}

	void add(const Res & x, const Res & y)
	{
		const Complex * const zx = x._z;
		const Complex * const zy = y._z;
		Complex * const z = _z;
		for (size_t k = 0; k < _n; ++k) z[k] = zx[k] + zy[k];
	}

	void sub(const Res & x, const Res & y)
	{
		const Complex * const zx = x._z;
		const Complex * const zy = y._z;
		Complex * const z = _z;
		for (size_t k = 0; k < _n; ++k) z[k] = zx[k] - zy[k];
	}

	static void addsub(Res & x, Res & y)
	{
		Complex * const zx = x._z;
		Complex * const zy = y._z;
		for (size_t k = 0; k < _n; ++k) { const Complex u = zx[k], v = zy[k]; zx[k] = u + v; zy[k] = u - v; }
	}

	void sqr()
	{
		const Complex * const w123 = mainPool.w123;

		forward4_0(_n / 4, _z);
		for (size_t j = 0; j < 4; ++j)
		{
			Complex * const z = &_z[j * _n / 4];

			const Complex * const w = &w123[3 * (4 + j)];
			const Complex * const w4 = &w123[3 * 4 * (4 + j)];
			const Complex * const ws = &w123[3 * 2 * 4 * (4 + j)];
			forward4(8, z, w);
			for (size_t i = 0; i < 4; ++i) forward4(2, &z[8 * i], &w4[3 * i]);
			for (size_t i = 0; i < 4; ++i) square2(&z[8 * i], ws[3 * 2 * i]);
			for (size_t i = 0; i < 4; ++i) backward4(2, &z[8 * i], &w4[3 * i]);
			backward4(8, z, w);
		}
		backward4_0(_n / 4, _z);

		norm(_z, _n);
	}

	void mul(const Res & x)
	{
		Complex r_z[_n];
		for (size_t k = 0; k < _n; ++k) r_z[k] = x._z[k];

		const Complex * const w123 = mainPool.w123;

		forward4_0(_n / 4, _z);
		forward4_0(_n / 4, r_z);
		for (size_t j = 0; j < 4; ++j)
		{
			Complex * const z = &_z[j * _n / 4];
			Complex * const zr = &r_z[j * _n / 4];

			const Complex * const w = &w123[3 * (4 + j)];
			const Complex * const w4 = &w123[3 * 4 * (4 + j)];
			const Complex * const ws = &w123[3 * 2 * 4 * (4 + j)];
			forward4(8, z, w);
			forward4(8, zr, w);
			for (size_t i = 0; i < 4; ++i) forward4(2, &z[8 * i], &w4[3 * i]);
			for (size_t i = 0; i < 4; ++i) forward4(2, &zr[8 * i], &w4[3 * i]);
			for (size_t i = 0; i < 4; ++i) mul2(&z[8 * i], &zr[8 * i], ws[3 * 2 * i]);
			for (size_t i = 0; i < 4; ++i) backward4(2, &z[8 * i], &w4[3 * i]);
			backward4(8, z, w);
		}
		backward4_0(_n / 4, _z);

		norm(_z, _n);
	}

	void mul(const Res & x, const Res & y)
	{
		set(x); mul(y);
	}
};

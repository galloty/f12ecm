/*
Copyright 2021, Yves Gallot

f12ecm is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include "pool.h"
#include "transform.h"

#include <gmpxx.h>

inline mpz_class mpz(const uint64_t n)
{
	mpz_class z;
	mp_limb_t * const p_limb = mpz_limbs_write(z.get_mpz_t(), 1);
	p_limb[0] = n;
	mpz_limbs_finish(z.get_mpz_t(), 1);
	return z;
}

template<typename VComplex>
class Res : public Transform<VComplex>
{
	using base = Transform<VComplex>;

private:
	static const size_t _n = 256 / 2;	// F_12 = 65536^256 + 1
	VComplex * _z;

public:
	Res() : _z(nullptr) {}	// init(thread_index) must be called
	Res(const size_t thread_index) { init(thread_index); }
	Res(const Res &) = delete;
	Res & operator=(const Res &) = delete;

	void init(const size_t thread_index)
	{
		_z = static_cast<VComplex *>(mainPool.alloc(thread_index, _n * sizeof(VComplex)));
	}

	mpz_class get_z(const size_t i) const
	{
		const VComplex * const z = _z;

		mpz_class r = 0;
		for (size_t k = 0; k < _n; ++k)
		{
			const int h = int(z[_n - 1 - k].imag(i));
			r = (r << 16) + h;
		}
		for (size_t k = 0; k < _n; ++k)
		{
			const int l = int(z[_n - 1 - k].real(i));
			r = (r << 16) + l;
		}
		return r;
	}

	void set_z(const size_t i, const mpz_class & r)
	{
		const size_t r_size = mpz_size(r.get_mpz_t());
		// if ((r < 0) ||(r_size * sizeof(uint64_t) > 2 * _n * sizeof(uint16_t)))
		// {
		// 	throw std::runtime_error("set_z");
		// }

		uint16_t s[2 * _n];
		const mp_limb_t * r_ptr = mpz_limbs_read(r.get_mpz_t());
		for (size_t k = 0; k < r_size; ++k)
		{
			const uint64_t r_k = r_ptr[k];
			for (size_t j = 0; j < 4; ++j) s[4 * k + j] = uint16_t(r_k >> (16 * j));
		}
		for (size_t k = 4 * r_size; k < 2 * _n; ++k) s[k] = 0;

		VComplex * const z = _z;
		for (size_t k = 0; k < _n; ++k)
		{
			z[k].real(i) = double(s[k]);
			z[k].imag(i) = double(s[k + _n]);
		}
	}

	void set1()
	{
		VComplex * const z = _z;
		z[0] = VComplex(1.0, 0.0);
		for (size_t k = 1; k < _n; ++k) z[k] = VComplex(0.0, 0.0);
	}

	void set1(const size_t i)
	{
		VComplex * const z = _z;
		z[0].real(i) = double(1.0);
		z[0].imag(i) = double(0.0);
		for (size_t k = 1; k < _n; ++k)
		{
			z[k].real(i) = 0.0;
			z[k].imag(i) = 0.0;
		}
	}

	void set(const Res & x)
	{
		const VComplex * const zr = x._z;
		VComplex * const z = _z;
		for (size_t k = 0; k < _n; ++k) z[k] = zr[k];
	}

	void swap(Res & x)
	{
		std::swap(_z, x._z);
	}

	void add(const Res & x, const Res & y, const double m = 1)
	{
		const VComplex * const zx = x._z;
		const VComplex * const zy = y._z;
		VComplex * const z = _z;
		for (size_t k = 0; k < _n; ++k) z[k] = zx[k] + zy[k] * m;
	}

	void sub(const Res & x, const Res & y, const double m = 1)
	{
		const VComplex * const zx = x._z;
		const VComplex * const zy = y._z;
		VComplex * const z = _z;
		for (size_t k = 0; k < _n; ++k) z[k] = zx[k] - zy[k] * m;
	}

	static void addsub(Res & x, Res & y, const double m = 1)
	{
		VComplex * const zx = x._z;
		VComplex * const zy = y._z;
		for (size_t k = 0; k < _n; ++k) { const VComplex u = zx[k], v = zy[k]; zx[k] = u + v * m; zy[k] = u - v * m; }
	}

private:
	finline static void forward4_8(VComplex * const z, const Complex * const w123)
	{
		for (size_t j = 0; j < 4; ++j) base::forward4(8, &z[j * _n / 4], &w123[12 - 11 + 3 * j]);
	}

	finline static void backward4_8(VComplex * const z, const Complex * const w123)
	{
		for (size_t j = 0; j < 4; ++j) base::backward4(8, &z[j * _n / 4], &w123[12 - 11 + 3 * j]);
	}

	finline static void forward4_2(VComplex * const z, const Complex * const w123)
	{
		for (size_t j = 0; j < 4; ++j)
		{
			for (size_t i = 0; i < 4; ++i) base::forward4(2, &z[j * _n / 4 + 8 * i], &w123[48 - 11 + 3 * (4 * j + i)]);
		}
	}

	finline static void backward4_2(VComplex * const z, const Complex * const w123)
	{
		for (size_t j = 0; j < 4; ++j)
		{
			for (size_t i = 0; i < 4; ++i) base::backward4(2, &z[j * _n / 4 + 8 * i], &w123[48 - 11 + 3 * (4 * j + i)]);
		}
	}

	finline static void square2_4(VComplex * const z, const Complex * const w123)
	{
		for (size_t j = 0; j < 4; ++j)
		{
			for (size_t i = 0; i < 4; ++i) base::square2(&z[j * _n / 4 + 8 * i], &w123[96 - 11 + 6 * (4 * j + i)]);
		}
	}

	finline static void mul2_4(VComplex * const z, const VComplex * const zr, const Complex * const w123)
	{
		for (size_t j = 0; j < 4; ++j)
		{
			for (size_t i = 0; i < 4; ++i) base::mul2(&z[j * _n / 4 + 8 * i], &zr[j * _n / 4 + 8 * i], &w123[96 - 11 + 6 * (4 * j + i)]);
		}
	}

	finline static void _to_m(VComplex * const z)
	{
		const Complex * const w123 = mainPool.w123;

		base::forward4_0(_n / 4, z, w123);
		forward4_8(z, w123);
		forward4_2(z, w123);
	}

	finline static void _sqr(VComplex * const z)
	{
		const Complex * const w123 = mainPool.w123;

		base::forward4_0(_n / 4, z, w123);
		forward4_8(z, w123);
		forward4_2(z, w123);
		square2_4(z, w123);
		backward4_2(z, w123);
		backward4_8(z, w123);
		base::backward4_0(_n / 4, z, w123);
	}

	finline static void _mul(VComplex * const z, VComplex * const zr)
	{
		const Complex * const w123 = mainPool.w123;

		base::forward4_0(_n / 4, z, w123);
		base::forward4_0(_n / 4, zr, w123);
		forward4_8(z, w123);
		forward4_8(zr, w123);
		forward4_2(z, w123);
		forward4_2(zr, w123);
		mul2_4(z, zr, w123);
		backward4_2(z, w123);
		backward4_8(z, w123);
		base::backward4_0(_n / 4, z, w123);
	}

	finline static void _mul_m(VComplex * const z, const VComplex * const zr)
	{
		const Complex * const w123 = mainPool.w123;

		base::forward4_0(_n / 4, z, w123);
		forward4_8(z, w123);
		forward4_2(z, w123);
		mul2_4(z, zr, w123);
		backward4_2(z, w123);
		backward4_8(z, w123);
		base::backward4_0(_n / 4, z, w123);
	}

	finline static void _mul_t(VComplex * const z, VComplex * const zr)
	{
		const Complex * const w123 = mainPool.w123;

		base::forward4_0(_n / 4, zr, w123);
		forward4_8(zr, w123);
		forward4_2(zr, w123);
		mul2_4(z, zr, w123);
		backward4_2(z, w123);
		backward4_8(z, w123);
		base::backward4_0(_n / 4, z, w123);
	}

	finline static void _mul_mm(VComplex * const z, const VComplex * const zr)
	{
		const Complex * const w123 = mainPool.w123;

		mul2_4(z, zr, w123);
		backward4_2(z, w123);
		backward4_8(z, w123);
		base::backward4_0(_n / 4, z, w123);
	}

public:
	// this = T(this)
	void to_m()
	{
		_to_m(_z);
	}

	// this = this^2
	void sqr()
	{
		VComplex * const z = _z;
		_sqr(z);
		base::norm(z, _n, 1);
	}

	// this = 2 * this^2
	void sqr_add()
	{
		VComplex * const z = _z;
		_sqr(z);
		base::norm(z, _n, 2);
	}

	// this *= x, x = T(x)
	void mul(Res & x)
	{
		VComplex * const z = _z;
		VComplex * const zr = x._z;
		_mul(z, zr);
		base::norm(z, _n, 1);
	}

	// x is T(x), this *= x
	void mul_m(const Res & x)
	{
		VComplex * const z = _z;
		const VComplex * const zr = x._z;
		_mul_m(z, zr);
		base::norm(z, _n, 1);
	}

	// this is T(this), this *= x, x = T(x)
	void mul_t(const Res & x)
	{
		VComplex * const z = _z;
		VComplex * const zr = x._z;
		_mul_t(z, zr);
		base::norm(z, _n, 1);
	}

	// this is T(this), x is T(x), this *= x
	void mul_mm(const Res & x)
	{
		VComplex * const z = _z;
		const VComplex * const zr = x._z;
		_mul_mm(z, zr);
		base::norm(z, _n, 1);
	}

	// this is T(this), x is T(x), this *= 2 * x
	void mul_mm_add(const Res & x)
	{
		VComplex * const z = _z;
		const VComplex * const zr = x._z;
		_mul_mm(z, zr);
		base::norm(z, _n, 2);
	}
};

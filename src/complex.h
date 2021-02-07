/*
Copyright 2021, Yves Gallot

f12ecm is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>
#include <cmath>

#include <immintrin.h>

class Complex
{
private:
	double _re, _im;

public:
	Complex() {}
	explicit Complex(const double real, const double imag) : _re(real), _im(imag) {}

	double real() const { return _re; }
	double imag() const { return _im; }

	static Complex exp_2iPi(const size_t a, const size_t b)
	{
		const long double C2PI = 6.2831853071795864769252867665590057684L;
		const long double alpha = C2PI * (long double)(a) / (long double)(b);
		const double cs = double(cosl(alpha)), sn = double(sinl(alpha));
		return Complex(cs, sn / cs);
	}
};

#if defined (__AVX512F__)

class Complex8
{
};

typedef Complex8 VComplex;
#define VCOMPLEX_SIZE	8

#elif defined (__AVX__)

class Complex4
{
};

typedef Complex4 VComplex;
#define VCOMPLEX_SIZE	4

#else	// SSE4

class Complex2
{
private:
	__v2df _re, _im;

public:
	Complex2() {}
	explicit Complex2(const __v2df re, const __v2df im) : _re(re), _im(im) {}
	explicit Complex2(const double real, const double imag) : _re(_mm_set1_pd(real)), _im(_mm_set1_pd(imag)) {}
	Complex2(const Complex2 & rhs) : _re(rhs._re), _im(rhs._im) {}
	Complex2 & operator=(const Complex2 & rhs) { _re = rhs._re; _im = rhs._im; return *this; }

	const __v2df real() const { return _re; }
	const __v2df imag() const { return _im; }

	const double & real(const size_t i) const { return _re[i]; }
	const double & imag(const size_t i) const { return _im[i]; }

	double & real(const size_t i) { return _re[i]; }
	double & imag(const size_t i) { return _im[i]; }

	bool is_zero() const { return (_mm_movemask_pd(_mm_or_pd(_mm_cmpneq_pd(_re, _mm_setzero_pd()), _mm_cmpneq_pd(_im, _mm_setzero_pd()))) == 0); }

	Complex2 & operator+=(const Complex2 & rhs) { _re += rhs._re; _im += rhs._im; return *this; }
	Complex2 & operator-=(const Complex2 & rhs) { _re -= rhs._re; _im -= rhs._im; return *this; }
	Complex2 & operator*=(const double & f) { _re *= f; _im *= f; return *this; }

	Complex2 operator+(const Complex2 & rhs) const { return Complex2(_re + rhs._re, _im + rhs._im); }
	Complex2 operator-(const Complex2 & rhs) const { return Complex2(_re - rhs._re, _im - rhs._im); }
	Complex2 addi(const Complex2 & rhs) const { return Complex2(_re - rhs._im, _im + rhs._re); }
	Complex2 subi(const Complex2 & rhs) const { return Complex2(_re + rhs._im, _im - rhs._re); }

	Complex2 operator*(const Complex2 & rhs) const { return Complex2(_re * rhs._re - _im * rhs._im, _im * rhs._re + _re * rhs._im); }
	Complex2 operator*(const double & f) const { return Complex2(_re * f, _im * f); }

	Complex2 sqr() const { return Complex2(_re * _re - _im * _im, (_re + _re) * _im); }

	Complex2 mulW(const Complex & rhs) const { return Complex2((_re - _im * rhs.imag()) * rhs.real(), (_im + _re * rhs.imag()) * rhs.real()); }
	Complex2 mulWconj(const Complex & rhs) const { return Complex2((_re + _im * rhs.imag()) * rhs.real(), (_im - _re * rhs.imag()) * rhs.real()); }

	Complex2 round() const { return Complex2(_mm_round_pd(_re, _MM_FROUND_TO_NEAREST_INT), _mm_round_pd(_im, _MM_FROUND_TO_NEAREST_INT)); }
	Complex2 abs() const { const __v2df t = _mm_set1_pd(-0.0); return Complex2(_mm_andnot_pd(t, _re), _mm_andnot_pd(t, _im)); }
	Complex2 max(const Complex2 & rhs) const { return Complex2(_mm_max_pd(_re, rhs._re), _mm_max_pd(_im, rhs._im)); }

	double max() const { const __v2df t = _mm_max_pd(_re, _im); return std::fmax(((const double*)&t)[0], ((const double*)&t)[1]); }

	Complex2 rotate() const { return Complex2(-_im, _re); }	// f x^n = -f
};

typedef Complex2 VComplex;
#define VCOMPLEX_SIZE	2

#endif

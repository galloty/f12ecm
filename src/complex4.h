/*
Copyright 2021, Yves Gallot

f12ecm is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>
#include <cmath>

#include <immintrin.h>

#include "complex.h"

class Complex4
{
private:
	__v4df _re, _im;

public:
	Complex4() {}
	explicit Complex4(const __v4df re, const __v4df im) : _re(re), _im(im) {}
	explicit Complex4(const double real, const double imag) : _re(_mm256_set1_pd(real)), _im(_mm256_set1_pd(imag)) {}
	Complex4(const Complex4 & rhs) : _re(rhs._re), _im(rhs._im) {}
	Complex4 & operator=(const Complex4 & rhs) { _re = rhs._re; _im = rhs._im; return *this; }

	const __v4df real() const { return _re; }
	const __v4df imag() const { return _im; }

	const double & real(const size_t i) const { return _re[i]; }
	const double & imag(const size_t i) const { return _im[i]; }

	double & real(const size_t i) { return _re[i]; }
	double & imag(const size_t i) { return _im[i]; }

	bool is_zero() const { return (_mm256_movemask_pd(_mm256_or_pd(_mm256_cmp_pd(_re, _mm256_setzero_pd(), _CMP_NEQ_OQ), _mm256_cmp_pd(_im, _mm256_setzero_pd(), _CMP_NEQ_OQ))) == 0); }

	Complex4 & operator+=(const Complex4 & rhs) { _re += rhs._re; _im += rhs._im; return *this; }
	Complex4 & operator-=(const Complex4 & rhs) { _re -= rhs._re; _im -= rhs._im; return *this; }
	Complex4 & operator*=(const double & f) { _re *= f; _im *= f; return *this; }

	Complex4 operator+(const Complex4 & rhs) const { return Complex4(_re + rhs._re, _im + rhs._im); }
	Complex4 operator-(const Complex4 & rhs) const { return Complex4(_re - rhs._re, _im - rhs._im); }
	Complex4 addi(const Complex4 & rhs) const { return Complex4(_re - rhs._im, _im + rhs._re); }
	Complex4 subi(const Complex4 & rhs) const { return Complex4(_re + rhs._im, _im - rhs._re); }

	Complex4 operator*(const Complex4 & rhs) const { return Complex4(_re * rhs._re - _im * rhs._im, _im * rhs._re + _re * rhs._im); }
	Complex4 operator*(const double & f) const { return Complex4(_re * f, _im * f); }

	Complex4 sqr() const { return Complex4(_re * _re - _im * _im, (_re + _re) * _im); }

	Complex4 mulW(const Complex & rhs) const { return Complex4((_re - _im * rhs.imag()) * rhs.real(), (_im + _re * rhs.imag()) * rhs.real()); }
	Complex4 mulWconj(const Complex & rhs) const { return Complex4((_re + _im * rhs.imag()) * rhs.real(), (_im - _re * rhs.imag()) * rhs.real()); }

	Complex4 round() const { return Complex4(_mm256_round_pd(_re, _MM_FROUND_TO_NEAREST_INT), _mm256_round_pd(_im, _MM_FROUND_TO_NEAREST_INT)); }
	Complex4 abs() const { const __v4df t = _mm256_set1_pd(-0.0); return Complex4(_mm256_andnot_pd(t, _re), _mm256_andnot_pd(t, _im)); }
	Complex4 max(const Complex4 & rhs) const { return Complex4(_mm256_max_pd(_re, rhs._re), _mm256_max_pd(_im, rhs._im)); }

	Complex4 rotate() const { return Complex4(-_im, _re); }	// f x^n = -f
};

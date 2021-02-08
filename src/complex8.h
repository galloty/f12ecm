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

class Complex8
{
private:
	__v8df _re, _im;

public:
	Complex8() {}
	explicit Complex8(const __v8df re, const __v8df im) : _re(re), _im(im) {}
	explicit Complex8(const double real, const double imag) : _re(_mm512_set1_pd(real)), _im(_mm512_set1_pd(imag)) {}
	Complex8(const Complex8 & rhs) : _re(rhs._re), _im(rhs._im) {}
	Complex8 & operator=(const Complex8 & rhs) { _re = rhs._re; _im = rhs._im; return *this; }

	const __v8df real() const { return _re; }
	const __v8df imag() const { return _im; }

	const double & real(const size_t i) const { return _re[i]; }
	const double & imag(const size_t i) const { return _im[i]; }

	double & real(const size_t i) { return _re[i]; }
	double & imag(const size_t i) { return _im[i]; }

	bool is_zero() const
	{
		const __mmask8 mr = _mm512_cmp_pd_mask(_re, _mm512_setzero_pd(), _CMP_NEQ_OQ);
		const __mmask8 mi = _mm512_cmp_pd_mask(_im, _mm512_setzero_pd(), _CMP_NEQ_OQ);
		return (mr | mi) == 0;
	}

	Complex8 & operator+=(const Complex8 & rhs) { _re += rhs._re; _im += rhs._im; return *this; }
	Complex8 & operator-=(const Complex8 & rhs) { _re -= rhs._re; _im -= rhs._im; return *this; }
	Complex8 & operator*=(const double & f) { _re *= f; _im *= f; return *this; }

	Complex8 operator+(const Complex8 & rhs) const { return Complex8(_re + rhs._re, _im + rhs._im); }
	Complex8 operator-(const Complex8 & rhs) const { return Complex8(_re - rhs._re, _im - rhs._im); }
	Complex8 addi(const Complex8 & rhs) const { return Complex8(_re - rhs._im, _im + rhs._re); }
	Complex8 subi(const Complex8 & rhs) const { return Complex8(_re + rhs._im, _im - rhs._re); }

	Complex8 operator*(const Complex8 & rhs) const { return Complex8(_re * rhs._re - _im * rhs._im, _im * rhs._re + _re * rhs._im); }
	Complex8 operator*(const double & f) const { return Complex8(_re * f, _im * f); }

	Complex8 sqr() const { return Complex8(_re * _re - _im * _im, (_re + _re) * _im); }

	Complex8 mulW(const Complex & rhs) const { return Complex8((_re - _im * rhs.imag()) * rhs.real(), (_im + _re * rhs.imag()) * rhs.real()); }
	Complex8 mulWconj(const Complex & rhs) const { return Complex8((_re + _im * rhs.imag()) * rhs.real(), (_im - _re * rhs.imag()) * rhs.real()); }

	Complex8 round() const
	{
		return Complex8(_mm512_roundscale_pd(_re, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC),
		                _mm512_roundscale_pd(_im, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC));
	}
	Complex8 abs() const { return Complex8(_mm512_abs_pd(_re), _mm512_abs_pd(_im)); }
	Complex8 max(const Complex8 & rhs) const { return Complex8(_mm512_max_pd(_re, rhs._re), _mm512_max_pd(_im, rhs._im)); }

	Complex8 rotate() const { return Complex8(-_im, _re); }	// f x^n = -f
};

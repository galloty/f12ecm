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

class Complex2
{
private:
	__v2df _re, _im;

public:
	finline Complex2() {}
	finline explicit Complex2(const __v2df re, const __v2df im) : _re(re), _im(im) {}
	finline explicit Complex2(const double real, const double imag) : _re(_mm_set1_pd(real)), _im(_mm_set1_pd(imag)) {}
	finline Complex2(const Complex2 & rhs) : _re(rhs._re), _im(rhs._im) {}
	finline Complex2 & operator=(const Complex2 & rhs) { _re = rhs._re; _im = rhs._im; return *this; }

	finline const __v2df real() const { return _re; }
	finline const __v2df imag() const { return _im; }

	finline const double & real(const size_t i) const { return _re[i]; }
	finline const double & imag(const size_t i) const { return _im[i]; }

	finline double & real(const size_t i) { return _re[i]; }
	finline double & imag(const size_t i) { return _im[i]; }

	finline bool is_zero() const { return (_mm_movemask_pd(_mm_or_pd(_mm_cmpneq_pd(_re, _mm_setzero_pd()), _mm_cmpneq_pd(_im, _mm_setzero_pd()))) == 0); }

	finline Complex2 & operator+=(const Complex2 & rhs) { _re += rhs._re; _im += rhs._im; return *this; }
	finline Complex2 & operator-=(const Complex2 & rhs) { _re -= rhs._re; _im -= rhs._im; return *this; }
	finline Complex2 & operator*=(const double & f) { _re *= f; _im *= f; return *this; }

	finline Complex2 operator+(const Complex2 & rhs) const { return Complex2(_re + rhs._re, _im + rhs._im); }
	finline Complex2 operator-(const Complex2 & rhs) const { return Complex2(_re - rhs._re, _im - rhs._im); }
	finline Complex2 addi(const Complex2 & rhs) const { return Complex2(_re - rhs._im, _im + rhs._re); }
	finline Complex2 subi(const Complex2 & rhs) const { return Complex2(_re + rhs._im, _im - rhs._re); }

	finline Complex2 operator*(const Complex2 & rhs) const { return Complex2(_re * rhs._re - _im * rhs._im, _im * rhs._re + _re * rhs._im); }
	finline Complex2 operator*(const double & f) const { return Complex2(_re * f, _im * f); }

	finline Complex2 sqr() const { return Complex2(_re * _re - _im * _im, (_re + _re) * _im); }

	finline Complex2 mulW(const Complex & rhs) const { return Complex2((_re - _im * rhs.imag()) * rhs.real(), (_im + _re * rhs.imag()) * rhs.real()); }
	finline Complex2 mulWconj(const Complex & rhs) const { return Complex2((_re + _im * rhs.imag()) * rhs.real(), (_im - _re * rhs.imag()) * rhs.real()); }

	finline Complex2 round() const
	{
		return Complex2(_mm_round_pd(_re, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC),
		                _mm_round_pd(_im, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC));
	}
	finline Complex2 abs() const { const __v2df t = _mm_set1_pd(-0.0); return Complex2(_mm_andnot_pd(t, _re), _mm_andnot_pd(t, _im)); }
	finline Complex2 max(const Complex2 & rhs) const { return Complex2(_mm_max_pd(_re, rhs._re), _mm_max_pd(_im, rhs._im)); }

	finline double max() const { const __v2df t = _mm_max_pd(_re, _im); return std::fmax(((const double*)&t)[0], ((const double*)&t)[1]); }

	finline Complex2 rotate() const { return Complex2(-_im, _re); }	// f x^n = -f
};

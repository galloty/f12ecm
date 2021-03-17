/*
Copyright 2021, Yves Gallot

f12ecm is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>
#include <cmath>

#define finline		__attribute__((always_inline))

class Complex
{
private:
	double _re, _im;

public:
	finline Complex() {}
	finline explicit Complex(const double real, const double imag) : _re(real), _im(imag) {}

	finline double real() const { return _re; }
	finline double imag() const { return _im; }

	static Complex exp_2iPi(const size_t a, const size_t b)
	{
		const long double C2PI = 6.2831853071795864769252867665590057684L;
		const long double alpha = C2PI * (long double)(a) / (long double)(b);
		const double cs = double(cosl(alpha)), sn = double(sinl(alpha));
		return Complex(cs, sn / cs);
	}
};

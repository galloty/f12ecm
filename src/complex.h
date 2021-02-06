/*
Copyright 2021, Yves Gallot

f12ecm is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>
#include <cmath>

class Complex
{
private:
	double re, im;

public:
	Complex() {}
	explicit Complex(double real, double imag) : re(real), im(imag) {}
	Complex(const Complex & rhs) : re(rhs.re), im(rhs.im) {}
	Complex & operator=(const Complex & rhs) { re = rhs.re; im = rhs.im; return *this; }

	double real() const { return re; }
	double imag() const { return im; }

	bool is_zero() const { return ((re == 0.0) & (im == 0.0)); }

	Complex & operator+=(const Complex & rhs) { re += rhs.re; im += rhs.im; return *this; }
	Complex & operator-=(const Complex & rhs) { re -= rhs.re; im -= rhs.im; return *this; }
	Complex & operator*=(const double & f) { re *= f; im *= f; return *this; }

	Complex operator+(const Complex & rhs) const { return Complex(re + rhs.re, im + rhs.im); }
	Complex operator-(const Complex & rhs) const { return Complex(re - rhs.re, im - rhs.im); }
	Complex addi(const Complex & rhs) const { return Complex(re - rhs.im, im + rhs.re); }
	Complex subi(const Complex & rhs) const { return Complex(re + rhs.im, im - rhs.re); }

	Complex operator*(const Complex & rhs) const { return Complex(re * rhs.re - im * rhs.im, im * rhs.re + re * rhs.im); }
	Complex operator*(const double & f) const { return Complex(re * f, im * f); }

	Complex sqr() const { return Complex(re * re - im * im, (re + re) * im); }

	Complex mulW(const Complex & rhs) const { return Complex((re - im * rhs.im) * rhs.re, (im + re * rhs.im) * rhs.re); }
	Complex mulWconj(const Complex & rhs) const { return Complex((re + im * rhs.im) * rhs.re, (im - re * rhs.im) * rhs.re); }

	Complex round() const { return Complex(::round(re), ::round(im)); }

	Complex max_abs(const Complex & rhs) const { return Complex(::fmax(re, ::fabs(rhs.re)), ::fmax(im, ::fabs(rhs.im))); }

	Complex rotate() const { return Complex(-im, re); }	// f x^n = -f

	static Complex exp_2iPi(const size_t a, const size_t b)
	{
		const long double C2PI = 6.2831853071795864769252867665590057684L;
		const long double alpha = C2PI * (long double)(a) / (long double)(b);
		const double cs = double(cosl(alpha)), sn = double(sinl(alpha));
		return Complex(cs, sn / cs);
	}
};

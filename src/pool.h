/*
Copyright 2021, Yves Gallot

f12ecm is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include "complex.h"
#include <vector>

inline size_t bit_rev(const size_t i, const size_t n)
{
	size_t r = 0;
	for (size_t k = n, j = i; k != 1; k /= 2, j /= 2) r = (2 * r) | (j % 2);
	return r;
}

class MainPool
{
public:
	Complex w123[256];

private:
	std::vector<char *> _mem;
	std::vector<size_t> _offset;

public:
	void init(const size_t D, const size_t thread_count)
	{
		for (size_t s = 1; s <= 128 / 4; s *= 2)
		{
			Complex * const w123_s = &w123[3 * s];
			for (size_t i = 0; i < s; ++i)
			{
				const size_t r = bit_rev(i, 4 * s) + 1;
				w123_s[3 * i + 0] = Complex::exp_2iPi(r, 8 * s);
				w123_s[3 * i + 1] = Complex::exp_2iPi(r, 2 * 8 * s);
				w123_s[3 * i + 2] = Complex::exp_2iPi(3 * r, 2 * 8 * s);
			}
		}

		for (size_t i = 0; i < thread_count; ++i)
		{
			// ec: 3 points and 2 res: 8
			// P, T, C: 3 points: 6
			// S, beta: D points and D res: 3*D
			// g, alpha, t1, t2: 4 res: 4
			char * const ptr = static_cast<char *>(::_aligned_malloc(sizeof(Complex) * (256 / 2) * (3 * D + 18), 1024));
			_mem.push_back(ptr);
			_offset.push_back(0);
		}
	}

	void * alloc(const size_t thread_index, const size_t n)
	{
		const size_t offset = _offset[thread_index];
		_offset[thread_index] = offset + n;
		return _mem[thread_index] + offset;
	}

	void free(const size_t thread_index)
	{
		_offset[thread_index] = 0;
	}
};

static MainPool mainPool;

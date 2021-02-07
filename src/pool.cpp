/*
Copyright 2021, Yves Gallot

f12ecm is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include "pool.h"

inline size_t bit_rev(const size_t i, const size_t n)
{
	size_t r = 0;
	for (size_t k = n, j = i; k != 1; k /= 2, j /= 2) r = (2 * r) | (j % 2);
	return r;
}

void MainPool::init(const size_t D, const size_t vec_size, const size_t thread_count)
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
		char * const ptr = static_cast<char *>(::_aligned_malloc(vec_size * (256 / 2) * (3 * D + 18), 1024));
		_mem.push_back(ptr);
		_offset.push_back(0);
	}
}

void MainPool::release()
{
	for (char * const ptr : _mem) ::_aligned_free(ptr);
	_mem.clear();
	_offset.clear();
}

void * MainPool::alloc(const size_t thread_index, const size_t n)
{
	const size_t offset = _offset[thread_index];
	_offset[thread_index] = offset + n;
	return _mem[thread_index] + offset;
}

void MainPool::free(const size_t thread_index)
{
	_offset[thread_index] = 0;
}

MainPool mainPool;
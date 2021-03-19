/*
Copyright 2021, Yves Gallot

f12ecm is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include <stdexcept>

#include "pool.h"

inline size_t bit_rev(const size_t i, const size_t n)
{
	size_t r = 0;
	for (size_t k = n, j = i; k != 1; k /= 2, j /= 2) r = (2 * r) | (j % 2);
	return r;
}

void MainPool::init(const size_t D, const size_t W, const size_t vec_size, const size_t thread_count, size_t & size1, size_t & size2)
{
	w123[0] = Complex(0.92387953251128675612818318939678828682, 0.41421356237309504880168872420969807857);

	for (size_t s = 4; s <= 128 / 4; s *= 2)
	{
		Complex * const w123_s = &w123[3 * s - 11];
		for (size_t i = 0; i < s; ++i)
		{
			const size_t r = bit_rev(i, 4 * s) + 1;
			w123_s[3 * i + 0] = Complex::exp_2iPi(r, 8 * s);
			w123_s[3 * i + 1] = Complex::exp_2iPi(r, 2 * 8 * s);
			w123_s[3 * i + 2] = Complex::exp_2iPi(3 * r, 2 * 8 * s);
		}
	}

	// stage 1
	// ec_e:
	//  _d, _D, _E, _F: 4 res: 4
	//  _Pi: W/4 points: 3*W/4
	// ecm:
	//  _e:
	//  _P: 1 point: 3

	// stage 2
	// ec_m:
	//  _A2_4, _t: 2 res: 2
	//  _A, _B, _C, _T1, _T2, _T, _Tm: 7 points: 14
	// ecm:
	//  _m:
	//  P, Se, T, R, Rm: 5 points: 10
	//  S: D points: 2*D
	//  g, t1, t2: 3 res: 3

	size1 = vec_size * (256 / 2) * (3 * (W / 4) + 7);
	size2 = vec_size * (256 / 2) * (2 * D + 29);
	_size = size1 + size2;

	for (size_t i = 0; i < thread_count; ++i)
	{
		char * const ptr = static_cast<char *>(::_aligned_malloc(_size, 1024));
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
	if (offset + n > _size) throw std::runtime_error("MainPool::alloc");
	_offset[thread_index] = offset + n;
	return _mem[thread_index] + offset;
}

void MainPool::free(const size_t thread_index)
{
	_offset[thread_index] = 0;
}

MainPool mainPool;

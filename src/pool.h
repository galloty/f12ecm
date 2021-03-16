/*
Copyright 2021, Yves Gallot

f12ecm is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include "complex.h"
#include <vector>

class MainPool
{
public:
	Complex w123[256];

private:
	size_t _size;
	std::vector<char *> _mem;
	std::vector<size_t> _offset;

public:
	size_t init(const size_t D, const size_t W, const size_t vec_size, const size_t thread_count);
	void release();
	void * alloc(const size_t thread_index, const size_t n);
	void free(const size_t thread_index);
};

extern MainPool mainPool;
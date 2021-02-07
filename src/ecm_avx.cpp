/*
Copyright 2021, Yves Gallot

f12ecm is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include "ecm_i.h"

#include "Complex4.h"

#define VCOMPLEX_SIZE	4

#include "ecm.h"

void ECM_avx::run(const uint64_t B1, const uint64_t B2, const uint64_t sigma_0, const size_t thread_count)
{
	ECM<Complex4> & ecm = ECM<Complex4>::getInstance();
	ecm.run(B1, B2, sigma_0, thread_count);
}

void ECM_avx::quit()
{
	ECM<Complex4> & ecm = ECM<Complex4>::getInstance();
	ecm.quit();
}

/*
Copyright 2021, Yves Gallot

f12ecm is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include "ecm_i.h"
#include "ecm.h"
#include "Complex8.h"

void ECM_avx512::run(const uint64_t B1, const uint64_t B2, const uint64_t sigma_0, const bool isEdwards, const size_t thread_count, const std::string & ext)
{
	ECM<Complex8> & ecm = ECM<Complex8>::getInstance();
	ecm.run(B1, B2, sigma_0, isEdwards, thread_count, ext);
}

void ECM_avx512::quit()
{
	ECM<Complex8> & ecm = ECM<Complex8>::getInstance();
	ecm.quit();
}

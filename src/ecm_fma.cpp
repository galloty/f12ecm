/*
Copyright 2021, Yves Gallot

f12ecm is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include "ecm_i.h"
#include "ecm.h"
#define Complex4	Complex4_FMA
#include "Complex4.h"

void ECM_fma::run(const uint64_t B1, const uint64_t B2, const uint64_t sigma_0, const size_t thread_count, const std::string & ext)
{
	ECM<Complex4_FMA> & ecm = ECM<Complex4_FMA>::getInstance();
	ecm.run(B1, B2, sigma_0, thread_count, ext);
}

void ECM_fma::quit()
{
	ECM<Complex4_FMA> & ecm = ECM<Complex4_FMA>::getInstance();
	ecm.quit();
}

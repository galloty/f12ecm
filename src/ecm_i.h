/*
Copyright 2021, Yves Gallot

f12ecm is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>
#include <string>

struct ECM_i
{
	ECM_i() {}
	virtual ~ECM_i() {}

	virtual void run(const uint64_t B1, const uint64_t B2, const uint64_t sigma_0, const bool isEdwards, const size_t thread_count, const std::string & ext) = 0;
	virtual void quit() = 0;
};

struct ECM_sse4 : public ECM_i
{
	void run(const uint64_t B1, const uint64_t B2, const uint64_t sigma_0, const bool isEdwards, const size_t thread_count, const std::string & ext) override;
	void quit() override;
};

struct ECM_avx : public ECM_i
{
	void run(const uint64_t B1, const uint64_t B2, const uint64_t sigma_0, const bool isEdwards, const size_t thread_count, const std::string & ext) override;
	void quit() override;
};

struct ECM_fma : public ECM_i
{
	void run(const uint64_t B1, const uint64_t B2, const uint64_t sigma_0, const bool isEdwards, const size_t thread_count, const std::string & ext) override;
	void quit() override;
};

struct ECM_avx512 : public ECM_i
{
	void run(const uint64_t B1, const uint64_t B2, const uint64_t sigma_0, const bool isEdwards, const size_t thread_count, const std::string & ext) override;
	void quit() override;
};

/*
Copyright 2021, Yves Gallot

f12ecm is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include <cstdint>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <random>
#include <iostream>
#include <sstream>
#include <thread>

#if defined (_WIN32)
#include <Windows.h>
#else
#include <signal.h>
#endif

#include "ecm_i.h"

class application
{
private:
	struct deleter { void operator()(const application * const p) { delete p; } };

private:
	ECM_i * _ecm = nullptr;

private:
	static void quit(int) { getInstance().quit(); }
	void quit() { if (_ecm != nullptr) _ecm->quit(); }

private:
#if defined (_WIN32)
	static BOOL WINAPI HandlerRoutine(DWORD)
	{
		quit(1);
		return TRUE;
	}
#endif

public:
	application()
	{
#if defined (_WIN32)	
		SetConsoleCtrlHandler(HandlerRoutine, TRUE);
#else
		signal(SIGTERM, quit);
		signal(SIGINT, quit);
#endif
	}

	virtual ~application() {}

	static application & getInstance()
	{
		static std::unique_ptr<application, deleter> pInstance(new application());
		return *pInstance;
	}

private:
	static std::string header()
	{
		const char * const sys =
#if defined(_WIN64)
			"win64";
#elif defined(_WIN32)
			"win32";
#elif defined(__linux__)
#ifdef __x86_64
			"linux64";
#else
			"linux32";
#endif
#elif defined(__APPLE__)
			"macOS";
#else
			"unknown";
#endif

		std::ostringstream ssc;
#if defined(__GNUC__)
		ssc << "gcc-" << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__;
#elif defined(__clang__)
		ssc << "clang-" << __clang_major__ << "." << __clang_minor__ << "." << __clang_patchlevel__;
#endif

		std::ostringstream ss;
		ss << "f12ecm 0.1.0 " << sys << " " << ssc.str() << std::endl;
		ss << "Copyright (c) 2021, Yves Gallot" << std::endl;
		ss << "f12ecm is free source code, under the MIT license." << std::endl;
		ss << std::endl;
		return ss.str();
	}

private:
	static std::string usage()
	{
		std::ostringstream ss;
		ss << " Usage: f12ecm [options]  options may be specified in any order" << std::endl;
		ss << "   -t <n>      number of threads (default: 3)" << std::endl;
		ss << "   -b <n>      bound of stage 1 (default: B1 = 1,000,000)" << std::endl;
		ss << "   -B <n>      bound of stage 2 (default: B2 = 100*B1)" << std::endl;
		ss << "   -m <n>      Edwards curve: (s, t) = m * (-2, 4) in the Montgomery construction (default: random)" << std::endl;
		ss << "   -s <n>      Montgomery curve: sigma in Suyama construction" << std::endl;
		// ss << "   -c <n>      checking code using B1 = 10,000,000 and sigma = 16414" << std::endl;
		ss << std::endl;
		return ss.str();
	}

public:
	void run(int argc, char * argv[])
	{
		std::vector<std::string> args;
		for (int i = 1; i < argc; ++i) args.push_back(argv[i]);

		std::cout << header();

		if (args.empty()) std::cout << usage();

		std::random_device rd; std::uniform_int_distribution<uint64_t> dist(6, uint64_t(-1));

		// size_t thread_count = 1;	//3;
		// uint64_t B1 = 31*31, B2 = B1;	// 10000000
		// // uint64_t B1 = 9733667, B2 = 800000000;	// 10000000
		// uint64_t sigma_0 = 5553475;	//dist(rd);	// 16414
		// bool isEdwards = true;

		size_t thread_count = 3;
		uint64_t B1 = 1000000, B2 = 0;	// 10000000
		uint64_t sigma_0 = dist(rd);	// 16414
		bool isEdwards = false;

		// parse args
		for (size_t i = 0, size = args.size(); i < size; ++i)
		{
			const std::string & arg = args[i];

			if (arg.substr(0, 2) == "-t")
			{
				const std::string t = ((arg == "-t") && (i + 1 < size)) ? args[++i] : arg.substr(2);
				thread_count = std::atoi(t.c_str());
			}
			else if (arg.substr(0, 2) == "-b")
			{
				const std::string b = ((arg == "-b") && (i + 1 < size)) ? args[++i] : arg.substr(2);
				B1 = std::max(std::atoll(b.c_str()), 1000ll);
			}
			else if (arg.substr(0, 2) == "-B")
			{
				const std::string B = ((arg == "-B") && (i + 1 < size)) ? args[++i] : arg.substr(2);
				B2 = std::max(std::atoll(B.c_str()), 1000ll);
			}
			else if (arg.substr(0, 2) == "-m")
			{
				const std::string s = ((arg == "-m") && (i + 1 < size)) ? args[++i] : arg.substr(2);
				sigma_0 = std::max(std::atoll(s.c_str()), 2ll);
				isEdwards = true;
			}
			else if (arg.substr(0, 2) == "-s")
			{
				const std::string s = ((arg == "-s") && (i + 1 < size)) ? args[++i] : arg.substr(2);
				sigma_0 = std::max(std::atoll(s.c_str()), 6ll);
				isEdwards = false;
			}
		}

		if (B2 == 0) B2 = 100 * B1;
		if (B2 < B1) B2 = B1;

		__builtin_cpu_init();
		std::string ext;
		if (__builtin_cpu_supports("avx512f"))
		{
			_ecm = new ECM_avx512();
			ext = "avx-512";
		}
		else if (__builtin_cpu_supports("fma"))
		{
			_ecm = new ECM_fma();
			ext = "avx-fma";
		}
		else if (__builtin_cpu_supports("avx"))
		{
			_ecm = new ECM_avx();
			ext = "avx";
		}
		else
		{
			_ecm = new ECM_sse4();
			ext = "sse4";
		}
		_ecm->run(B1, B2, sigma_0, isEdwards, thread_count, ext);
		delete _ecm;
	}
};

int main(int argc, char * argv[])
{
	try
	{
		application & app = application::getInstance();
		app.run(argc, argv);
	}
	catch (const std::runtime_error & e)
	{
		std::ostringstream ss; ss << std::endl << "error: " << e.what() << "." << std::endl;
		std::cerr << ss.str();
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

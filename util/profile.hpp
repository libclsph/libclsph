#ifndef _PROFILE_HPP_
#define _PROFILE_HPP_

#include <chrono>
#include <functional>
#include <iostream>
#include <iomanip>

#define profile(name, fun) if(profile_block p = profile_block(name, fun))

enum steps_enum
{
    STEP1,
    STEP2,
    SORT,
    MEMORY_TRANSFERS
};

std::map<steps_enum,std::string> step_names    = { { STEP1, "Sph step 1" }, { STEP2, "Sph step 2" }, { SORT, "Sort" }, { MEMORY_TRANSFERS, "Memory transfers" } };
std::map<steps_enum,long long> step_run_length = { { STEP1, 0 },            { STEP2, 0 },            { SORT, 0 },      { MEMORY_TRANSFERS, 0 } };

class profile_block {
public:
	profile_block(std::string name, std::function<void(std::string, long long)> f = COUT_LOG)
	: block_name(name),
	  func(f) {
		start = std::chrono::high_resolution_clock::now();
	}

	~profile_block() {
		auto elapsed = std::chrono::high_resolution_clock::now() - start;
		long long microseconds = 
			std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
		if(func) {
			func(block_name, microseconds);
		}
	}

	operator bool() const {
		return true;
	}

	static std::function<void(std::string, long long)> COUT_LOG;

private:
	std::chrono::high_resolution_clock::time_point start;
	std::string block_name;
	std::function<void(std::string, long long)> func;
};

std::function<void(std::string, long long)> profile_block::COUT_LOG = 
	[] (std::string block_name, long long elapsed) {
		std::cout << std::fixed << std::setprecision(3) << 
			"[PROFILE] " << block_name	<< 
			" completed in " << ((float)elapsed) / 1000.f << 
			" milliseconds" << std::endl;
	};

#endif
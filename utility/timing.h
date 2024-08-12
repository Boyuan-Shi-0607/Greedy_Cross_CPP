#ifndef GREEDY_CROSS_CPP_TIMING_H
#define GREEDY_CROSS_CPP_TIMING_H

#include <iostream>
#include <chrono>

auto start_timing() {
	auto start = std::chrono::high_resolution_clock::now();
	return start;
}

void end_timing(std::chrono::time_point<std::chrono::high_resolution_clock>& start, std::string info, std::string unit) {
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::nano> duration = end - start;
	if (unit == "ns") {
		std::cout << "time in nano-seconds for " << info << " is " << duration.count() << std::endl;
	} else if (unit == "ms") {
		std::cout << "time in milli-seconds for " << info << " is " << duration.count() / std::pow(10, 6) << std::endl;
	} else if (unit == "mus") {
		std::cout << "time in micro-seconds for " << info << " is " << duration.count() / std::pow(10, 3) << std::endl;
	} else if (unit == "s") {
		std::cout << "time in seconds for " << info << " is " << duration.count() / std::pow(10, 9) << std::endl;
	}
}

#endif //GREEDY_CROSS_CPP_TIMING_H

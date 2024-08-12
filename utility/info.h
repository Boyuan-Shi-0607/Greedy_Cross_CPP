#ifndef GREEDY_CROSS_CPP_INFO_H
#define GREEDY_CROSS_CPP_INFO_H

#include <iostream>
#include <vector>

template <typename T1>
void info_vec(const std::string& information, const T1& data) {
	std::cout << information << " ";
	for (const auto& element : data) {
		std::cout << element << " ";
	}
	std::cout << std::endl;
}

template <typename T1>
void info_scalar(const std::string& information, const T1& data) {
	std::cout << information << " " << data << std::endl;
}


#endif //GREEDY_CROSS_CPP_INFO_H

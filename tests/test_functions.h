#ifndef GREEDY_CROSS_CPP_TEST_FUNCTIONS_H
#define GREEDY_CROSS_CPP_TEST_FUNCTIONS_H

#include <iostream>
#include "../utility/change_vector_type.h"
#include "../include/parameters.h"
#include <cmath>
#include <cstdlib>
#include <ctime>

double func_0(const std::vector<int>& y) {
	// low rank
	std::vector<double> x = change_vector_type(y);
	return std::exp(- x[0] * std::cos(x[1] / dim[1] + x[2] / dim[2]) / (static_cast<double>(dim[0]) / 10));
}

double func_1(const std::vector<int>& x) {
	return matrix(x[0], x[1]); // rank: rank in parameters
}

double func_2(const std::vector<int>& x) {
	return std::cos(x[0] + x[1]); // rank 2
}

double func_3(const std::vector<int>& x) {
	return std::cos(x[0] + x[1] + x[2]); // rank (2, 2)
}

#endif //GREEDY_CROSS_CPP_TEST_FUNCTIONS_H

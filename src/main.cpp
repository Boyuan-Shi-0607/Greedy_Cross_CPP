#include <iostream>
#include "../include/greedy_cross.h"
#include "../utility/change_vector_type.h"
#include <cmath>

double func(std::vector<int> y) {
	std::vector<double> x = change_vector_type(y);
	return std::exp(x[0] * std::cos(x[1] / dim[1] + x[2] / dim[2]) / (static_cast<double>(dim[0]) / 10));
	// return std::cos(static_cast<double>(x[1]) / dim[1] + static_cast<double>(x[2]) / dim[2] + static_cast<double>(x[0]) / dim[0]);
}

int main () {
	dim = {1000, 50, 50};
	n_swp = 10;
	tol = {1e-7, 1e-7, 1e-7};
	verbose = true;
	cross_data results = greedy_cross(func);
	std::vector<double> eps = results.eps;
	return 0;
}

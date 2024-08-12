#ifndef GREEDY_CROSS_CPP_SWEEP_H
#define GREEDY_CROSS_CPP_SWEEP_H

#include <iostream>
#include "maxvol.h"
#include "parameters.h"
#include <random>
#include <tuple>

double sweep(std::vector<std::vector<std::vector<int>>>& I,
             std::vector<std::vector<std::vector<int>>>& J,
			 int& alpha,
			 const std::function<double(std::vector<int>)>& func,
			 const double& local_tol,
			 Eigen::MatrixXd& inverse_pivot_mat) {
	
	std::vector<std::vector<int>> I1;
	if (alpha == 0) {
		I1 = {{}};
	} else {
		I1 = I[alpha-1];
	}
	std::vector<std::vector<int>> J1;
	if (alpha == dim.size() - 2) {
		J1 = {{}};
	} else {
		J1 = J[alpha+1];
	}
	
	std::vector<std::vector<int>> pivots_x;
	std::vector<std::vector<int>> pivots_y;
	std::vector<double> res;
	double eps_full_rank = maxvol(func, alpha, I[alpha], J[alpha],
								  I1, J1, pivots_x,
								  pivots_y, res,
								  inverse_pivot_mat);
	
	// std::cout << "eps full rank " << eps_full_rank << std::endl;
	if (pivots_x.empty()) {
		return eps_full_rank;
	} else {
		auto maxElementIter = std::max_element(res.begin(), res.end());
		auto maxIndex = std::distance(res.begin(), maxElementIter);
		std::vector<int> new_pivot_x = pivots_x[maxIndex];
		std::vector<int> new_pivot_y = pivots_y[maxIndex];
		if (eps_full_rank > local_tol) {
			I[alpha].push_back(new_pivot_x);
			J[alpha].push_back(new_pivot_y);
		}
		return eps_full_rank;
	}
}

#endif //GREEDY_CROSS_CPP_SWEEP_H

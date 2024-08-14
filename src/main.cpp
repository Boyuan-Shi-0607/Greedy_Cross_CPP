#include <iostream>
#include "../include/greedy_cross.h"
#include "../utility/change_vector_type.h"
#include "../tests/test_functions.h"
#include <Eigen/Dense>
#include <cmath>
#include <cstdlib>
#include <ctime>

int main () {
	/*
	dim = {1000, 50, 50};
	n_swp = 10;
	tol = {1e-1, 1e-1, 1e-1};
	verbose = true;
	cross_data results = greedy_cross(func);
	// std::vector<double> eps = results.eps;
	 */
	
	std::random_device rd_mat;
	std::mt19937 gen_mat(rd_mat());
	std::uniform_real_distribution<double> dist_mat(-10., 10.);

	for (int r = 0; r < rank; ++r) {
		Eigen::VectorXd u(rows);
		Eigen::VectorXd v(cols);
		
		for (int i = 0; i < rows; ++i) {
			u(i) = dist_mat(gen_mat);
		}
		for (int j = 0; j < cols; ++j) {
			v(j) = dist_mat(gen_mat);
		}
		matrix = u * v.transpose();
	}
	
	dim = {20, 20};
	n_swp = 10;
	tol = {1e-10, 1e-10};
	verbose = true;
	cross_data results = greedy_cross(func_1);
	return 0;
	
}

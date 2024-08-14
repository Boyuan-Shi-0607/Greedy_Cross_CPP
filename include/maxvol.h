#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <boost/multi_array.hpp>
#include "../utility/vector_concatenation.h"
#include "parameters.h"
#include "../utility/info.h"

double maxvol(const std::function<double(std::vector<int>)>& func,
			const int& alpha,
			std::vector<std::vector<std::vector<int>>>& I,
			std::vector<std::vector<std::vector<int>>>& J,
			const double& local_tol) {
	
	if (I[alpha].size() != J[alpha].size()) {
		throw std::runtime_error("Ranks must be equal across one pivot");
	}
	
	std::vector<std::vector<int>> pivots_x;
	std::vector<std::vector<int>> pivots_y;
	std::vector<double> res;
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
	
	std::vector<double> res_conservative;
	Eigen::MatrixXd pivot_mat(I[alpha].size(), I[alpha].size());
	
	int row = 0;
	for (const std::vector<int>& i: I[alpha]) {
		int col = 0;
		for (const std::vector<int>& j: J[alpha]) {
			pivot_mat(row, col) = func(concatenateVectors(i, j));
			col ++;
		}
		row ++;
	}
	
	/*
	if (std::abs(pivot_mat.determinant()) < 1e-12) {
		throw std::runtime_error("Pivot matrix has too small determinant, <1e-12.");
	} else if (std::abs(pivot_mat.determinant()) > 1e12) {
		throw std::runtime_error("Pivot matrix has too large determinant, >1e12.");
	}
	*/
	
	Eigen::MatrixXd inverse_pivot_mat = pivot_mat.inverse();
	
	// Search for new pivots.
	// Can be done in alternative search.
	for (const std::vector<int>& x0_b: I1) {
		for (int x0_t=0; x0_t<dim[alpha]; ++x0_t) {
			for (int y0_h=0; y0_h<dim[alpha+1]; ++y0_h) {
				for (const std::vector<int>& y0_b: J1) {
					std::vector<int> x0 = x0_b;
					x0.push_back(x0_t);
					std::vector<int> y0 = y0_b;
					y0.insert(y0.begin(), y0_h);
					
					double func_val = func(concatenateVectors(x0, y0));
					// std::cout << func_val << std::endl;
					Eigen::VectorXd vec_1(J[alpha].size());
					Eigen::VectorXd vec_2(I[alpha].size());
					int index = 0;
					for (const std::vector<int>& j:J[alpha]) {
						vec_1(index) = func(concatenateVectors(x0, j));
						index ++;
					}
					index = 0;
					for (const std::vector<int>& i:I[alpha]) {
						vec_2(index) = func(concatenateVectors(i, y0));
						index ++;
					}
					// std::cout << "res " << std::abs(func_val - vec_1.dot(inverse_pivot_mat * vec_2)) << std::endl;
					res_conservative.push_back(std::abs(func_val - vec_1.dot(inverse_pivot_mat * vec_2)));
					
					Eigen::MatrixXd new_pivot_mat(I[alpha].size()+1, I[alpha].size()+1);
					for (int i=0; i<I[alpha].size(); ++i) {
						for (int j=0; j<I[alpha].size(); ++j) {
							new_pivot_mat(i, j) = pivot_mat(i, j);
						}
					}
					
					for (int i=0; i<I[alpha].size(); ++i) {
						new_pivot_mat(static_cast<int>(I[alpha].size()), i) = vec_1(i);
						new_pivot_mat(i, static_cast<int>(I[alpha].size())) = vec_2(i);
					}
					
					new_pivot_mat(static_cast<int>(I[alpha].size()), static_cast<int>(I[alpha].size())) = func_val;
					double fac = func_val - vec_1.dot(inverse_pivot_mat * vec_2);
					// Eigen::JacobiSVD<Eigen::MatrixXd> svd(new_pivot_mat);
					// svd.setThreshold(1e-14);
					// int rank = static_cast<int>(svd.nonzeroSingularValues());
					
					if (std::abs(fac) > 1e-14) {
						res.push_back(std::abs(func_val - vec_1.dot(inverse_pivot_mat * vec_2)));
						// Set of all the pivots that can be added,
						// i.e., those pivots that would not lead to
						// non-full-ranked pivot matrix.
						pivots_x.push_back(x0);
						pivots_y.push_back(y0);
					}
				}
			}
		}
	}
	double eps_max = *std::max_element(res_conservative.begin(), res_conservative.end());
	
	// If there is such pivot that when adding it, the new pivot matrix remains full-ranked,
	// then return the error estimate during the search.
	// (A more conservative error estimate)
	if (pivots_x.empty()) {
		return eps_max;
	// Otherwise, find the new pivot based on the maxvol principle among
	// legitimate candidates;
	} else {
		auto maxElementIter = std::max_element(res.begin(), res.end());
		auto maxIndex = std::distance(res.begin(), maxElementIter);
		std::vector<int> new_pivot_x = pivots_x[maxIndex];
		std::vector<int> new_pivot_y = pivots_y[maxIndex];
		// update I and J if error estimate is greater than tolerance
		if (eps_max > local_tol) {
			I[alpha].push_back(new_pivot_x);
			J[alpha].push_back(new_pivot_y);
		}
		return eps_max;
	}
}

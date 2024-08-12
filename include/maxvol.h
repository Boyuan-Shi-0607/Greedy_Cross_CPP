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
			const std::vector<std::vector<int>>& I_alpha,
			const std::vector<std::vector<int>>& J_alpha_p_1,
			const std::vector<std::vector<int>>& I_alpha_m_1,
			const std::vector<std::vector<int>>& J_alpha_p_2,
			std::vector<std::vector<int>>& pivots_x,
			std::vector<std::vector<int>>& pivots_y,
			std::vector<double>& res,
			Eigen::MatrixXd& inverse_pivot_mat) {
	
	if (I_alpha.size() != J_alpha_p_1.size()) {
		throw std::runtime_error("Ranks must be equal across one pivot");
	}
	
	std::vector<double> res_conservative;
	Eigen::MatrixXd pivot_mat(I_alpha.size(), I_alpha.size());
	
	int row = 0;
	for (const std::vector<int>& i: I_alpha) {
		int col = 0;
		for (const std::vector<int>& j: J_alpha_p_1) {
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
	
	// Eigen::MatrixXd inverse_pivot_mat = pivot_mat.inverse();
	
	for (const std::vector<int>& x0_b: I_alpha_m_1) {
		for (int x0_t=0; x0_t<dim[alpha]; ++x0_t) {
			for (int y0_h=0; y0_h<dim[alpha+1]; ++y0_h) {
				for (const std::vector<int>& y0_b: J_alpha_p_2) {
					std::vector<int> x0 = x0_b;
					x0.push_back(x0_t);
					std::vector<int> y0 = y0_b;
					y0.insert(y0.begin(), y0_h);
					
					double func_val = func(concatenateVectors(x0, y0));
					// std::cout << func_val << std::endl;
					Eigen::VectorXd vec_1(J_alpha_p_1.size());
					Eigen::VectorXd vec_2(I_alpha.size());
					int index = 0;
					for (const std::vector<int>& j:J_alpha_p_1) {
						vec_1(index) = func(concatenateVectors(x0, j));
						index ++;
					}
					index = 0;
					for (const std::vector<int>& i:I_alpha) {
						vec_2(index) = func(concatenateVectors(i, y0));
						index ++;
					}
					// std::cout << "res " << std::abs(func_val - vec_1.dot(inverse_pivot_mat * vec_2)) << std::endl;
					res_conservative.push_back(std::abs(func_val - vec_1.dot(inverse_pivot_mat * vec_2)));
					
					Eigen::MatrixXd new_pivot_mat(I_alpha.size()+1, I_alpha.size()+1);
					for (int i=0; i<I_alpha.size(); ++i) {
						for (int j=0; j<I_alpha.size(); ++j) {
							new_pivot_mat(i, j) = pivot_mat(i, j);
						}
					}
					for (int i=0; i<I_alpha.size(); ++i) {
						new_pivot_mat(static_cast<int>(I_alpha.size()), i) = vec_1(i);
						new_pivot_mat(i, static_cast<int>(I_alpha.size())) = vec_2(i);
					}
					
					new_pivot_mat(static_cast<int>(I_alpha.size()), static_cast<int>(I_alpha.size())) = func_val;
					
					Eigen::JacobiSVD<Eigen::MatrixXd> svd(new_pivot_mat);
					// svd.setThreshold(1e-14);
					int rank = static_cast<int>(svd.nonzeroSingularValues());
					// std::cout << rank - I_alpha.size() << std::endl;
					if (rank == I_alpha.size() + 1) {
						res.push_back(std::abs(func_val - vec_1.dot(inverse_pivot_mat * vec_2)));
						pivots_x.push_back(x0);
						pivots_y.push_back(y0);
						
						// update inverse pivot matrix.
						double beta = std::pow((func_val - vec_1.dot(inverse_pivot_mat * vec_2)), -1);
						Eigen::VectorXd inv_p_vec2 = inverse_pivot_mat * vec_2;
						Eigen::VectorXd vec1_T_inv_p = (inverse_pivot_mat.transpose() * vec_1).transpose();
						inverse_pivot_mat += beta * inv_p_vec2 * vec1_T_inv_p;
						inverse_pivot_mat.resize(static_cast<int>(I_alpha.size())+1, static_cast<int>(I_alpha.size())+1);
						inverse_pivot_mat.bottomLeftCorner(1, I_alpha.size()) = - beta * vec1_T_inv_p;
						inverse_pivot_mat.topRightCorner(I_alpha.size(), 1) = -beta * inv_p_vec2;
						inverse_pivot_mat(static_cast<int>(I_alpha.size()), static_cast<int>(I_alpha.size())) = beta;
					}
				}
			}
		}
	}
	// std::cout << "max element " << *std::max_element(res_conservative.begin(), res_conservative.end()) << std::endl;
	return *std::max_element(res_conservative.begin(), res_conservative.end());
}

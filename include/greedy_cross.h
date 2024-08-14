#ifndef GREEDY_CROSS_CPP_GREEDY_CROSS_H
#define GREEDY_CROSS_CPP_GREEDY_CROSS_H
#include <iostream>
#include "maxvol.h"
#include <random>
#include <tuple>
#include "../utility/info.h"
#include "../utility/timing.h"
#include "../utility/vector_concatenation.h"
#include "parameters.h"

struct cross_data {
	std::vector<std::vector<std::vector<int>>> I;
	std::vector<std::vector<std::vector<int>>> J;
	std::vector<double> eps;
};

cross_data greedy_cross(const std::function<double(std::vector<int>)>& func) {
	
	std::vector<double> eps(dim.size()-1);
	std::vector<int> y0(dim.size());
	for (int i=0; i<dim.size();++i) {
		std::uniform_int_distribution<int> dis(0, dim[i]-1);
		y0[i] = dis(gen);
	}
	std::vector<std::vector<std::vector<int>>> I;
	std::vector<std::vector<std::vector<int>>> J;
	
	// nested condition
	for (int alpha=0; alpha<dim.size()-1;++alpha) {
		std::vector<int> y0_alpha_left;
		y0_alpha_left.reserve(alpha+1);
		for (int k=0; k<alpha+1;++k) {
			y0_alpha_left.push_back(y0[k]);
		}
		I.push_back({y0_alpha_left});
	}
	
	for (int alpha=1; alpha<dim.size();++alpha) {
		std::vector<int> y0_alpha_right;
		y0_alpha_right.reserve(dim.size()-alpha-1);
		for (int k=alpha; k<dim.size();++k) {
			y0_alpha_right.push_back(y0[k]);
		}
		J.push_back({y0_alpha_right});
	}
	
	// ensure the initial rank-1 1 * 1 pivot matrix is not singular.
	for (int alpha=0; alpha<dim.size()-1;++alpha) {
		for (const std::vector<int> &i: I[alpha]) {
			for (const std::vector<int> &j: J[alpha]) {
				if (std::abs(func(concatenateVectors(i, j))) < std::pow(10, -12)) {
					throw std::runtime_error("Initial rank-1 choice leads to "
											 "singular 1 * 1 pivot matrices.");
				}
			}
		}
	}
	
	int i = 0;
	while (i < n_swp) {
		auto start = start_timing();
		
		int counter = 0;
		// forward sweep
		for (int alpha=0; alpha<dim.size()-1;++alpha) {
			eps[alpha] = maxvol(func, alpha, I, J, tol[alpha]);
		}
		
		// backward sweep
		for (int alpha=static_cast<int>(dim.size())-2; alpha>=0;--alpha) {
			eps[alpha] = maxvol(func, alpha, I, J, tol[alpha]);
			if (eps[alpha] < tol[alpha]) {
				counter += 1;
			}
		}
		if (verbose) {
			std::cout << "---------------------------------------" << std::endl;
			end_timing(start, "sweep " + std::to_string(i), "s");
			std::cout << "number of sweeps " << i+1 << std::endl;
			std::vector<int> remaining_node;
			std::vector<double> remaining_eps;
			std::vector<int> ranks;
			for (int alpha=0;alpha<dim.size()-1;++alpha) {
				if (eps[alpha] > tol[alpha]) {
					remaining_node.push_back(alpha);
					remaining_eps.push_back(eps[alpha]);
				}
			}
			ranks.reserve(I.size());
			for (const auto & item : I) {
				ranks.push_back(static_cast<int>(item.size()));
			}
			info_scalar("total number of nodes: ", static_cast<int>(dim.size())-1);
			info_vec("max_dx: ", eps);
			info_vec("ranks: ", ranks);
			if (not remaining_node.empty()) {
				info_scalar("number of remaining nodes: ", static_cast<int>(remaining_node.size()));
				info_vec("remaining nodes: ", remaining_node);
				info_vec("remaining eps: ", remaining_eps);
			}
		}
		
		if (counter == dim.size()-1) {
			std::cout << "Sweep finished. Tol achieved." << std::endl;
			break;
		}
		
		++i;
	}
	return cross_data{I, J, eps};
}

#endif //GREEDY_CROSS_CPP_GREEDY_CROSS_H

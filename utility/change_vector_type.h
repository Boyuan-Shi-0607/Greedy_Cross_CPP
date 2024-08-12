#ifndef GREEDY_CROSS_CPP_CHANGE_VECTOR_TYPE_H
#define GREEDY_CROSS_CPP_CHANGE_VECTOR_TYPE_H

#include <iostream>
#include <vector>

std::vector<double> change_vector_type(const std::vector<int>& x) {
	std::vector<double> converted_vec;
	converted_vec.reserve(x.size());
	for (int val:x) {
		converted_vec.push_back(static_cast<double>(val));
	}
	return converted_vec;
}


#endif //GREEDY_CROSS_CPP_CHANGE_VECTOR_TYPE_H

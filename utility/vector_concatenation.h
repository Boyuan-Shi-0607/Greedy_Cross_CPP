#ifndef GREEDY_CROSS_CPP_VECTOR_CONCATENATION_H
#define GREEDY_CROSS_CPP_VECTOR_CONCATENATION_H

template <typename T1>
T1 concatenateVectors(const T1& v1, const T1& v2) {
	T1 result = v1;
	result.insert(result.end(), v2.begin(), v2.end());
	return result;
}

#endif //GREEDY_CROSS_CPP_VECTOR_CONCATENATION_H

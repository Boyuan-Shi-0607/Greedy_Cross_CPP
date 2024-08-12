#ifndef GREEDY_CROSS_CPP_PARAMETERS_H
#define GREEDY_CROSS_CPP_PARAMETERS_H

#include <iostream>
#include <vector>
#include <random>

std::vector<int> dim;
int n_swp;
std::vector<double> tol;
bool verbose = true;
std::random_device rd;
std::mt19937 gen(rd());

#endif //GREEDY_CROSS_CPP_PARAMETERS_H

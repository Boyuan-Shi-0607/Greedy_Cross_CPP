#ifndef GREEDY_CROSS_CPP_PARAMETERS_H
#define GREEDY_CROSS_CPP_PARAMETERS_H

#include <iostream>
#include <vector>
#include <random>
#include <Eigen/Dense>

std::vector<int> dim;
int n_swp;
std::vector<double> tol;
bool verbose = true;
std::random_device rd;
std::mt19937 gen(rd());

const int rows = 20;
const int cols = 20;
const int rank = 7;

Eigen::MatrixXd matrix(rows, cols);

#endif //GREEDY_CROSS_CPP_PARAMETERS_H

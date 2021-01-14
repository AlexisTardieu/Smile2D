#ifndef EIGEN_H
#define EIGEN_H

#include <Eigen/Sparse>

typedef Eigen::SparseMatrix<double, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;
typedef Eigen::Triplet<double> T;

#endif // EIGEN_H

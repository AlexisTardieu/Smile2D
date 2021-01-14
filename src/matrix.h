#ifndef MATRIX_H
#define MATRIX_H

#include "data.h"
#include "eigen.h"
#include <vector>

struct Mesh;
struct Point;

void buildMatrixElasticity (Matrix * mat, Mesh * mesh);

void buildSecMember (Point (*f) (double, double), Vector & secMember, Mesh * mesh);

void imposeDirichlet (Matrix * mat, Vector & secMember, Mesh * mesh,
                      Point (*g) (double, double), std::vector<int> & list);

#endif // MATRIX_H

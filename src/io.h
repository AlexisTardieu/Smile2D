#ifndef IO_H
#define IO_H

#include "data.h"
#include "eigen.h"
#include <string>

struct Mesh;

void writeForPlot (std::string filename, Mesh * mesh, Vector & sol);

void readDataFile (std::string filename, Mesh * mesh);

void printMatrixBlock (Matrix * mat, int minR, int maxR, int minC, int maxC);

void printMatrix (Matrix * mat);

#endif // IO_H

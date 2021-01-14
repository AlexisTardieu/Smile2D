#ifndef BORDER_H
#define BORDER_H

#include <vector>

struct Mesh;
struct Point;

// Les conditions aux limites possibles

Point zeroFunction (double x, double y);
Point oneFunction (double x, double y);
Point gravityFunction (double x, double y);
Point trinomFunction (double x, double y);

// Les différents ensembles de points où imposer des conditions aux limites

std::vector<int> borderLeft (Mesh * mesh);
std::vector<int> borderRight (Mesh * mesh);
std::vector<int> borderBottom (Mesh * mesh);
std::vector<int> borderTop (Mesh * mesh);
std::vector<int> borderOmega (Mesh * mesh);
std::vector<int> outsidePoints (Mesh * mesh);

#endif // BORDER_H

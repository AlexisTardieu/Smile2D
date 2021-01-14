#ifndef SOLVER_H
#define SOLVER_H

#include "eigen.h"

void solveCG (Matrix * mat, Vector & sol, Vector & secMember);

#endif // SOLVER_H

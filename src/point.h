#ifndef POINT_H
#define POINT_H

#include <vector>
#include "data.h"

struct Point
{
    // Coordonnées
    double x;
    double y;

    // Localisation : INSIDE (phi < 0) ou OUTSIDE (phi >= 0)
    Loc loc;

    // Voisins dans l'ordre : [Gauche, Droit, Bas, Haut]
    std::vector<Point*> neighbours;
};

double Distance (Point * pP, Point * pQ);
double prodScal (Point * pA, Point * pB);

void searchMin (Point * pP, std::vector<Point *> & gamma, int * rang, double * dist);

#endif // POINT_H

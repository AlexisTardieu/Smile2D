#include "point.h"

#include <iostream>
#include <cmath>

// Distance euclidienne entre 2 points

double Distance (Point * pP, Point * pQ)
{
    return std::sqrt ((pP->x - pQ->x) * (pP->x - pQ->x) + (pP->y - pQ->y) * (pP->y - pQ->y));
}


// Produit scalaire euclidien entre 2 vecteurs (vus comme des points car 2 coordonnées)

double prodScal (Point * pA, Point * pB)
{
    return pA->x * pB->x + pA->y * pB->y;
}


/* Cherche le point pQ de la liste gamma qui est le + proche du point pP,
ainsi que la distance associée (sans signe) */

void searchMin (Point * pP, std::vector<Point *> & gamma, int * rang, double * dist)
{
    int taille = gamma.size ();

    for (int m = 0; m < taille; m++) {

        Point * pQ = gamma [m];
        double d = Distance (pP, pQ);

        // std::cout << pP->x << " " << pP->y << '\n';

        if (d < *dist) {
            *dist = d;
            *rang = m;
        }
    }

    return;
}

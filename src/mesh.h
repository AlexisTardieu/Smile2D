#ifndef MESH_H
#define MESH_H

#include <vector>
#include "data.h"

struct Mesh {
    // Nombres de points de discrétisation
    int Nx, Ny;

    // Dimensions du grand domaine
    double xmin, xmax;
    double ymin, ymax;
    double Lx, Ly;

    // Pas de discrétisation
    double hx, hy;

    // Forme choisie pour la frontière irrégulière immergée
    int geom;

    // Coefficients de Lamé
    double mu, lambda;
    double eps; // Tolérance pour le lissage de la fonction de Heaviside au bord

    // Informations sur les points de grille
    std::vector<double> phiOnGrid; // Valeurs de la level-set phi (domaine irrégulier)
    std::vector<double> HeavisideOnGrid;
    std::vector<double> muOnGrid; // Valeurs de mu via la fonction de Heaviside
    std::vector<double> lambdaOnGrid; // Valeurs de lambda via la fonction de Heaviside
    std::vector<Loc> loc; // Point intérieur : 1 (phi < 0), extérieur = 0 (phi >= 0)
};

void makeVectorPhi (Mesh * mesh);
void makeVectorHeaviside (Mesh * mesh);
void makeVectorMuLambda (Mesh * mesh);
#endif // MESH_H

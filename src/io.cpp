#include "io.h"

#include <iostream>
#include <fstream>
#include "mesh.h"
#include "point.h"


// Fichier de sortie
// fichier : # x y u v [COULEUR]

/* Ecriture dans un fichier de sortie, au format : #X #Y #U #V #COLOR, où
[X,Y] sont les coordonnées d'un point, [U,V] le déplacement solution, et
COLOR la couleur dans laquelle on souhaite afficher le point avec gnuplot */

void writeForPlot (std::string filename, Mesh * mesh, Vector & sol)
{
    std::cout << "STEP 5 : Writing in the file : " << filename << "..." << std::endl;

    std::ofstream fichier (filename);

    // Vérification que l'ouverture a fonctionné
    if (! fichier)
    {
        std::cout << "Impossible to open the file : " << filename << "." << std::endl;

        return;
    }

    fichier << SPC "#X" << SPC "#Y" << SPC "#U" << SPC "#V"
            << SPC "#COLOR" << SPC "#MU" << SPC "#LAMBDA" << SPC "#PHI" << std::endl;

    int Nx = mesh->Nx;
    int Ny = mesh->Ny;
    int Ntot = Nx * Ny;

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {

            // Indice du point (i,j) dans la numérotation globale
            int idx = j * mesh->Nx + i;

            // Coordonnées [x,y] du point
            double x = mesh->xmin + static_cast<double> (i) * mesh->hx;
            double y = mesh->ymin + static_cast<double> (j) * mesh->hy;

            // Le déplacement [u,v] qu'il subit, calculé par résolution de AX = B
            double u = sol [idx];
            double v = sol [idx + Ntot];

            // La couleur pour affichage Gnuplot, codée par un entier
            int colour = mesh->loc [idx];

            // Les coefficients de Lamé et la valeur de la level-set phi
            double mu = mesh->muOnGrid [idx];
            double lambda = mesh->lambdaOnGrid [idx];
            double phi = mesh->phiOnGrid [idx];

            fichier << SPC x << SPC y << SPC u << SPC v << SPC colour << SPC mu << SPC lambda << SPC phi << std::endl;
            // fichier << x << "\t" << y << "\t" << u << "\t" << v << std::endl;
        }
    }

    std::cout << "STEP 5 : done !" << std::endl << std::endl;
    return;
}

// Lecture du fichier de données INPUT.dat

void readDataFile (std::string filename, Mesh * mesh)
{
    std::cout << "STEP 1 : Reading the file : " << filename << "..." << std::endl;

    std::ifstream fichier (filename);

    // Vérification que l'ouverture a fonctionné
    if (! fichier)
    {
        std::cout << "Impossible to open the file : " << filename << "." << std::endl;

        return;
    }

    std::string mot = "", name = "", star = "";

    while (! fichier.eof ())
    {
        fichier >> mot >> std::ws;

        if (mot == "GRILLE")
        {
            fichier >> name >> star >> mesh->Nx >> std::ws;
            fichier >> name >> star >> mesh->Ny >> std::ws;
            fichier >> name >> star >> mesh->xmin >> std::ws;
            fichier >> name >> star >> mesh->xmax >> std::ws;
            fichier >> name >> star >> mesh->ymin >> std::ws;
            fichier >> name >> star >> mesh->ymax >> std::ws;
        }

        else if (mot == "CARACTERISTIQUES")
        {
            fichier >> name >> star >> mesh->geom >> std::ws;
            fichier >> name >> star >> mesh->mu >> std::ws;
            fichier >> name >> star >> mesh->lambda >> std::ws;
        }

        else
            break;
    }

    mesh->Lx = mesh->xmax - mesh->xmin;
    mesh->Ly = mesh->ymax - mesh->ymin;

    mesh->hx = mesh->Lx / static_cast<double> (mesh->Nx - 1);
    mesh->hy = mesh->Ly / static_cast<double> (mesh->Ny - 1);

    mesh->eps = 0. * std::max (mesh->hx, mesh->hy);

    makeVectorPhi (mesh);
    makeVectorMuLambda (mesh);

    for (int i = 0; i < mesh->Nx; i++) {
        for (int j = 0; j < mesh->Ny; j++) {

            int idx = j * mesh->Nx + i;
            double phi = mesh->phiOnGrid [idx];

            if (phi < 0)
                mesh->loc [idx] = INSIDE;

            else
                mesh->loc [idx] = OUTSIDE;
        }
    }

    // std::cout << "hx = " << mesh->hx << " et hy = " << mesh->hy << std::endl;
    // std::cout << "mu = " << mesh->mu << " et lambda = " << mesh->lambda << std::endl;

    std::cout << "STEP 1 : done !" << std::endl << std::endl;

    return;
}


/* Affichage d'un bloc d'une matrice : entre les lignes minR et maxR (inclues),
et les colonnes minC et maxC (inclues), en commençant la numérotation à 0 */

void printMatrixBlock (Matrix * mat, int minR, int maxR, int minC, int maxC)
{
    int nbRows = mat->rows ();
    int nbCols = mat->cols ();

    if (minR < 0 || maxR >= nbRows || minC < 0 || maxC >= nbCols)
    {
        std::cout << "Slicing out of the matrix size." << std::endl;
        return;
    }

    std::cout << std::endl;

    for (int i = minR; i <= maxR; i++) {

        std::cout << "[ ";

        for (int j = minC; j < maxC; j++)
            std::cout << SPCl mat->coeffRef (i, j) << ", ";

        std::cout << mat->coeffRef (i, maxC) << " ]" << std::endl;
    }

    std::cout << std::endl;

    return;
}

// Affichage des 4 blocs de la matrice A = [ [A_11, A_12], [A_21, A_22]]

void printMatrix (Matrix * mat)
{
    // La matrice est carrée et de taille 2*Nx*Ny en pratique
    int Ntot = mat->rows ();
    int taille = Ntot / 2;

    std::cout << "Block (1,1) of the matrix :" << std::endl;
    printMatrixBlock (mat, 0, taille - 1, 0, taille - 1);

    std::cout << "Block (1,2) of the matrix :" << std::endl;
    printMatrixBlock (mat, 0, taille - 1, taille, Ntot - 1);

    std::cout << "Block (2,1) of the matrix :" << std::endl;
    printMatrixBlock (mat, taille, Ntot - 1, 0, taille - 1);

    std::cout << "Block (2,2) of the matrix :" << std::endl;
    printMatrixBlock (mat, taille, Ntot - 1, taille, Ntot - 1);

    return;
}

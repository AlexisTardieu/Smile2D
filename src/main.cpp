#include <iostream>
#include <fstream>

#include "mesh.h"
#include "point.h"
#include "matrix.h"
#include "levelset.h"
#include "io.h"
#include "border.h"
#include "eigen.h"
#include "solver.h"

int main (int argc, char **argv)
{
    std::cout << "***********************************************" << std::endl;
    std::cout << "**    Bienvenue dans Smile2D ! Un solveur    **" << std::endl;
    std::cout << "**  par différences finies des équations de  **" << std::endl;
    std::cout << "**  l'élasticité linéaire en 2 dimensions :  **" << std::endl;
    std::cout << "**           - div (sigma (u)) = f           **" << std::endl;
    std::cout << "***********************************************" << std::endl << std::endl;

    // std::cout << "Test des fonctions d'affichage avec la matrice identité (taille paire !) :" << std::endl << std::endl;
    // Matrix * mat = new Matrix (10, 10);
    // mat->setIdentity ();
    // printMatrix (mat);
    // delete mat;

    // Le maillage
    Mesh mesh;

    // Le fichier de données d'entrée
    std::string filename = "../INPUT.dat";

    // On lit le fichier de données et on renseigne les variables du maillage
    readDataFile (filename, & mesh);

    // Nombre d'inconnues : Nx*Ny pour u + Nx*Ny pour v
    int taille = 2 * mesh.Nx * mesh.Ny;

    // La matrice de discrétisation du problème de l'élasticité linéaire 2D
    Matrix mat (taille, taille);
    buildMatrixElasticity (& mat, & mesh);

    // std::cout << "Avant Dirichlet :" << std::endl;
    // printMatrix (& mat);

    // Le vecteur X = [U,V] représentant le déplacement solution sur la grille
    Vector sol (taille);

    // Le second membre B = -F = - [f1,f2]
    Vector secMember (taille);
    buildSecMember (gravityFunction, secMember, & mesh);

    // std::cout << secMember.transpose () << std::endl;

    // Liste des points où imposer les conditions aux limites
    // std::vector<int> listOfPoints = borderOmega (& mesh);

    // std::vector<int> outside = outsidePoints (& mesh);
    // imposeDirichlet (& mat, secMember, & mesh, zeroFunction, outside);

    std::vector<int> idx;
    for (int i = 1; i < mesh.Nx-1; ++i)
        for (int j = 1; j < mesh.Ny-1; ++j)
        {
            int id = j * mesh.Nx + i;

            if (mesh.loc [id + 1] == OUTSIDE && mesh.loc [id - 1] == INSIDE)
                idx.push_back(id);
            else if (mesh.loc [id - 1] == OUTSIDE && mesh.loc [id + 1] == INSIDE)
                idx.push_back(id);
            else if (mesh.loc [id - mesh.Nx] == OUTSIDE && mesh.loc [id + mesh.Nx] == INSIDE)
                idx.push_back(id);
            else if (mesh.loc [id + mesh.Nx] == OUTSIDE && mesh.loc [id - mesh.Nx] == INSIDE)
                idx.push_back(id);
        }

    std::vector<int> MYCUSTOMPOINTS;
    for (int i = 0; i < idx.size () / 1.0; ++i)
        MYCUSTOMPOINTS.push_back(idx [i]);

    // imposeDirichlet (& mat, secMember, & mesh, oneFunction, MYCUSTOMPOINTS);

    // imposeDirichlet (& mat, secMember, & mesh, oneFunction, idx);

    // imposeDirichlet (& mat, secMember, & mesh, zeroFunction, MYCUSTOMPOINTS);



    std::vector<int> top = borderTop (& mesh);
    std::vector<int> bottom = borderBottom (& mesh);
    std::vector<int> left = borderLeft (& mesh);
    std::vector<int> right = borderRight (& mesh);

    // imposeDirichlet (& mat, secMember, & mesh, zeroFunction, top);
    imposeDirichlet (& mat, secMember, & mesh, zeroFunction, left);
    // imposeDirichlet (& mat, secMember, & mesh, zeroFunction, bottom);
    // imposeDirichlet (& mat, secMember, & mesh, zeroFunction, right);

    // std::vector<int> outsidepoints = outsidePoints (& mesh);
    // imposeDirichlet (& mat, secMember, & mesh, zeroFunction, outsidepoints);


    // Conditions de Dirichlet
    // imposeDirichlet (& mat, secMember, & mesh, oneFunction, listOfPoints);


    // std::cout << "Après Dirichlet :" << std::endl;
    // printMatrix (& mat);
    // std::cout << secMember.transpose () << std::endl;

    // Résolution du problème matriciel AX = B avec le gradient conjugué
    solveCG (& mat, sol, secMember);

    // for (int i = 0; i < mesh.Nx * mesh.Ny; ++i)
    // {
    //     sol.coeffRef(i) = mesh.muOnGrid [i];
    //     sol.coeffRef(i + mesh.Nx * mesh.Ny) = mesh.muOnGrid [i];
    //
    //     // sol.coeffRef(i) = mesh.phiOnGrid [i];
    // }
    /* Le fichier de sortie : #X #Y #U #V #COLOUR
    où [X,Y] = coordonnées du point, et [U,V] = déplacement solution */
    std::string output = "solution.dat";

    // Ecriture de la solution dans un fichier
    writeForPlot (output, & mesh, sol);

    return 0;
}

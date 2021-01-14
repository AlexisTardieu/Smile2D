#include "matrix.h"

#include <iostream>
#include <cmath>
#include "mesh.h"
#include "point.h"

// Construit la matrice A de discrétisation DF du système de l'élasticité linéaire

void buildMatrixElasticity (Matrix * mat, Mesh * mesh)
{
    std::cout << "STEP 2 : Building the matrix A for the linear elasticity problem..." << std::endl;

    // Nombres de points de discrétisation et pas d'espace

    int Nx = mesh->Nx;
    int Ny = mesh->Ny;
    // int Ntot = 2 * Nx * Ny;

    double hx = mesh->hx;
    double hy = mesh->hy;

    std::vector<T> listOfTriplets;
    listOfTriplets.reserve (10 * 2 * Nx * Ny);

    int Ntot = Nx * Ny;

    auto idx_u = [Ntot, Nx] (int i, int j) -> int {return j * Nx + i;};
    auto idx_v = [Ntot, Nx] (int i, int j) -> int {return j * Nx + i + Ntot;};

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {

            // Indice du point (i,j) dans la numérotation globale
            int idx = j * Nx + i;

            // Coefficients de Lamé
            double mu = mesh->muOnGrid [idx];
            double lambda = mesh->lambdaOnGrid [idx];

            // Coefficients de la matrice
            double a = (2 * mu + lambda) / (hx * hx);
            double b = mu / (hy * hy);
            double c = mu / (hx * hx);
            double d = (2 * mu + lambda) / (hy * hy);
            double e = (lambda + mu) / (4 * hx * hy);

            double diagU = - (a + b);// + 1e0 * (1. - mesh->HeavisideOnGrid[idx]);
            double diagV = - (c + d) ;//+ 1e0 * (1. - mesh->HeavisideOnGrid[idx]);

            // Coefficient diagonal du 1er bloc = -(a+b)
            // listOfTriplets.push_back (T (idx, idx, - (a + b)));
            listOfTriplets.push_back (T (idx_u(i, j), idx_u(i, j), diagU));

            // Coefficient diagonal du 2ème bloc = -(c+d)
            // listOfTriplets.push_back (T (idx + Ntot, idx + Ntot, - (c + d)));
            listOfTriplets.push_back (T (idx_v(i, j), idx_v(i, j), diagV));


            if (j < Ny - 1)
            {
                // Coefficient décalé de Nx du 1er bloc = b
                // listOfTriplets.push_back (T (idx, idx + Nx, b));
                listOfTriplets.push_back (T (idx_u(i, j), idx_u(i, j+1), b));

                // Coefficient décalé de Nx du 2ème bloc = d
                // listOfTriplets.push_back (T (idx + Ntot, idx + Ntot + Nx, d));
                listOfTriplets.push_back (T (idx_v(i, j), idx_v(i, j+1), d));

            }

            if (i < Nx - 1)
            {
                // Coefficient décalé de 1 du 1er bloc = a
                // listOfTriplets.push_back (T (idx, idx + 1, a));
                listOfTriplets.push_back (T (idx_u(i, j), idx_u(i+1, j), a));

                // Coefficient décalé de 1 du 2ème bloc = c
                // listOfTriplets.push_back (T (idx + Ntot, idx + Ntot + 1, c));
                listOfTriplets.push_back (T (idx_v(i, j), idx_v(i+1, j), c));


                if (j < Ny - 1)
                {
                    // Coefficient décalé de Nx*Ny + Nx + 1 = e
                    // listOfTriplets.push_back (T (idx, idx + Ntot + Nx + 1, e));
                    listOfTriplets.push_back (T (idx_u(i, j), idx_v(i+1, j+1), e));

                }

                if (j > 0)
                {
                    // Coefficient décalé de Nx*Ny - Nx + 1 = -e
                    // listOfTriplets.push_back (T (idx, idx + Ntot - Nx + 1, - e));
                    listOfTriplets.push_back (T (idx_u(i, j), idx_v(i+1, j-1), -e));

                }
            }

            if (i > 0)
            {
                if (j< Ny - 1)
                {
                    // Coefficient décalé de Nx*Ny + Nx - 1 = -e
                    // listOfTriplets.push_back (T (idx, idx + Ntot + Nx - 1, - e));
                    listOfTriplets.push_back (T (idx_u(i, j), idx_v(i-1, j+1), -e));

                }

                if (j > 0)
                {
                    // Coefficient décalé de Nx*Ny - Nx - 1 = e
                    // listOfTriplets.push_back (T (idx, idx + Ntot - Nx - 1, e));
                    listOfTriplets.push_back (T (idx_u(i, j), idx_v(i-1, j-1), e));

                }
            }
        }
    }

    // for (unsigned int id = 0; id < listOfTriplets.size (); ++id)
    // {
    //     auto it = listOfTriplets [id];
    //     if (it.row () < 0 || it.row () >= 20000 || it.col () < 0 || it.col () >= 20000)
    //     std::cout << "triplet invalid = " << id << " row " << it.row() << " col " << it.col() << " value " << it.value() << std::endl;
    //
    //
    //     if (it.col () < it.row())
    //         std::cout << "!! triplet invalid = " << id << " row " << it.row() << " col " << it.col() << " value " << it.value() << std::endl;
    //
    // }

    mat->setFromTriplets (listOfTriplets.begin (), listOfTriplets.end ());

    /* La matrice mat est symétrique, donc on a construit sa partie
    triangulaire supérieure juste avant, et on lui ajoute sa transposée */
    Matrix * transp = new Matrix (mat->transpose ());
    * mat += * transp;
    delete transp;

    // On retire les éventuels zéros inutiles pour optimiser le stockage
    * mat = mat->pruned ();

    std::cout << "STEP 2 : done !" << std::endl << std::endl;

    return;
}


// Construit le vecteur de second membre initial : B = -F

void buildSecMember (Point (*f) (double, double), Vector & secMember, Mesh * mesh)
{
    std::cout << "STEP 3 : Building the second member B for the linear elasticity problem..." << std::endl;

    int Ntot = mesh->Nx * mesh->Ny;

    for (int i = 0; i < mesh->Nx; i++) {
        for (int j = 0; j < mesh->Ny; j++) {

            double x = mesh->xmin + static_cast<double> (i) * mesh->hx;
            double y = mesh->ymin + static_cast<double> (j) * mesh->hy;

            Point valF = f (x, y);

            double f1 = valF.x;
            double f2 = valF.y;

            int idx = j * mesh->Nx + i;

            secMember [idx] = - f1;
            secMember [idx + Ntot] = - f2;
        }
    }

    std::cout << "STEP 3 : done !" << std::endl << std::endl;

    return;
}


/* Impose des conditions de Dirichlet g(x,y) sur tous les points dont
les indices globaux sont dans list : la ligne de la matrice A est mise à
l'identité, et le 2nd membre B est mis à la valeur g(x,y) voulue */

void imposeDirichlet (Matrix * mat, Vector & secMember, Mesh * mesh,
                      Point (*g) (double, double), std::vector<int> & list)
{
    int Ntot = mesh->Nx * mesh->Ny;

    for (int k : list)
    {
        mat->row (k) *= 0.;
        mat->row (k + Ntot) *= 0.;
    }

    * mat = mat->transpose ();

    for (int k : list)
    {
        int i = k % mesh->Nx;
        int j = (k - i) / mesh->Nx;

        if (j == 0)
            std::cout << " j = " << j << "  and i = " << i << std::endl;
        // Coordonnées du point d'indice global k = j*Nx+i
        double x = mesh->xmin + static_cast<double> (i) * mesh->hx;
        double y = mesh->ymin + static_cast<double> (j) * mesh->hy;

        // std::cout << "k % Nx = " << x << " et k / Nx = " << y << std::endl;

        Point valG = g (x, y);

        double g1 = valG.x;
        double g2 = valG.y;

        Vector vec1 = mat->row (k).transpose ();
        Vector vec2 = mat->row (k + Ntot).transpose ();

        for (int i = 0; i < mesh->Nx; i++) {
            for (int j = 0; j < mesh->Ny; j++) {

                int idx = j * mesh->Nx + i;

                vec1 [idx] *= g1;
                vec1 [idx + Ntot] *= g2;

                vec2 [idx] *= g1;
                vec2 [idx + Ntot] *= g2;
            }
        }

        secMember -= vec1;
        secMember -= vec2;

        // secMember -= g (x, y) * mat->row (k).transpose ();

        mat->row (k) *= 0.;
        mat->row (k + Ntot) *= 0.;

        mat->coeffRef (k, k) = 1.;
        mat->coeffRef (k + Ntot, k + Ntot) = 1.;

        secMember.coeffRef (k) = g1;
        secMember.coeffRef (k + Ntot) = g2;
    }

    * mat = mat->transpose ();
    * mat = mat->pruned ();

    return;
}

#include "mesh.h"

#include <iostream>
#include <cmath>
#include "levelset.h"

// Construit le vecteur 1D des valeurs de la level-set phi sur la grille

void makeVectorPhi (Mesh * mesh)
{
    mesh->phiOnGrid.resize (mesh->Nx * mesh->Ny);
    mesh->HeavisideOnGrid.resize (mesh->Nx * mesh->Ny);
    mesh->muOnGrid.resize (mesh->Nx * mesh->Ny);
    mesh->lambdaOnGrid.resize (mesh->Nx * mesh->Ny);
    mesh->loc.resize (mesh->Nx * mesh->Ny);

    switch (mesh->geom)
    {
        case 0 :
            std::cout << "You have picked the whole domain." << std::endl;
            phi_border (mesh);
            break;

        case 1 :
            std::cout << "You have picked the circle." << std::endl;
            phi_circle (mesh);
            break;

        case 2 :
            std::cout << "You have picked the ellipse." << std::endl;
            phi_ellipse (mesh);
            break;

        case 3 :
            std::cout << "You have picked the square." << std::endl;
            phi_square (mesh);
            break;

        case 4 :
            std::cout << "You have picked the seastar." << std::endl;
            phi_seastar (mesh);
            break;

        case 5 :
            std::cout << "You have picked the face." << std::endl;
            phi_face (mesh);
            break;

        default :
            std::cout << "The form you picked doesn\'t exist." << std::endl;
            break;
    }

    return;
}


// Construit le vecteur 1D des valeurs de mu et de lambda sur la grille

void makeVectorMuLambda (Mesh * mesh)
{
    for (int i = 0; i < mesh->Nx; i++) {
        for (int j = 0; j < mesh->Ny; j++) {

            int idx = j * mesh->Nx + i;
            double phi_ij = mesh->phiOnGrid [idx];
            double eps = mesh->eps;

            /* Valeur de la fonction de Heaviside au point (i,j) :
             H_ij  =  1  si phi_ij <= -2eps
                      0.5*[1 - (1/eps)*(phi_ij + eps) - (1/pi)*sin((pi/eps)*(phi_ij + eps))]  sinon
                      0  si phi_ij >= 0 */
            double h_ij = 0.;

            // if (phi_ij <= - 2 * eps)
            //     h_ij = 1.;
            //
            // else if (phi_ij >= 0.)
            //     h_ij = 0.;
            //
            // else
            //     h_ij = 0.5* (1 - (1. / eps) * (phi_ij + eps) - (1. / PI) * std::sin ((PI / eps) * (phi_ij + eps)));
            //
            // if (phi_ij <= -  eps)
            //     h_ij = 1.;
            //
            // else if (phi_ij >= eps)
            //     h_ij = 0.;
            //
            // else
            //     h_ij = 0.5* (1 - (1. / eps) * phi_ij - (1. / PI) * std::sin ((PI / eps) * phi_ij));

            // if (phi_ij <= 0.)
            //     h_ij = 1.;
            //
            // else if (phi_ij >= 2. * eps)
            //     h_ij = 0.;
            //
            // else
            //     h_ij = 0.5* (1 - (1. / eps) * (phi_ij - eps) - (1. / PI) * std::sin ((PI / eps) * (phi_ij - eps)));

            h_ij = 1.;
            double mu = mesh->mu;
            double lambda = mesh->lambda;

            // mesh->muOnGrid [idx] = mu * h_ij;
            // mesh->lambdaOnGrid [idx] = lambda * h_ij;

            mesh->HeavisideOnGrid [idx] = h_ij;
            mesh->muOnGrid [idx] = mu * h_ij;
            mesh->lambdaOnGrid [idx] = lambda * h_ij;
        }
    }

    return;
}

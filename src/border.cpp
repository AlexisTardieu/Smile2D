#include "border.h"
#include "mesh.h"
#include "data.h"
#include "point.h"

// Fonction nulle

Point zeroFunction (double x, double y)
{
    Point p = {0., 0.};

    return p;
}

Point trinomFunction (double x, double y)
{
    // Point p = {0., (y == 0 ? +1. : -1.) * 0.5 * x * (1 - x)};

    return {0., - 0.5 * x * (1 - x)};
}

Point oneFunction (double x, double y)
{
    return {0.05, 0.0};
    // return {(x == 0) ? 0.1 : 0.0, (y == 1) ? -0.1 : 0.0};

    // Point p = {x * x + y * y, x * x + y * y};

    // return p;
}

Point gravityFunction (double x, double y)
{
    return {0., -260.};
}


// Liste des indices globaux des points du bord gauche

std::vector<int> borderLeft (Mesh * mesh)
{
    int Nx = mesh->Nx;
    int Ny = mesh->Ny;

    std::vector<int> list;
    list.reserve (Ny);

    for (int j = 0; j < Ny; j++)
        list.push_back (j * Nx);

    return list;
}


// Liste des indices globaux des points du bord droit

std::vector<int> borderRight (Mesh * mesh)
{
    int Nx = mesh->Nx;
    int Ny = mesh->Ny;

    std::vector<int> list;
    list.reserve (Ny);

    for (int j = 0; j < Ny; j++)
        list.push_back (j * Nx + (Nx - 1));

    return list;
}


// Liste des indices globaux des points du bord bas

std::vector<int> borderBottom (Mesh * mesh)
{
    int Nx = mesh->Nx;

    std::vector<int> list;
    list.reserve (Nx);

    // for (int i = 1; i < Nx - 1; i++)
    for (int i = 0; i < Nx; i++)
        list.push_back (i);

    return list;
}


// Liste des indices globaux des points du bord haut

std::vector<int> borderTop (Mesh * mesh)
{
    int Nx = mesh->Nx;
    int Ny = mesh->Ny;

    std::vector<int> list;
    list.reserve (Nx);

    // for (int i = 1; i < Nx - 1; i++)
    for (int i = 0; i < Nx; i++)
        list.push_back ((Ny - 1) * Nx + i);

    return list;
}


// Liste des indices globaux des points du bord du domaine rectangulaire Omega

std::vector<int> borderOmega (Mesh * mesh)
{
    int Nx = mesh->Nx;
    int Ny = mesh->Ny;

    int taille = 2 * (Nx + Ny) - 4;

    std::vector<int> list;
    list.reserve (taille);

    // j * Nx + i
    for (int i = 0; i < Nx; i++)
    {
        list.push_back (i); // Bord du bas
        list.push_back (Nx * (Ny - 1) + i); // Bord du haut
    }

    for (int j = 1; j < Ny - 1; j++)
    {
        list.push_back (j * Nx);
        list.push_back (j * Nx + (Nx - 1));
    }

    return list;
}


// Liste des indices globaux des points extérieurs au bord irrégulier

std::vector<int> outsidePoints (Mesh * mesh)
{
    int Nx = mesh->Nx;
    int Ny = mesh->Ny;

    int taille = Nx * Ny;

    std::vector<int> list;
    list.reserve (taille);

    std::vector<Loc> loc = mesh->loc;

    for (int i = 0; i < Nx; i ++) {
        for (int j = 0; j < Ny; j++) {

            int idx = j * Nx + i;

            if (loc [idx] == OUTSIDE)
                list.push_back (idx);
        }
    }

    return list;
}

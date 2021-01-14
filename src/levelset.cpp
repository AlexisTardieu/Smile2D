#include "levelset.h"

#include <iostream>
#include <cmath>
#include "point.h"
#include "mesh.h"
#include "data.h"

/* Bloc élémentaire : segment [A,B] de vecteur normal n unitaire
tel que det(AB,n) > 0 si s = 1, et det(AB,n) < 0 si s = -1 */

void segment (Point * pA, Point * pB, double s, std::vector<Point *> & gamma, std::vector<Point *> & normalVec, Mesh * mesh)
{
    int Np = 10 * floor (Distance (pA, pB) * (std::max (mesh->Nx, mesh->Ny) / std::min (mesh->Lx, mesh->Ly)));

    double norm = Distance (pA, pB);

    double absVec = s * (pA->y - pB->y) / norm;
    double ordVec = s * (pB->x - pA->x) / norm;

    for (int k = 0; k < Np; k++) {

        double coeff = static_cast<double> (k) / static_cast<double> (Np);

        double abs = pA->x + coeff * (pB->x - pA->x);
        double ord = pA->y + coeff * (pB->y - pA->y);

        Point * newPoint = new Point ({abs, ord});
        gamma.push_back (newPoint);

        Point * newVec = new Point ({absVec, ordVec});
        normalVec.push_back (newVec);
    }
}


/* Bloc élémentaire : arc de cercle de centre C et rayon r,
entre les angles th_init et th_end, de vecteur normal n unitaire
pointant vers l'extérieur si s = 1, et vers l'intérieur si s = -1 */

void arcCercle (Point * pC, double r, double th_init, double th_end, double s, std::vector<Point *> & gamma, std::vector<Point *> & normalVec, Mesh * mesh)
{
    int Np = 10 * floor ((th_end - th_init) * r * (std::max (mesh->Nx, mesh->Ny) / std::min (mesh->Lx, mesh->Ly)));

    for (int k = 0; k < Np; k++) {

        double coeff = static_cast<double> (k) / static_cast<double> (Np - 1);

        double angle = th_init + coeff * (th_end - th_init);

        double abs = pC->x + r * cos (angle);
        double ord = pC->y + r * sin (angle);

        Point * newPoint = new Point ({abs, ord});
        gamma.push_back (newPoint);

        abs = s * cos (angle);
        ord = s * sin (angle);

        Point * newVec = new Point ({abs, ord});
        normalVec.push_back (newVec);
    }
}


/* Le domaine cartésien lui-même (lorsqu'on ne veut pas de forme intérieure)...
... privé de ses quelques premières rangées de mailles en partant de l'extérieur */

void phi_border (Mesh * mesh)
{
    std::vector<Point *> gamma; // Liste des points de bord
    std::vector<Point *> normalVec; // Liste des vecteurs normaux associés
    gamma.reserve (10 * std::max (mesh->Nx, mesh->Ny));
    normalVec.reserve (10 * std::max (mesh->Nx, mesh->Ny));

    Point basGauche = {mesh->xmin, mesh->ymin};
    Point hautGauche = {mesh->xmin, mesh->ymax};
    Point hautDroite = {mesh->xmax, mesh->ymax};
    Point basDroite = {mesh->xmax, mesh->ymin};

    segment (& basGauche, & hautGauche, 1., gamma, normalVec, mesh);
    segment (& hautGauche, & hautDroite, 1., gamma, normalVec, mesh);
    segment (& hautDroite, & basDroite, 1., gamma, normalVec, mesh);
    segment (& basDroite, & basGauche, 1., gamma, normalVec, mesh);

    for (int i = 0; i < mesh->Nx; i++) {
        for (int j = 0; j < mesh->Ny; j++) {

            double x = mesh->xmin + i * mesh->hx;
            double y = mesh->ymin + j * mesh->hy;

            Point p = {x, y};

            int rankOfMin = -1;
            double distToGamma = 1000.;

            searchMin (& p, gamma, &rankOfMin, &distToGamma);
            double valPhi = + distToGamma; // On suppose + par défaut, si c'est - on change après

            if (rankOfMin < 0)
            {
                std::cout << " Error " << gamma.size () << "\n";
                std::cout << " Error " << rankOfMin << "\n";
                std::cout << " Error " << distToGamma << "\n\n";

            }

            Point * q = gamma [rankOfMin]; // Point qui minimise la distance
            Point * n = normalVec [rankOfMin]; // Vecteur normal en q
            Point vec_qp = {p.x - q->x, p.y - q->y}; // Vecteur qp

            if (prodScal (n, & vec_qp) < 0) // Signe - si nécessaire
                valPhi = - distToGamma;

            mesh->phiOnGrid [j * mesh->Nx + i] = valPhi;
        }
    }

    for (Point * p : gamma)
        delete p;

    for (Point * p : normalVec)
        delete p;

    return;
}


// Cercle centré dans le domaine

void phi_circle (Mesh * mesh)
{
    double xC = (mesh->xmin + mesh->xmax) / 2.;
    double yC = (mesh->ymin + mesh->ymax) / 2.;
    double r = 0.9 * std::min (mesh->Lx / 2., mesh->Ly / 2.);

    for (int i = 0; i < mesh->Nx; i++) {
        for (int j = 0; j < mesh->Ny; j++) {

            double x = mesh->xmin + i * mesh->hx;
            double y = mesh->ymin + j * mesh->hy;

            mesh->phiOnGrid [j * mesh->Nx + i] = std::sqrt ((x - xC) * (x - xC) + (y - yC) * (y - yC)) - r;
        }
    }

    return;
}


// Ellipse centrée dans le domaine, 2 fois + allongée en direction x

void phi_ellipse (Mesh * mesh)
{
    double xC = (mesh->xmin + mesh->xmax) / 2.;
    double yC = (mesh->ymin + mesh->ymax) / 2.;
    double r = 0.9 * std::min ((mesh->xmax - mesh->xmin) / 2., (mesh->ymax - mesh->ymin) / 2.);

    for (int i = 0; i < mesh->Nx; i++) {
        for (int j = 0; j < mesh->Ny; j++) {

            double x = mesh->xmin + i * mesh->hx;
            double y = mesh->ymin + j * mesh->hy;

            mesh->phiOnGrid [j * mesh->Nx + i] = std::sqrt ((x - xC)*(x - xC) + 4*(y - yC)*(y - yC)) - r;
        }
    }

    return;
}


// Carré centré dans le domaine

void phi_square (Mesh * mesh)
{
    std::vector<Point *> gamma; // Liste des points de bord
    std::vector<Point *> normalVec; // Liste des vecteurs normaux associés
    gamma.reserve (10 * std::max (mesh->Nx, mesh->Ny));
    normalVec.reserve (10 * std::max (mesh->Nx, mesh->Ny));

    double xC = (mesh->xmin + mesh->xmax) / 2.;
    double yC = (mesh->ymin + mesh->ymax) / 2.;
    double r = 0.6 * std::min (mesh->Lx / 2., mesh->Ly / 2.);

    Point basGauche = {xC - r, yC - r};
    Point hautGauche = {xC - r, yC + r};
    Point hautDroite = {xC + r, yC + r};
    Point basDroite = {xC + r, yC - r};

    segment (& basGauche, & hautGauche, 1., gamma, normalVec, mesh);
    segment (& hautGauche, & hautDroite, 1., gamma, normalVec, mesh);
    segment (& hautDroite, & basDroite, 1., gamma, normalVec, mesh);
    segment (& basDroite, & basGauche, 1., gamma, normalVec, mesh);

    for (int i = 0; i < mesh->Nx; i++) {
        for (int j = 0; j < mesh->Ny; j++) {

            double x = mesh->xmin + i * mesh->hx;
            double y = mesh->ymin + j * mesh->hy;

            Point p = {x, y};

            int rankOfMin = -1;
            double distToGamma = 1000.;

            searchMin (& p, gamma, &rankOfMin, &distToGamma);
            double valPhi = + distToGamma; // On suppose + par défaut, si c'est - on change après

            if (rankOfMin < 0)
            {
                std::cout << " Error " << gamma.size () << "\n";
                std::cout << " Error " << rankOfMin << "\n";
                std::cout << " Error " << distToGamma << "\n\n";

            }

            Point * q = gamma [rankOfMin]; // Point qui minimise la distance
            Point * n = normalVec [rankOfMin]; // Vecteur normal en q
            Point vec_qp = {p.x - q->x, p.y - q->y}; // Vecteur qp

            if (prodScal (n, & vec_qp) < 0) // Signe - si nécessaire
                valPhi = - distToGamma;

            mesh->phiOnGrid [j * mesh->Nx + i] = valPhi;

            mesh->loc [j * mesh->Nx + i] = OUTSIDE;
            if (valPhi < 0)
                mesh->loc [j * mesh->Nx + i] = INSIDE;
        }
    }

    for (Point * p : gamma)
        delete p;

    for (Point * p : normalVec)
        delete p;

    return;
}


// Etoile de mer à 5 branches centrée dans le domaine

void phi_seastar (Mesh * mesh)
{
    std::vector<Point *> gamma; // Liste des points de bord
    std::vector<Point *> normalVec; // Liste des vecteurs normaux associés
    gamma.reserve (10 * std::max (mesh->Nx, mesh->Ny));
    normalVec.reserve (10 * std::max (mesh->Nx, mesh->Ny));

    // x(t) = xC + [0.5 + 0.2*sin(5t)]*cos(t)

    double xC = (mesh->xmin + mesh->xmax) / 2.;
    double yC = (mesh->ymin + mesh->ymax) / 2.;
    double semi = std::min (mesh->Lx / 2., mesh->Ly / 2.);
    double r = 0.5;
    double dr = 0.2;

    int Np = 2 * floor (2 * PI * r * (std::max (mesh->Nx, mesh->Ny) / std::min (mesh->Lx, mesh->Ly)));

    for (int k = 0; k < Np; k++) {

        double coeff = static_cast<double> (k) / static_cast<double> (Np - 1);

        double theta = coeff * 2 * PI;

        double abs = xC + (r + dr * sin (5 * theta)) * semi * cos (theta);
        double ord = yC + (r + dr * sin (5 * theta)) * semi * sin (theta);

        Point * newPoint = new Point ({abs, ord});
        gamma.push_back (newPoint);

        double xprim = std::cos (5 * theta) * std::cos (theta) - (r + dr * std::sin (5 * theta)) * std::sin (theta);
        double yprim = std::cos (5 * theta) * std::sin (theta) + (r + dr * std::sin (5 * theta)) * std::cos (theta);

        double norm = std::sqrt (xprim * xprim + yprim * yprim);

        abs = yprim / norm;
        ord = - xprim / norm;

        Point * newVec = new Point ({abs, ord});
        normalVec.push_back (newVec);
    }

    for (int i = 0; i < mesh->Nx; i++) {
        for (int j = 0; j < mesh->Ny; j++) {

            double x = mesh->xmin + i * mesh->hx;
            double y = mesh->ymin + j * mesh->hy;

            Point p = {x, y};

            int rankOfMin = -1;
            double distToGamma = 1000.;

            searchMin (& p, gamma, &rankOfMin, &distToGamma);
            double valPhi = + distToGamma; // On suppose + par défaut, si c'est - on change après

            Point * q = gamma [rankOfMin]; // Point qui minimise la distance
            Point * n = normalVec [rankOfMin]; // Vecteur normal en q
            Point vec_qp = {p.x - q->x, p.y - q->y}; // Vecteur qp

            if (prodScal (n, & vec_qp) < 0) // Signe - si nécessaire
                valPhi = - distToGamma;

            mesh->phiOnGrid [j * mesh->Nx + i] = valPhi;
        }
    }

    for (Point * p : gamma)
        delete p;

    for (Point * p : normalVec)
        delete p;

    return;
}


// Visage centré dans le domaine (bien proportionné si le domaine est carré)

void phi_face (Mesh * mesh)
{
    std::vector<Point *> gamma; // Liste des points de bord
    std::vector<Point *> normalVec; // Liste des vecteurs normaux associés
    gamma.reserve (10 * std::max (mesh->Nx, mesh->Ny));
    normalVec.reserve (10 * std::max (mesh->Nx, mesh->Ny));

    return;
}

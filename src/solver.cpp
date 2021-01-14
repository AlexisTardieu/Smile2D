#include "solver.h"

#include <iostream>
#include <Eigen/IterativeLinearSolvers>

void solveCG (Matrix * mat, Vector & sol, Vector & secMember)
{
    std::cout << "STEP 4 : Solving the linear system AX = B to get the solution X..." << std::endl;

    // Verification de la compatibilité des tailles
    if (mat->cols () != sol.size () || mat->rows () != secMember.size ())
    {
        std::cout << "Problem of size compatibility for the system AX = B." << std::endl;
        return;
    }

    // Si la matrice est hermitienne, on utilise le gradient conjugué
    if (mat->isApprox (mat->adjoint ()))
    {
        // Déclare la fonction gradient conjugé
        Eigen::ConjugateGradient<Matrix, Eigen::Lower | Eigen::Upper> CGsolver;

        CGsolver.setTolerance (1e-18); // Précision machine par défaut
        CGsolver.compute (* mat);

        // Calcul de la solution avec le gradient conjugué
        sol = CGsolver.solve (secMember);

        std::cout << "Number of iterations :" << CGsolver.iterations () << std::endl;
        std::cout << "Estimation of error  :" << CGsolver.error ()      << std::endl;
    }

    // Si la matrice n'est pas hermitienne, on utilise le bi-gradient conjugué stabilisé
    else
    {
        Eigen::BiCGSTAB<Matrix> BiCGSTABsolver;

        BiCGSTABsolver.setTolerance (1e-18); // Précision machine par défaut
        BiCGSTABsolver.compute (* mat);

        // Calcul de la solution avec le bi-gradient conjugué stabilisé
        sol = BiCGSTABsolver.solve (secMember);

        std::cout << "Number of iterations <CUICUI>:" << BiCGSTABsolver.iterations () << std::endl;
        std::cout << "Estimation of error  :" << BiCGSTABsolver.error ()      << std::endl;
    }

    std::cout << "STEP 4 : done !" << std::endl << std::endl;

    return;
}

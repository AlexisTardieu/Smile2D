#ifndef LEVELSET_H
#define LEVELSET_H

#include <vector>

struct Point;
struct Mesh;

void segment (Point * pA, Point * pB, double s,
              std::vector<Point *> & gamma, std::vector<Point *> & normalVec, Mesh * mesh);

void arcCercle (Point * pC, double r, double th_init, double th_end, double s,
                std::vector<Point *> & gamma, std::vector<Point *> & normalVec, Mesh * mesh);

void phi_border  (Mesh * mesh);
void phi_circle  (Mesh * mesh);
void phi_ellipse (Mesh * mesh);
void phi_square  (Mesh * mesh);
void phi_seastar (Mesh * mesh);
void phi_face    (Mesh * mesh);

#endif // LEVELSET_H

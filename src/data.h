#ifndef DATA_H
#define DATA_H

#include <iomanip>

#define PI 3.14159265

#define SPC std::setw(15) <<
#define SPCl std::setw(6) <<


typedef enum
{
    INSIDE = 1, // Point intérieur au domaine irrégulier : phi < 0
    OUTSIDE = 0, // Point extérieur au domaine irrégulier : phi >= 0
} Loc;

#endif // DATA_H

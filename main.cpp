#include "Voronoi/Voronoi.hpp"
#include <iostream>
using namespace Vor;

int main(int argc, char ** argv)
{
    std::vector<Point *> points(3);
    points[0] = new Point(0.0, 1.0);
    points[1] = new Point(1.0, 1.0);
    points[2] = new Point(0.0, 2.0);
    Voronoi * v = new Voronoi();

    std::vector<Edge *> edges = v->Generate(points, 5, 5);
    for (auto it = edges.begin(); it != edges.end(); ++it)
    {
        std::cout << "Start: (" << (*it)->start->x << ", " << (*it)->start->y << ")\n";
        std::cout << "End: (" << (*it)->end->x << ", " << (*it)->end->y << ")\n\n";
    }

    return 0;
}
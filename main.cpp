#include "Voronoi/Voronoi.hpp"
using namespace Vor;

int main(int argc, char ** argv)
{
    std::vector<Point *> points(3);
    points[0] = new Point(0.0, 1.0);
    points[1] = new Point(1.0, 1.0);
    points[2] = new Point(0.0, 2.0);
    Voronoi * v = new Voronoi();

    v->Generate(points, 5, 5);

    return 0;
}
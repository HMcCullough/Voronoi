#include "Voronoi/Voronoi.hpp"
using namespace Vor;

int main(int argc, char ** argv)
{
    std::vector<Point> points(3);
    points[0] = { 0.0, 1.0 };
    points[1] = { 1.0, 1.0 };
    points[2] = { 0.0, 2.0 };
    Voronoi * v = new Voronoi(points);

    return 0;
}
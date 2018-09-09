#include <vector>
#include <queue>

namespace Vor
{
    #pragma region Point Class Prototype
    struct Point
    {
        double x;
        double y;
        Point();
        Point(double x, double y);
    };
    #pragma endregion

    #pragma region Event Class Prototype
    typedef enum { SiteEvent, CircleEvent } EventType;
    class Event
    {
    public:
        EventType type;
        Point site;

        Event(Point &p, const EventType &type);

        struct EComparer
        {
            bool operator()(Event * e1, Event * e2)
            {
                return ((e1->site.y < e2->site.y) ||
                       ((e1->site.y == e2->site.y) && (e1->site.x < e2->site.x)));
            }
        };
    };
    #pragma endregion

    #pragma region Edge Class Prototype
    class Edge
    {
        Point start, end;
        Edge(Point &start, Point &end);
    };
    #pragma endregion

    #pragma region Voronoi Class Prototype
    class Voronoi
    {
    public:
        Voronoi();
        Voronoi(std::vector<Point> &points);

        std::vector<Edge> Generate();

    private:
        double sweepLineY;
        std::priority_queue<Event *, std::vector<Event *>, Event::EComparer> eventQ;

        void ProcessSiteEvent(Event * event);
        void ProcessCircleEvent(Event * event);
    };
    #pragma endregion
}

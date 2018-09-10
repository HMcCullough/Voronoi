#include <vector>
#include <queue>
#include <set>
#include <list>

namespace Vor
{
    static const double EPSILON = 0.1;
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
        Point * site;

        Event(Point * p, const EventType &type);

        struct EComparer
        {
            bool operator()(const Event * e1, const Event * e2)
            {
                return ((e1->site->y < e2->site->y) ||
                       ((e1->site->y - e2->site->y < EPSILON) && (e1->site->x < e2->site->x)));
            }
        };
    };
    #pragma endregion

    #pragma region Edge Class Prototype
    class Edge
    {
    public:
        Point * start, * left, * right;
        Edge(Point * start, Point * left, Point * right);
    };
    #pragma endregion

    #pragma region BeachLine Class Prototype
    // Per Stephen Fortune's algorithm, the beach-line sequence will be represented as a binary tree where the leaf nodes
    // contain the parabola information and the inner nodes contain edge information.
    class Parabola
    {
    public:
        friend class BeachLine;
        Parabola();
        Parabola(Point * point);

        // Getters
        bool isLeaf() const { return (_left == nullptr && _right == nullptr); }

        Point * getSite() const { return _site; }
        Parabola * getLeft() const { return _left; }
        Parabola * getRight() const { return _right; }

        Edge * getEdge() const { return _edge; }
        Event * getCircleEvent() const { return _event; }

        // Setters
        void setLeft(Parabola * par) { _left = par; }
        void setRight(Parabola * par) { _right = par; }
        void setEdge(Edge * edge) { _edge = edge; }
        void setCircleEvent(Event * cEvent) { _event = cEvent; }
    private:
        Edge * _edge;
        Event * _event;
        Point * _site;
        Parabola * _left, * _right, * _parent;

        bool _isLeaf;
    };

    class BeachLine
    {
    public:
        BeachLine();

        bool isEmpty() const;
        Parabola * GetParabolaByX(double x);

        Parabola * getRoot() const { return _root; }
        void setRoot(Parabola * par) { _root = par; }

    private:
        Parabola * _root;
    };
    #pragma endregion

    #pragma region Voronoi Class Prototype
    class Voronoi
    {
    public:
        Voronoi();

        std::vector<Edge *> Generate(std::vector<Point *> &points, int width, int height);

    private:
        double sweepLineY;
        std::priority_queue<Event *, std::vector<Event *>, Event::EComparer> eventQ;
        std::set<Event *> deletedEvents;
        std::vector<Edge *> edges;
        std::list<Point *> points;
        BeachLine beachline;

        int width, height;

        void InsertParabola(Point * point);
        void RemoveParabola(Event * event);
        void FixEdges();
    };
    #pragma endregion
}

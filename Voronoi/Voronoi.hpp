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
    class Parabola;
    typedef enum { SiteEvent, CircleEvent } EventType;
    class Event
    {
    public:
        Parabola * parabola;
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
        Parabola * parabola;
        Point * start, * left, * right, * end;
        Edge * neighbor;
        Edge(Point * start, Point * left, Point * right);

        Point * GetIntersection(Edge * other);
    private:
        Point GetDirection();
    };
    #pragma endregion

    #pragma region Parabola Class Prototype
    // Per Stephen Fortune's algorithm, the beach-line sequence will be represented as a binary tree where the leaf nodes
    // contain the parabola information and the inner nodes contain edge information.
    class Parabola
    {
    public:
        friend class BeachLine;

        Parabola();
        Parabola(Point * point);

        // Member Functions
        double GetY(double x, double sweepLineY);
        Parabola * GetLeftParent();
        Parabola * GetRightParent();
        Parabola * GetLeftChild();
        Parabola * GetRightChild();

        // Getters
        bool isLeaf() const { return (_left == nullptr && _right == nullptr); }

        Point * getSite() const { return _site; }
        Parabola * getLeft() { return _left; }
        Parabola * getRight() { return _right; }
        Parabola * getParent() const { return _parent; }

        Edge * getEdge() const { return _edge; }
        Event * getCircleEvent() const { return _event; }

        // Setters
        void setLeft(Parabola * par)
        {
            _left = par;
            par->_parent = this;
        }
        void setRight(Parabola * par)
        {
            _right = par;
            par->_parent = this;
        }
        void setEdge(Edge * edge) { _edge = edge; }
        void setCircleEvent(Event * cEvent) { _event = cEvent; }
        
    private:
        Edge * _edge;
        Event * _event;
        Point * _site;
        Parabola * _left, * _right, * _parent;

        // Useful for calculating parabola intersections and specific coordinates on parabola
        typedef struct
        {
            // In the form f(x)=ax^2+bx+c
            double a, b, c;
        } ParabolaProperties;

        double GetXOfEdge(double sweepLineY);
        ParabolaProperties GetProperties(double sweepLineY);
    };
    #pragma endregion

    #pragma region BeachLine Class Prototype
    class BeachLine
    {
    public:
        BeachLine();

        bool isEmpty() const;
        Parabola * GetParabolaByX(double x, double sweepLineY);

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
        double _sweepLineY;
        std::priority_queue<Event *, std::vector<Event *>, Event::EComparer> _eventQ;
        std::set<Event *> _deletedEvents;
        std::vector<Edge *> _edges;
        std::list<Point *> _points;
        BeachLine _beachline;

        int width, height;

        void InsertParabola(Point * point);
        void RemoveParabola(Event * event);
        void CheckCircle (Parabola * par);
        void FixEdges();
    };
    #pragma endregion
}

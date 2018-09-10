/*
For inserting new events we will need a priority queue. This is because we need to keep the sites
sorted "lexicograpically" according to the algroithm and a priority queue will do this quickly.
*/
#include "Voronoi.hpp"
namespace Vor
{
    #pragma region Point Class Definition
    Point::Point() : x(0.0), y(0.0) {}
    Point::Point(double x, double y) : x(x), y(y) {}
    #pragma endregion

    #pragma region Event Class Definition
    Event::Event(Point * p, const EventType &type) : site(p), type(type) {}
    #pragma endregion

    #pragma region Edge Class Definition
    Edge::Edge(Point * start, Point * left, Point * right) : start(start), left(left), right(right) {}
    #pragma endregion

    #pragma region BeachLine Class Definition
    // Parabola definitions
    #pragma region Parabola Class Definition
    Parabola::Parabola()
    {
        _isLeaf = false;
        _edge = nullptr;
        _event = nullptr;
        _site = nullptr;
        _left = _right = _parent = nullptr;
    }

    Parabola::Parabola(Point * point)
    {
        _isLeaf = false;
        _edge = nullptr;
        _event = nullptr;
        _site = point;
        _left = _right = _parent = nullptr;
    }
    #pragma endregion

    // BeachLine definitions
    BeachLine::BeachLine() : _root(nullptr) {}
    bool BeachLine::isEmpty() const { return (_root == nullptr); }

    Parabola * BeachLine::GetParabolaByX(double x)
    {
        Parabola * current = _root;
        while (!current->isLeaf())
        {
            if (current->_site->x > x)
                current = current->_left;
            else
                current = current->_right;
        }
        return current;
    }
    #pragma endregion

    #pragma region Voronoi Class Definition
    Voronoi::Voronoi() {}

    std::vector<Edge *> Voronoi::Generate(std::vector<Point *> &points, int width, int height)
    {
        this->width = width;
        this->height = height;

        sweepLineY = 0.0;
        for (auto it = points.begin(); it != points.end(); ++it)
            eventQ.push(new Event(*it, EventType::SiteEvent));

        // Here we process the event queue until it is empty
        Event * event = nullptr;
        while (!eventQ.empty())
        {
            // We start by peeking at the event queue and setting the sweepLineY to the event site's y coord
            // This is guaranteed to be the highest point since the y-coord is the priority in eventQ.
            event = eventQ.top();
            eventQ.pop();
            sweepLineY = event->site->y;

            // First we check to see if the event has already been removed. This occurs when circle events are encountered.
            if (deletedEvents.find(event) != deletedEvents.end())
            {
                deletedEvents.erase(event);
                delete event;
                continue;
            }

            // Process Site Event
            if (event->type == EventType::SiteEvent)
                InsertParabola(event->site);
            // Process Circle Event
            else if (event->type == EventType::CircleEvent)
                RemoveParabola(event);
            // Dispose of the event
            delete event;

            // Whether we encounter a site or circle event, we will always check to see if we have added a new triple of parabolic
            // arcs on the beach line that may lead to a future circle event. It is also possible that we will remove a future
            // circle event from consideration if the triple of parabolic arcs that led to it no longer exists after the event.
            // This step occurs in process functions
        }

        // Fortune's Algorithm produces a "distorted" but topographically equivalent version of the Voronoi diagram. This is directly
        // taken from CMSC 754 lecture notes at http://people.math.gatech.edu/~randall/Algs07/mount.pdf but I could not find the author's name.
        FixEdges();

        return edges;
    }

    // If the next event that the sweep line encounters is a site event, we simply insert the new site into our list
    // of sites in the order in which the corresponding parabolic arc appears on the beach line. We then record the
    // fact that we have encountered a new edge in the diagram.
    void Voronoi::InsertParabola(Point * point)
    {
        if (!beachline.getRoot())
        {
            beachline.setRoot(new Parabola(point));
            return;
        }

        // This is a degenerate case where the first two nodes are at the same y-coordinate
        // In this case we will be making the root node an inner node with two leaf children for each site.
        if (beachline.getRoot()->isLeaf() && beachline.getRoot()->getSite()->y - point->y < EPSILON)
        {
            // _root is no longer a leaf
            Point * firstSite = beachline.getRoot()->getSite();

            // Add a child foreach site to _root
            beachline.getRoot()->setLeft(new Parabola(firstSite));
            beachline.getRoot()->setRight(new Parabola(point));

            // If these sites are at the same y-coord then the edge beginning will be the mid-x starting at the
            // highest possible point (height)
            Point * midPoint = new Point((firstSite->x + point->x) / 2.0, height);
            points.push_back(midPoint);

            // We do not need to determine which is on the left or right because the priority queue sorts
            // lexicographically which means if they are at the same y-coord they will pop out with the lower
            // x-coord first
            beachline.getRoot()->setEdge(new Edge(midPoint, firstSite, point));
            edges.push_back(beachline.getRoot()->getEdge());
            return;
        }

        Parabola * parabolaAboveX = beachline.GetParabolaByX(point->x);
        if (parabolaAboveX->getCircleEvent())
        {
            deletedEvents.insert(parabolaAboveX->getCircleEvent());
            parabolaAboveX->setCircleEvent(nullptr);
        }
    }

    // If the next event that the sweep line encounters is a circle event, we record the fact that we have encountered
    // a vertex in the diagram and that this vertex is the endpoint of the edges corresponding to the two breakpoints
    // that have come together. We also record the new edge corresponding to the new breakpoint that results from the
    // circle event.
    void Voronoi::RemoveParabola(Event * event)
    {
    }

    void FixEdges()
    {
    }
    #pragma endregion
}
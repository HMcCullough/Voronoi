/*
The algorithm

We now have everything we need to know to understand Fortune's algorithm. To detect edges and vertices in the diagram, it is enough to find the appearance and disappearance
of parabolic arcs in the beach line. We will therefore keep track of the beach line by imagining that we walk along it from left to right and record the order of the sites
that produce its constituent parabolic arcs. We know that this order will not change until the sweep line reaches either a site event or circle event. Also, the breakpoints
are implicitly recorded by noting the adjacency of the parabolic arcs on the beach line.

In this way, the diagram is constructed by considering the finite sequence of events. Shown below is the sequence of events that computes the Voronoi diagram of the
collection of sites shown. Circle events are indicated by green dots.
*/

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
    Edge::Edge(Point &start, Point &end) : start(start), end(end) {}
    #pragma endregion

    #pragma region BeachLine Class Definition
    // Parabola definitions
    Parabola::Parabola()
    {
        _isLeaf = false;
        _edge = nullptr;
        _event = nullptr;
        _point = nullptr;
        _left = _right = _parent = nullptr;
    }

    Parabola::Parabola(Point * point)
    {
        _isLeaf = false;
        _edge = nullptr;
        _event = nullptr;
        _point = point;
        _left = _right = _parent = nullptr;
    }

    // BeachLine definitions
    BeachLine::BeachLine() : _root(nullptr) {}
    bool BeachLine::isEmpty() const { return (_root == nullptr); }

    void BeachLine::insert(Event * event)
    {
    }
    #pragma endregion

    #pragma region Voronoi Class Definition
    Voronoi::Voronoi()
    {
    }

    std::vector<Edge> Voronoi::Generate(std::vector<Point *> &points)
    {
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
                ProcessSiteEvent(event);
            // Process Circle Event
            else if (event->type == EventType::CircleEvent)
                ProcessCircleEvent(event);
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
    void Voronoi::ProcessSiteEvent(Event * event)
    {
    }

    // If the next event that the sweep line encounters is a circle event, we record the fact that we have encountered
    // a vertex in the diagram and that this vertex is the endpoint of the edges corresponding to the two breakpoints
    // that have come together. We also record the new edge corresponding to the new breakpoint that results from the
    // circle event.
    void Voronoi::ProcessCircleEvent(Event * event)
    {
    }

    void FixEdges()
    {
    }
    #pragma endregion
}
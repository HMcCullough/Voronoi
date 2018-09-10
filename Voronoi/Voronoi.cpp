/*
For inserting new events we will need a priority queue. This is because we need to keep the sites
sorted "lexicograpically" according to the algroithm and a priority queue will do this quickly.
*/
#include "Voronoi.hpp"
#include <cmath>
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

    double Parabola::GetXOfEdge(double sweepLineY)
    {
        // left is the rightmost child in the left subtree of this
        // right is the leftmost child in the right subtree of this
        // These represent the parabolas that form the current node's edge
        Parabola * left = nullptr, * right = nullptr;
        for (left = _left; left->_right != nullptr; left = left->_right);
        for (right = _right; right->_left != nullptr; right = right->_left);

        typedef struct
        {
            // In the form f(x)=ax^2+bx+c
            double a, b, c;
        } ParabolaProperties;
        // NOTE : This math comes directly from the directrix/focus relationship of parabolas:
        //      f(x)=((x - p_x)^2) / (2 * (p_y - sweepLineY))+(p_y + sweepLineY) / 2
        // dp is used to "simplify" the calculations
        double dp = 2.0 * (left->_site->y - sweepLineY);
        ParabolaProperties leftPara = {
            1.0 / dp,
            -2.0 * left->_site->x / dp,
            sweepLineY + dp / 4.0 + left->_site->x * left->_site->x / dp
        };
        
        dp = 2.0 * (right->_site->y - sweepLineY);
        ParabolaProperties rightPara = {
            1.0 / dp,
            -2.0 * right->_site->x / dp,
            sweepLineY + dp / 4.0 + right->_site->x * right->_site->x / dp
        };

        // The roots of the avgPara will be the x-coordinates of the intersections of the left and right parabolas.
        ParabolaProperties avgPara = { leftPara.a - rightPara.a, leftPara.b - rightPara.b, leftPara.c - rightPara.c };
        double discriminant = std::sqrt((avgPara.b * avgPara.b) - 4.0 * avgPara.a * avgPara.c);
        double x0 = (-avgPara.b + discriminant) / (2 * avgPara.a),
               x1 = (-avgPara.b - discriminant) / (2 * avgPara.a);

        // Since the parabola have two intersects we need to select the one between the two parabola
        // the intersects will (roughly) be in between the parabola and one on the left or right depending on the y-coords of the sites.
        // If the left parabola is lower then it will be on the left of both parabola, meaning the in-between intersection will have a higher x-coord.
        // Otherwise, it will be on the right meaning the in-between intersection will have a smaller x-coord.
        // Play here to see for yourself: https://www.desmos.com/calculator/tuyioagkxh
        return (left->_site->y < right->_site->y) ? std::max(x0, x1) : std::min(x0, x1);
    }
    #pragma endregion

    // BeachLine definitions
    BeachLine::BeachLine() : _root(nullptr) {}
    bool BeachLine::isEmpty() const { return (_root == nullptr); }

    Parabola * BeachLine::GetParabolaByX(double x, double sweepLineY)
    {
        Parabola * current = _root;
        // We search for the parabola we are under by following the edge connections
        // in the inner nodes until we hit a leaf.
        while (!current->isLeaf())
        {
            double edgeX = current->GetXOfEdge(sweepLineY);
            if (edgeX > x)
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

        _sweepLineY = 0.0;
        for (auto it = points.begin(); it != points.end(); ++it)
            _eventQ.push(new Event(*it, EventType::SiteEvent));

        // Here we process the event queue until it is empty
        Event * event = nullptr;
        while (!_eventQ.empty())
        {
            // We start by peeking at the event queue and setting the sweepLineY to the event site's y coord
            // This is guaranteed to be the highest point since the y-coord is the priority in eventQ.
            event = _eventQ.top();
            _eventQ.pop();
            _sweepLineY = event->site->y;

            // First we check to see if the event has already been removed. This occurs when circle events are encountered.
            if (_deletedEvents.find(event) != _deletedEvents.end())
            {
                _deletedEvents.erase(event);
                delete event;
                continue;
            }

            // Whether we encounter a site or circle event, we will always check to see if we have added a new triple of parabolic
            // arcs on the beach line that may lead to a future circle event. It is also possible that we will remove a future
            // circle event from consideration if the triple of parabolic arcs that led to it no longer exists after the event.
            // This step occurs in process functions
            // Process Site Event
            if (event->type == EventType::SiteEvent)
                InsertParabola(event->site);
            // Process Circle Event
            else if (event->type == EventType::CircleEvent)
                RemoveParabola(event);
            // Dispose of the event
            delete event;
        }

        // Fortune's Algorithm produces a "distorted" but topographically equivalent version of the Voronoi diagram. This is directly
        // taken from CMSC 754 lecture notes at http://people.math.gatech.edu/~randall/Algs07/mount.pdf but I could not find the author's name.
        FixEdges();

        return _edges;
    }

    // If the next event that the sweep line encounters is a site event, we simply insert the new site into our list
    // of sites in the order in which the corresponding parabolic arc appears on the beach line. We then record the
    // fact that we have encountered a new edge in the diagram.
    void Voronoi::InsertParabola(Point * point)
    {
        if (!_beachline.getRoot())
        {
            _beachline.setRoot(new Parabola(point));
            return;
        }

        // This is a degenerate case where the first two nodes are at the same y-coordinate
        // In this case we will be making the root node an inner node with two leaf children for each site.
        if (_beachline.getRoot()->isLeaf() && _beachline.getRoot()->getSite()->y - point->y < EPSILON)
        {
            // _root is no longer a leaf
            Point * firstSite = _beachline.getRoot()->getSite();

            // Add a child foreach site to _root
            _beachline.getRoot()->setLeft(new Parabola(firstSite));
            _beachline.getRoot()->setRight(new Parabola(point));

            // If these sites are at the same y-coord then the edge beginning will be the mid-x starting at the
            // highest possible point (height)
            Point * midPoint = new Point((firstSite->x + point->x) / 2.0, height);
            _points.push_back(midPoint);

            // We do not need to determine which is on the left or right because the priority queue sorts
            // lexicographically which means if they are at the same y-coord they will pop out with the lower
            // x-coord first
            _beachline.getRoot()->setEdge(new Edge(midPoint, firstSite, point));
            _edges.push_back(_beachline.getRoot()->getEdge());
            return;
        }

        // We get the parabola above this point's x-coord and check if it has a future circle event.
        Parabola * parabolaAboveX = _beachline.GetParabolaByX(point->x, _sweepLineY);
        if (parabolaAboveX->getCircleEvent())
        {
            _deletedEvents.insert(parabolaAboveX->getCircleEvent());
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
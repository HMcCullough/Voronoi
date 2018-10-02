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
    Edge::Edge(Point * start, Point * left, Point * right) : start(start), left(left), right(right), end(nullptr), neighbor(nullptr), parabola(nullptr) {}

    Point Edge::GetDirection()
    {
        return { -(left->y - right->y), left->x - right->x };
    }

    Point * Edge::GetIntersection(Edge * other)
    {
        Point thisDir = GetDirection(), otherDir = other->GetDirection(), * intersection = nullptr;
        double u = (start->y * otherDir.x + otherDir.y * other->start->x - other->start->y * otherDir.x - otherDir.y * start->x ) / (thisDir.x * otherDir.y - thisDir.y * otherDir.x);
        double v = (start->x + thisDir.x * u - other->start->x) / otherDir.x;
        
        if (u >= 0 && v >= 0)
            intersection = new Point(start->x + thisDir.x * u, start->y + thisDir.y * u);
        return intersection;
    }
    #pragma endregion

    #pragma region Parabola Class Definition
    Parabola::Parabola()
    {
        _edge = nullptr;
        _event = nullptr;
        _site = nullptr;
        _left = _right = _parent = nullptr;
    }

    Parabola::Parabola(Point * point)
    {
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
        Parabola * left = GetLeftChild(), * right = GetRightChild();

        ParabolaProperties leftPara = left->GetProperties(sweepLineY);
        ParabolaProperties rightPara = right->GetProperties(sweepLineY);

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

    // Return the y-coordinate of parabola at x defined by sweepLineY as the directrix and _site as the focus
    double Parabola::GetY(double x, double sweepLineY)
    {
        ParabolaProperties para = GetProperties(sweepLineY);
        return para.a * x * x + para.b * x + para.c;
    }

    Parabola::ParabolaProperties Parabola::GetProperties(double sweepLineY)
    {
        // NOTE : This math comes directly from the directrix/focus relationship of parabolas:
        //      f(x)=((x - p_x)^2) / (2 * (p_y - sweepLineY))+(p_y + sweepLineY) / 2
        // dp is just to "simplify" the equations
        double dp = 2.0 * (_site->y - sweepLineY);
        return {
            1.0 / dp,
            -2.0 * _site->x / dp,
            sweepLineY + dp / 4.0 + _site->x * _site->x / dp
        };
    }

    Parabola * Parabola::GetLeftParent()
    {
        Parabola * leftParent = this->_parent, * last = this;
        while (leftParent->getLeft() == last)
        {
            if (leftParent->getParent() == nullptr)
                return nullptr;
            last = leftParent;
            leftParent = leftParent->getParent();
        }
        return leftParent;
    }

    Parabola * Parabola::GetRightParent()
    {
        Parabola * rightParent = this->_parent, * last = this;
        while (rightParent->getRight() == last)
        {
            if (rightParent->getParent() == nullptr)
                return nullptr;
            last = rightParent;
            rightParent = rightParent->getParent();
        }
        return rightParent;
    }

    Parabola * Parabola::GetLeftChild()
    {
        Parabola * left = nullptr;
        for (left = _left; !left->isLeaf(); left = left->_right);
        return left;
    }

    Parabola * Parabola::GetRightChild()
    {
        Parabola * right = nullptr;
        for (right = _right; !right->isLeaf(); right = right->_left);
        return right;
    }
    #pragma endregion

    #pragma region BeachLine Class Definition
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
        this->_width = width;
        this->_height = height;

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
        FixEdges(_beachline.getRoot());

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
            Point * midPoint = new Point((firstSite->x + point->x) / 2.0, _height);
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

        Point * start = new Point(point->x, parabolaAboveX->GetY(point->x, _sweepLineY));
        _points.push_back(start);

        Edge * leftEdge = new Edge(start, parabolaAboveX->getSite(), point),
             * rightEdge = new Edge(start, point, parabolaAboveX->getSite());

        leftEdge->neighbor = rightEdge;
        _edges.push_back(leftEdge);

        Parabola * pLeft = new Parabola(parabolaAboveX->getSite()),
                 * pNew = new Parabola(point),
                 * pRight = new Parabola(parabolaAboveX->getSite());
        
        /*
                (parabolaAboveX)      ===>      (rightEdge)
                                                /         \
                                        (leftEdge)      (pRight)
                                        /        \
                                    (pLeft)     (pNew)
        */
        // The above diagram shows what happens to parabolaAboveX
        parabolaAboveX->setRight(pRight);
        parabolaAboveX->setLeft(new Parabola());
        parabolaAboveX->getLeft()->setLeft(pLeft);
        parabolaAboveX->getLeft()->setRight(pNew);

        parabolaAboveX->setEdge(rightEdge);
        parabolaAboveX->getLeft()->setEdge(leftEdge);

        CheckCircle(pLeft);
        CheckCircle(pRight);
    }

    // If the next event that the sweep line encounters is a circle event, we record the fact that we have encountered
    // a vertex in the diagram and that this vertex is the endpoint of the edges corresponding to the two breakpoints
    // that have come together. We also record the new edge corresponding to the new breakpoint that results from the
    // circle event.
    void Voronoi::RemoveParabola(Event * event)
    {
        Parabola * parabolaToRemove = event->parabola;
        Parabola * leftParent = parabolaToRemove->GetLeftParent(),
                 * rightParent = parabolaToRemove->GetRightParent();
        Parabola * leftSib = leftParent->GetLeftChild(),
                 * rightSib = rightParent->GetRightChild();

        if (leftSib->getCircleEvent())
        {
            _deletedEvents.insert(leftSib->getCircleEvent());
            leftSib->setCircleEvent(nullptr);
        }
        if (rightSib->getCircleEvent())
        {
            _deletedEvents.insert(rightSib->getCircleEvent());
            rightSib->setCircleEvent(nullptr);
        }

        // Note that the current point is the circumcenter of parabolaToRemove, leftSib, and rightSib (see CheckCircle)
        Point * circumcenter = new Point(event->site->x, parabolaToRemove->GetY(event->site->x, _sweepLineY));
        _points.push_back(circumcenter);

        // This is the new edge formed at the circle event since edges on either side will end here and "collide" into one
        Edge * newEdge = new Edge(circumcenter, leftSib->getSite(), rightSib->getSite());
        leftParent->getEdge()->end = circumcenter;
        rightParent->getEdge()->end = circumcenter;

        /*
                (leftParent)                            (leftParent)
                /          \                            /          \
            (lefSib)   (rightParent)      ===>     (leftSib)   (rightSib)
                       /           \
            (parabolatoRemove) (rightSib)
            
            OR

                            (rightParent)                   (rightParent)
                            /           \                   /           \
                    (leftParent)    (rightSib)   ===>   (leftSib)    (rightSib)
                    /          \
            (leftSib)   (parabolaToRemove)
        */
        // The above diagram shows what is happening in this function and as you can see we need a reference to the higher
        // of leftParent and rightParent so that it can store the new edge. NOTE : the diagram not the only case but it is the
        // simplest case to show, rightParent and leftParent might be farther up ther tree.
        Parabola * higher = nullptr, * par = parabolaToRemove;
        while (par != _beachline.getRoot())
        {
            par = par->getParent();
            if (par == leftParent)
                higher = leftParent;
            else if (par == rightParent)
                higher = rightParent;
        }
        higher->setEdge(newEdge);
        _edges.push_back(newEdge);

        /*
                    (rightParent)                               (rightParent)
                    /           \                               /           \
                 (....)      (rightSib)                     (....)       (rightSib)
                 /    \                         ===>        /    \
            (....)  (leftParent)                        (....)   (leftSib)
                    /          \
                (leftSib)  (parabolaToRemove)

                    (leftParent)                                (leftParent)
                    /          \                                /          \
                (leftSib)     (....)                        (leftSib)     (....)
                              /    \            ===>                      /    \
                    (rightParent)  (....)                           (rightSib) (....)
                    /           \
          (parabolaToRemove)  (rightSib)
        */
        // The above diagram shows the other 2 cases which are all represented by the code below. Think of the first diagrams as
        // cases 1.2 and 2.1 respectively and the second group of diagrams as 2.2 and 1.1 respectively.
        Parabola * gParent = parabolaToRemove->getParent()->getParent();
        if (parabolaToRemove->getParent()->getLeft() == parabolaToRemove)
        {
            if (gParent->getLeft() == parabolaToRemove->getParent())
                gParent->setLeft(parabolaToRemove->getParent()->getRight()); // Case 1.1
            else if (gParent->getRight() == parabolaToRemove->getParent())
                gParent->setRight(parabolaToRemove->getParent()->getRight()); // Case 1.2
        }
        else
        {
            if (gParent->getLeft() == parabolaToRemove->getParent())
                gParent->setLeft(parabolaToRemove->getParent()->getLeft()); // Case 2.1
            else if (gParent->getRight() == parabolaToRemove->getParent())
                gParent->setRight(parabolaToRemove->getParent()->getLeft()); // Case 2.2
        }

        delete parabolaToRemove->getParent();
        delete parabolaToRemove;

        CheckCircle(leftSib);
        CheckCircle(rightSib);
    }
    
    void Voronoi::CheckCircle (Parabola * par)
    {
        // This will get the parents to which par is the leftmost child in the right subtree and the rightmost child in the left subtree respectively
        // (leftParent => first parent on the left) and (rightParent => first parent on the right)
        // Topographically this represent the edge formed by intersection of the parabola to the left and right respectively 
        Parabola * leftParent = par->GetLeftParent(),
                 * rightParent = par->GetRightParent();
        // If the parent exists then grab its respective child. NOTE : leftParent->GetRightChild() == rightParent->GetLeftChild() == par
        Parabola * leftSib = leftParent ? leftParent->GetLeftChild() : nullptr,
                 * rightSib = rightParent ? rightParent->GetRightChild() : nullptr;

        if (!leftSib || !rightSib || leftSib->getSite() == rightSib->getSite())
            return;

        // The intersection point represents the points to which leftSib and rightSib squeeze par. This means that it is equidistant
        // from leftSib, par, and rightSib (since this point is the meeting of three cells in the diagram).
        Point * intersection = leftParent->getEdge()->GetIntersection(rightParent->getEdge());
        if (intersection)
            _points.push_back(intersection);
        else
            return;

        // Since this point is equidistant from all of the previously listed parabola, we simple need the distance from one to detect
        // when the future circle event will occur.
        double dx = leftSib->getSite()->x - intersection->x,
               dy = leftSib->getSite()->y - intersection->y;
        double d = std::sqrt(dx * dx + dy * dy);

        // We check to see if the circle event lies below the sweep-line and if so return
        if (intersection->y - d > _sweepLineY)
            return;

        // Otherwise we make the new event and push it to the queue
        Event * cEvent = new Event(new Point(intersection->x, intersection->y - d), EventType::CircleEvent);
        _points.push_back(cEvent->site);
        _eventQ.push(cEvent);

        par->setCircleEvent(cEvent);
        cEvent->parabola = par;
    }

    void Voronoi::FixEdges(Parabola * par)
    {
        // Delete leaf nodes since they will not make any more edges
        if (par && par->isLeaf())
        {
            delete par;
            return;
        }

        // We want to display the edges outside the bounds of the screen so that it looks like the edges continue past
        double mx = par->getEdge()->start->x;
        Point dir = par->getEdge()->GetDirection();
        if (dir.x != 0.0)
            mx = (dir.x > 0.0) ? mx = std::max(_width, par->getEdge()->start->x + 10) : mx = std::min(0.0, par->getEdge()->start->x - 10);

        // The end point will be the edge continued until it hits the edge of the screen
        Point * end = new Point(mx, (mx * dir.y + par->getEdge()->start->y));
        par->getEdge()->end = end;
        _points.push_back(end);

        // Recursively correct all other edges
        FixEdges(par->getLeft());
        FixEdges(par->getRight());
        delete par;
    }
    #pragma endregion
}
/*
The algorithm

We now have everything we need to know to understand Fortune's algorithm. To detect edges and vertices in the diagram, it is enough to find the appearance and disappearance
of parabolic arcs in the beach line. We will therefore keep track of the beach line by imagining that we walk along it from left to right and record the order of the sites
that produce its constituent parabolic arcs. We know that this order will not change until the sweep line reaches either a site event or circle event. Also, the breakpoints
are implicitly recorded by noting the adjacency of the parabolic arcs on the beach line.

If the next event that the sweep line encounters is a site event, we simply insert the new site into our list of sites in the order in which the corresponding parabolic arc
appears on the beach line. We then record the fact that we have encountered a new edge in the diagram.

If the next event that the sweep line encounters is a circle event, we record the fact that we have encountered a vertex in the diagram and that this vertex is the endpoint
of the edges corresponding to the two breakpoints that have come together. We also record the new edge corresponding to the new breakpoint that results from the circle event.

Whether we encounter a site or circle event, we will always check to see if we have added a new triple of parabolic arcs on the beach line that may lead to a future circle
event. It is also possible that we will remove a future circle event from consideration if the triple of parabolic arcs that led to it no longer exists after the event.

In this way, the diagram is constructed by considering the finite sequence of events. Shown below is the sequence of events that computes the Voronoi diagram of the
collection of sites shown. Circle events are indicated by green dots.
*/

/*
For inserting new events we will need a priority queue. This is because
*/

namespace Voronoi
{
    typedef struct
    {
        double x, y;
    } point;

    class Voronoi
    {
    public:
        Voronoi()
        {
        }

    private:
    };
}
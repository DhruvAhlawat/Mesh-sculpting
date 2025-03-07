
class HalfEdge
{
    public:
    HalfEdge *pair, *next;
    Vertex *head;
    Face *face;

};


class Vertex 
{
    public:
    HalfEdge *halfEdge;
};


class Face 
{
    public:
    HalfEdge *halfEdge;
};
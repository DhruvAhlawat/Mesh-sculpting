
#include <glm/glm.hpp>
#include <vector>
using namespace glm;


class mesh 
{
    public:
    std::vector<vec3> vertexPositions;
    std::vector<ivec3> triangles;
    std::vector<vec3> normals;
    std::vector<HalfEdge> halfEdges;
    std::vector<Vertex> verts;
    std::vector<Face> faces;
};


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
    int num_sides;
};
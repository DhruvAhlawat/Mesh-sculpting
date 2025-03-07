
#include <glm/glm.hpp>
#include <vector>
using namespace glm;


class mesh 
{
    public:
    std::vector<vec3> vertexPositions;
    std::vector<ivec3> triangles;
    std::vector<vec3> normals;
    std::vector<ivec2> edges; 

    std::vector<HalfEdge> halfEdges;
    std::vector<Vertex> verts;
    std::vector<Face> faces;
    //triangulates the mesh for rasterization. 
    void triangulateMesh();

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
    int id; //index of the vertex in the vertexPositions array.
    HalfEdge *halfEdge;
};

class Face 
{
    public:
    HalfEdge *halfEdge;
    int num_sides;
};
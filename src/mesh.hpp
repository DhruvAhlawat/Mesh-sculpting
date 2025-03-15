
#include <glm/glm.hpp>
#include <vector>
#include <cmath>

#include <iostream>
#include <fstream>
#include <sstream>

using namespace glm;


class Vertex;
class Face;

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

    Vertex(int ID)
    {
        id = ID;
    }

    Vertex()
    {
        id = -1;
    }
};

class Face 
{
    public:
    HalfEdge *halfEdge;
    int num_sides;
};


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
    void getMeshObj();
    // mesh(int total_verts);
};

mesh createGrid(int m, int n);
mesh generateSphere(int m, int n);
mesh generateCube(int m, int n, int o);
mesh loadOBJ(const std::string& filename);
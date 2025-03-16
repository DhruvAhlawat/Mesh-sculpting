
#include <glm/glm.hpp>
#include <vector>
#include <cmath>
#include <unordered_map>
#include<set>
#include <map>
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

    HalfEdge(HalfEdge *Pair, HalfEdge *Next, Vertex *Head, Face *Face)
    {
        pair = Pair;
        next = Next;
        head = Head;
        face = Face;
    }
    HalfEdge()
    {
        pair = nullptr;
        next = nullptr;
        head = nullptr;
        face = nullptr;
    }
};

class Vertex 
{
    public:
    int id; //index of the vertex in the vertexPositions array.
    HalfEdge *halfEdge;

    Vertex(int ID)
    {
        id = ID;
        halfEdge = nullptr;
    }

    Vertex()
    {
        id = -1;
        halfEdge = nullptr;
    }
};

class Face 
{
    public:
    HalfEdge *halfEdge;
    int num_sides;

    Face(int sides)
    {
        num_sides = sides;
        halfEdge = nullptr; //always initialize with nullptr.
    }
    Face()
    {
        num_sides = 0;
        halfEdge = nullptr;
    }
};


class Mesh 
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
    void recomputeVertexNormals();
    void clear(); //refreshes the mesh by clearing all the vectors.
    // mesh(int total_verts);
};

Mesh createGrid(int m, int n);
Mesh generateSphere(int m, int n);
Mesh generateCube(int m, int n, int o);
Mesh loadOBJ(const std::string& filename);
void getMeshFromVerts(Mesh &m, std::vector<vec3> &vertexPositions, std::vector<std::vector<int>> &faces, std::vector<vec3> normals = {});
void umbrellaSmooth(Mesh &m, float lambda, int iterations = 1);
HalfEdge *prev(HalfEdge *he);
// void getMeshFromVerts(Mesh &m, vector<vec3> &vertexPositions, vector<vector<int>> &faces);
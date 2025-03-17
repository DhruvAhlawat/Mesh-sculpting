
#include <glm/glm.hpp>
#include <vector>
#include <cmath>
#include <unordered_map>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace glm;

template <class T>
class vectorMap //class to mimic vector while retaining the pointers consistently.
{
    private:
    std::map<int, T> data;
    int tail = 0; //the last unallocated index. 

    public:
    vectorMap()
    {
        data.clear();
        tail = 0;
    }
    vectorMap(int size)
    {
        data.clear();
        tail = size; //size is now the last unallocated index.
    }

    T& push_back(T val)
    {
        data[tail] = val;
        tail++;
        return data[tail]; //also returning its reference. 
    }

    T& operator[](int index)
    {
        if(index > tail) {throw std::out_of_range("Index out of range. Accessed: index");}
        if(data.count(index) == 0) data[index] = T(); //create a default object if it doesn't exist
        return data[index];
    }

    int size()
    {
        return tail;
    }

    void clear()
    {
        data.clear();
        tail = 0;
    }

    void reserve(int size)
    {
        data.reserve(size);
    }

    void resize(int size)
    {
        tail = size;
    }

    vectorMap& operator=(const vectorMap& other) {
        if (this != &other) {
            data = other.data;
            tail = other.tail;
        }
        return *this;
    }

    bool operator==(const vectorMap& other) const {
        return data == other.data && tail == other.tail;
    }

    bool operator!=(const vectorMap& other) const {
        // return !(*this == other);
        return !(*this == other);
    }

};

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

    vectorMap<HalfEdge> halfEdges; // cant keep these as vectors as they have references and pointers to each others which would be messed up on r
    vectorMap<Vertex> verts;
    vectorMap<Face> faces;
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
HalfEdge *prev(HalfEdge *he);   // returns the previous half-edge in the face
void catmullClarkSubdivision(Mesh &m);
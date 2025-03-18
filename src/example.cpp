#include "viewer.hpp"
#include "mesh.hpp"

namespace V = COL781::Viewer;
using namespace glm;

int main() {

   
    V::Viewer v;
    if (!v.initialize("Mesh viewer", 640, 480)) {
        return EXIT_FAILURE;
    }
    // mesh sq = createGrid(40, 40);
    // Mesh sq = generateSphere(10, 10);
    Mesh sq = generateCube(3,3,3);
    // Mesh sq = loadOBJ("meshes/cube.obj");
    // Mesh sq = loadOBJ("meshes/spot_control_mesh.obj");
    
    // all center faces
    vector<int> frontfaces;
    frontfaces.push_back(sq.nearestFaceId(vec3(-0.49,0,0)));
    frontfaces.push_back(sq.nearestFaceId(vec3(0.51,0,0)));
    // frontfaces.push_back(sq.nearestFaceId(vec3(0,0.49,0)));
    // frontfaces.push_back(sq.nearestFaceId(vec3(0,-0.49,0)));
    // frontfaces.push_back(sq.nearestFaceId(vec3(0,0,0.49)));
    // frontfaces.push_back(sq.nearestFaceId(vec3(0,0,-0.49)));
    extrudeMultipleFaces(sq, 0.1, frontfaces);

    sq.triangulateMesh(); // required for rendering
    if(sq.normals.size() == 0)
        sq.recomputeVertexNormals();
    std::cout << "Mesh loaded" << std::endl;
    std::cout << "Vertices: " << sq.vertexPositions.size() << std::endl;
    std::cout << "Triangles: " << sq.triangles.size() << std::endl;
    std::cout << "Faces: " << sq.faces.size() << std::endl;
    std::cout << "Edges: " << sq.edges.size() << std::endl;
    std::cout << "Normals: " << sq.normals.size() << std::endl;
    v.setMesh(sq.vertexPositions.size(), sq.triangles.size(), sq.edges.size(), &sq.vertexPositions[0], &sq.triangles[0], &sq.edges[0], &sq.normals[0]);
    // v.setMesh(20, 10, 12, vertices, triangles, edges, normals);
    v.view();
}

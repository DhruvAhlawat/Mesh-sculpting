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

    int subdiv = 5;
    // Mesh sq = generateCube(subdiv, subdiv, subdiv);
    
    Mesh sq = generateSphere(20, 10);

    vector<int> topHandle = sq.nearestFacesToPoint(vec3(0,1,0), 40);
    extrudeMultipleFaces(sq, 10, topHandle, vec3(0,1,0));


    // Mesh sq = loadOBJ("meshes/cube.obj");
    // Mesh sq = loadOBJ("meshes/spot_control_mesh.obj");
    vector<int> frontfaces;

    frontfaces.push_back(sq.nearestFaceId(vec3(0,0,1)));
    frontfaces.push_back(sq.nearestFaceId(vec3(0,0,-1)));
    frontfaces.push_back(sq.nearestFaceId(vec3(1,0,0)));
    frontfaces.push_back(sq.nearestFaceId(vec3(-1,0,0)));

    for(int i = 1; i < subdiv/2; i++)
    {
        // now to get the faces it is pretty simple. we start off with the center face at 0.
        float pos = 2.0f/subdiv;
        vec3 p1 = vec3(pos,0,1.2);
        vec3 p2 = vec3(1.2,0,pos);

        vector<vec3> points;
        for(int i = -1; i < 2; i+= 2)
        {
            for(int j = -1; j < 2; j+= 2)
            {
                frontfaces.push_back(sq.nearestFaceId(vec3(i*pos,0,j*1.2)));
                frontfaces.push_back(sq.nearestFaceId(vec3(i*1.2,0,j*pos)));
            }
        }
        
        // frontfaces.push_back(sq.nearestFaceId(p2));
        // frontfaces.push_back(sq.nearestFaceId(negp2));
        // frontfaces.push_back(sq.nearestFaceId(p1));
        // frontfaces.push_back(sq.nearestFaceId(negp1));
        // frontfaces.push_back(sq.nearestFaceId(p2));

        // cout << "Face at : " << p2.x << ", " << p2y << ", " << p2.z << " is " << sq.nearestFaceId(p1) << endl;
    }
    vec3 upward = vec3(0,0.1,0);
    for(auto f: frontfaces)
    {
        revolveAndExtrude(sq, f, 7.5, 0.4, 43, upward);
    }
// 
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

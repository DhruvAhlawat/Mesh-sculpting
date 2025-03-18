#include "viewer.hpp"
#include "mesh.hpp"

namespace V = COL781::Viewer;
using namespace glm;

int main() {

    vec3 vertices[] = {
        vec3(-0.5, -0.5, -0.5),
        vec3( 0.5, -0.5, -0.5),
        vec3(-0.5,  0.5, -0.5),
        vec3( 0.5,  0.5, -0.5),

        vec3(-0.5, -0.5,  0.5),
        vec3( 0.5, -0.5,  0.5),
        vec3(-0.5,  0.5,  0.5),
        vec3( 0.5,  0.5,  0.5),

        vec3(-0.5, -0.5, -0.5),
        vec3(-0.5, -0.5,  0.5),
        vec3(-0.5,  0.5, -0.5),
        vec3(-0.5,  0.5,  0.5),

        vec3( 0.5, -0.5, -0.5),
        vec3( 0.5, -0.5,  0.5),
        vec3( 0.5,  0.5, -0.5),
        vec3( 0.5,  0.5,  0.5),

        vec3(-0.5, -0.5, -0.5),
        vec3( 0.5, -0.5, -0.5),
        vec3(-0.5, -0.5,  0.5),
        vec3( 0.5, -0.5,  0.5),
    };

    vec3 normals[] = {
        vec3(0.0, 0.0, -1.0),
        vec3(0.0, 0.0, -1.0),
        vec3(0.0, 0.0, -1.0),
        vec3(0.0, 0.0, -1.0),

        vec3(0.0, 0.0, 1.0),
        vec3(0.0, 0.0, 1.0),
        vec3(0.0, 0.0, 1.0),
        vec3(0.0, 0.0, 1.0),

        vec3(-1.0, 0.0, 0.0),
        vec3(-1.0, 0.0, 0.0),
        vec3(-1.0, 0.0, 0.0),
        vec3(-1.0, 0.0, 0.0),

        vec3(1.0, 0.0, 0.0),
        vec3(1.0, 0.0, 0.0),
        vec3(1.0, 0.0, 0.0),
        vec3(1.0, 0.0, 0.0),

        vec3(0.0, -1.0, 0.0),
        vec3(0.0, -1.0, 0.0),
        vec3(0.0, -1.0, 0.0),
        vec3(0.0, -1.0, 0.0),
    };

    ivec3 triangles[] = {
        ivec3(0, 1, 2),
        ivec3(1, 2, 3),

        ivec3(4, 5, 6),
        ivec3(5, 6, 7),

        ivec3(8, 9, 10),
        ivec3(9, 10, 11),

        ivec3(12, 13, 14),
        ivec3(13, 14, 15),

        ivec3(16, 17, 18),
        ivec3(17, 18, 19)
    };

    ivec2 edges[] = {
        ivec2(0, 1),
        ivec2(1, 3),
        ivec2(3, 2),
        ivec2(2, 0),

        ivec2(4, 5),
        ivec2(5, 7),
        ivec2(7, 6),
        ivec2(6, 4),

        ivec2(0, 4),
        ivec2(1, 5),
        ivec2(2, 6),
        ivec2(3, 7)
    };

    V::Viewer v;
    if (!v.initialize("Mesh viewer", 640, 480)) {
        return EXIT_FAILURE;
    }
    // mesh sq = createGrid(40, 40);
    // Mesh sq = generateSphere(10, 10);
    // Mesh sq = generateCube(2,2,2);
    // Mesh sq = loadOBJ("meshes/cube.obj");
    Mesh sq = loadOBJ("meshes/spot_control_mesh.obj");
    
    for (int i = 0; i < 3; i++) {
    catmullClarkSubdivision(sq);
    sq.triangulateMesh(); // required for rendering
    if(sq.normals.size() == 0)
        sq.recomputeVertexNormals();
    }
    std::cout << "Mesh loaded" << std::endl;
    std::cout << "Vertices: " << sq.vertexPositions.size() << std::endl;
    std::cout << "Triangles: " << sq.triangles.size() << std::endl;
    std::cout << "Edges: " << sq.edges.size() << std::endl;
    std::cout << "Normals: " << sq.normals.size() << std::endl;
    v.setMesh(sq.vertexPositions.size(), sq.triangles.size(), sq.edges.size(), &sq.vertexPositions[0], &sq.triangles[0], &sq.edges[0], &sq.normals[0]);
    // v.setMesh(20, 10, 12, vertices, triangles, edges, normals);
    v.view();
}

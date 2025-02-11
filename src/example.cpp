#include "viewer.hpp"

namespace V = COL781::Viewer;
using namespace glm;

int main() {
    // vec3 vertices[] = {
    //     vec3(-0.5, -0.5, 0.0),
    //     vec3( 0.5, -0.5, 0.0),
    //     vec3(-0.5,  0.5, 0.0),
    //     vec3( 0.5,  0.5, 0.0)
    // };
    // vec3 normals[] = {
    //     vec3(0.0, 0.0, 1.0),
    //     vec3(0.0, 0.0, 1.0),
    //     vec3(0.0, 0.0, 1.0),
    //     vec3(0.0, 0.0, 1.0)
    // };
    // ivec3 triangles[] = {
    //     ivec3(0, 1, 2),
    //     ivec3(1, 2, 3)
    // };
    // ivec2 edges[] = {
    //     ivec2(0, 1),
    //     ivec2(1, 3),
    //     ivec2(3, 2),
    //     ivec2(2, 0)
    // };

    vec3 vertices[] = {
        vec3(-0.5, -0.5, 0.0),
        vec3( 0.5, -0.5, 0.0),
        vec3(-0.5,  0.5, 0.0),
        vec3( 0.5,  0.5, 0.0),
        vec3(-1.5, -0.5, 0.0),
        vec3(-1.5,  0.5, 0.0),
        vec3(1.5,  0.5, 0.0),
        vec3(1.5, -0.5, 0.0)
    };

    vec3 normals[] = {
        vec3(0.0, 0.0, 1.0),
        vec3(0.0, 0.0, 1.0),
        vec3(0.0, 0.0, 1.0),
        vec3(0.0, 0.0, 1.0),
        vec3(0.0, 0.0, 1.0),
        vec3(0.0, 0.0, 1.0),
        vec3(0.0, 0.0, 1.0),
        vec3(0.0, 0.0, 1.0)
    };

    ivec3 triangles[] = {
        ivec3(0, 1, 2),
        ivec3(1, 2, 3),
        ivec3(4, 5, 2),
        ivec3(4, 0, 2),
        ivec3(1, 6, 7),
        ivec3(1, 6, 3)
    };

    ivec2 edges[] = {
        ivec2(0, 1),
        ivec2(1, 3),
        ivec2(3, 2),
        ivec2(2, 0),
        ivec2(4, 5),
        ivec2(6, 7),
        ivec2(4, 0),
        ivec2(5, 2),
        ivec2(6, 3),
        ivec2(7, 1)
    };

    V::Viewer v;
    if (!v.initialize("Mesh viewer", 640, 480)) {
        return EXIT_FAILURE;
    }
    v.setMesh(8, 6, 10, vertices, normals, triangles, edges);
    v.view();
}

#include "mesh.hpp"
#include <cmath>
using namespace std;

void mesh::triangulateMesh()
{
    //will require us to create more edges and hence more faces as well. so the current faces and edges array will have to be completely refreshed.
    vector<ivec3> newTriangles;
    vector<ivec2> newEdges; //this will get updated but. 
    vector<HalfEdge> newHalfEdges; 

    // we iterate over all the existing faces in our mesh.
    for(int i = 0; i < this->faces.size(); i++)
    {
        //else we will break it up into smaller pieces using the first vertex as the common vertex. 
        // triangle fanning procedure.
        HalfEdge *he = faces[i].halfEdge;
        int commonVert = he->head->id;
        he = he->next;
        //we will need to cover 2 faces at once and assign the newly created halfedges as their respective pairs.
        //iterate over the vertices of the face and create new triangles as well as halfedges.
        for(int j = 0; j < faces[i].num_sides - 2; j++)
        {
            ivec3 newFace;
            newFace[0] = commonVert; newFace[1] = he->head->id; 
            he = he->next;
            newFace[2] = he->head->id;
            newTriangles.push_back(newFace);
        }
    }
    this->triangles = newTriangles; //updates the triangles array for us.
    //dont really need to update the edges array as the edges will remain the same on displaying
}

// mesh::mesh(int total_verts)
// {
//     vertexPositions = vector
// }


vec3 spherePos(float theta, float phi)
{
    return vec3(cos(theta) * cos(phi), sin(theta) * cos(phi), sin(phi));
}


void createTriangle(Vertex v1, Vertex v2, Vertex v3)
{

}

mesh getSphere(int m, int n)
{
    //assumes that n >= 2 and m >= 3.
    int total_verts = (n-1)*m + 2; //since the poles are single vertex.
    mesh sphere;
    //create the vertices of the sphere.
    //create the angles phi and theta that we will be using. 
    vector<float> theta(m), phi(n+1);
    for(int i = 0; i < m; i++)
    {
        theta[i] = 2 * M_PI * i / m;
    }
    for(int i = 0; i < n+1; i++)
    {
        phi[i] = M_PI * (i / n - 0.5);
    }   

    //create the pole vertex first. 
    sphere.vertexPositions.push_back(spherePos(theta[0], phi[0])); //the bottom pole vertex.
    sphere.verts.push_back(Vertex(0)); //the pole vertex is thus created. 

    //now we need to connect the next layer with this pole vertex. 
    //at phi[1]
    sphere.vertexPositions.push_back(spherePos(theta[0], phi[1]));
    sphere.verts.push_back(Vertex(1)); 


    sphere.halfEdges.push_back(HalfEdge());
    sphere.halfEdges[0].head = &sphere.verts[0]; //EZ.
    sphere.verts[0].halfEdge = &sphere.halfEdges[0]; sphere.verts[1].halfEdge = &sphere.halfEdges[0]; //pointing both of the vertices to this halfEdge. 
    
    
    
    for(int i = 1; i < m; i++)
    {
        sphere.vertexPositions.push_back(spherePos(theta[i], phi[1]));
        sphere.verts.push_back(Vertex(i+1));

    }




}



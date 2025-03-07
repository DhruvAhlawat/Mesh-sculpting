#include "mesh.hpp"
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

mesh createGrid(int m, int n) {
    mesh sq;
    float dx = 1.0f / n;
    float dy = 1.0f / m;

    sq.vertexPositions.resize((m + 1) * (n + 1));
    sq.verts.resize((m + 1) * (n + 1));

    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= m; j++) {
            int idx = i * (m + 1) + j;
            sq.vertexPositions[idx] = vec3(i * dx, j * dy, 0.0f);
            sq.verts[idx].id = idx;
            sq.verts[idx].halfEdge = nullptr; // No half-edge assigned yet
        }
    }


    // Step 2: Create faces and half-edges
    sq.faces.resize(m * n);
    std::vector<HalfEdge *> edgePointers(4 * m * n, nullptr);

    int edgeIndex = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            int v0 = i * (m + 1) + j;
            int v1 = v0 + 1;
            int v2 = v1 + (m + 1);
            int v3 = v0 + (m + 1);
            
            sq.edges.push_back(ivec2(v0, v1));
            sq.edges.push_back(ivec2(v0, v3));
            if (v2%(m+1) != 0) {sq.edges.push_back(ivec2(v2, v1));}
            if(v2 > (m+1)*n-1) {sq.edges.push_back(ivec2(v2, v3));}

            int fIdx = i * m + j;
            sq.faces[fIdx].num_sides = 4;

            // Create 4 half-edges
            HalfEdge *e0 = new HalfEdge(), *e1 = new HalfEdge(), *e2 = new HalfEdge(), *e3 = new HalfEdge();
            e0->next = e3; e1->next = e0; e2->next = e1; e3->next = e2;
            e0->head = &sq.verts[v0]; 
            e1->head = &sq.verts[v1];
            e2->head = &sq.verts[v2];
            e3->head = &sq.verts[v3];
            e0->face = e1->face = e2->face = e3->face = &sq.faces[fIdx];

            sq.faces[fIdx].halfEdge = e0;

            // Store half-edges for twin assignment
            edgePointers[edgeIndex] = e0;
            edgePointers[edgeIndex + 1] = e1;
            edgePointers[edgeIndex + 2] = e2;
            edgePointers[edgeIndex + 3] = e3;
            edgeIndex += 4;

            // Store half-edges in mesh
            sq.halfEdges.push_back(*e0);
            sq.halfEdges.push_back(*e1);
            sq.halfEdges.push_back(*e2);
            sq.halfEdges.push_back(*e3);
        }
    }

    // Step 3: Assign twins
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            int fIdx = i * m + j;
            int e0 = fIdx * 4, e1 = e0 + 1, e2 = e1 + 1, e3 = e2 + 1;

            // Right neighbor (e1 -> e3 of right face)
            if (j < m - 1) {
                int rightFace = fIdx + 1;
                edgePointers[e1]->pair = edgePointers[rightFace * 4 + 3];
                edgePointers[rightFace * 4 + 3]->pair = edgePointers[e1];
            }

            // Top neighbor (e2 -> e0 of the top face)
            if (i < n - 1) {
                int topFace = fIdx + m;
                edgePointers[e2]->pair = edgePointers[topFace * 4];
                edgePointers[topFace * 4]->pair = edgePointers[e2];
            }
        }
    }
    return sq;
}

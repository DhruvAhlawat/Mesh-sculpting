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




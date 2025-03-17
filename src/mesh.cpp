#include "mesh.hpp"
using namespace std;

void Mesh::triangulateMesh()
{
    //will require us to create more edges and hence more faces as well. so the current faces and edges array will have to be completely refreshed.
    vector<ivec3> newTriangles;
    // we iterate over all the existing faces in our mesh.
    for(int i = 0; i < this->faces.size(); i++)
    {
        //else we will break it up into smaller pieces using the first vertex as the common vertex. 
        // triangle fanning procedure.
        HalfEdge *he = faces[i].halfEdge;
        if(he == nullptr)
        {
            cout << "NULLPTR ENCOUNTERED for face: " << i << endl;
        }
        int commonVert = he->head->id;
        he = he->next;
        //we will need to cover 2 faces at once and assign the newly created halfedges as their respective pairs.
        //iterate over the vertices of the face and create new triangles 

        for(int j = 0; j < faces[i].num_sides - 2; j++)
        {
            ivec3 newFace;
            newFace[0] = commonVert; 
            newFace[1] = he->head->id; 
            he = he->next;
            newFace[2] = he->head->id;
            newTriangles.push_back(newFace);
        }
    }
    this->triangles = std::move(newTriangles); //updates the triangles array for us.
    //dont really need to update the edges array as the edges will remain the same on displaying
}

void Mesh::clear()
{
    vertexPositions.clear();
    triangles.clear();
    normals.clear();
    edges.clear();
    halfEdges.clear();
    verts.clear();
    faces.clear();
}

vec3 spherePos(float theta, float phi)
{
    return vec3(cos(theta) * cos(phi), sin(theta) * cos(phi), sin(phi));
}

Mesh createGrid(int m, int n) {
    Mesh sq;
    float dx = 1.0f / n;
    float dy = 1.0f / m;

    sq.vertexPositions.resize((m + 1) * (n + 1));
    sq.verts.resize((m + 1) * (n + 1));

    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= m; j++) {
            int idx = i * (m + 1) + j;
            sq.vertexPositions[idx] = vec3(i * dx, j * dy, 0.0f);
            sq.verts[idx].id = idx;
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

Mesh generateSphere(int m, int n) {
    Mesh sphereMesh;

    // Generate vertices
    for (int j = 1; j < n; j++) { // Exclude poles
        float phi = M_PI * j / n;
        for (int i = 0; i < m; i++) {
            float theta = 2.0f * M_PI * i / m;
            float x = cos(theta) * sin(phi);
            float y = sin(theta) * sin(phi);
            float z = cos(phi);
            sphereMesh.vertexPositions.emplace_back(x, y, z);
        }
    }

    // Add poles
    int northPoleIndex = sphereMesh.vertexPositions.size();
    sphereMesh.vertexPositions.emplace_back(0.0f, 0.0f, 1.0f);
    int southPoleIndex = sphereMesh.vertexPositions.size();
    sphereMesh.vertexPositions.emplace_back(0.0f, 0.0f, -1.0f);

    vector<vector<int>> faces;
    
    // Middle quads (excluding poles)
    for (int j = 0; j < n - 2; j++) { // stacks
        for (int i = 0; i < m; i++) { // slices
            int nextI = (i + 1) % m;
            int currRow = j * m;
            int nextRow = (j + 1) * m;

            // quad face, anticlockwise
            faces.push_back({
                currRow + i,
                nextRow + i,
                nextRow + nextI,
                currRow + nextI
            });
        }
    }

    // Top cap (fan around north pole)
    for (int i = 0; i < m; i++) {
        int nextI = (i + 1) % m;
        faces.push_back({
            northPoleIndex,
            i,
            nextI
        });
    }

    // Bottom cap (fan around south pole)
    int bottomStart = (n - 2) * m;
    for (int i = 0; i < m; i++) {
        int nextI = (i + 1) % m;
        faces.push_back({
            southPoleIndex,
            bottomStart + nextI,
            bottomStart + i
        });
    }
    auto vertpos = sphereMesh.vertexPositions;
    getMeshFromVerts(sphereMesh, vertpos, faces);

    return sphereMesh;
}


Mesh generateCube(int m, int n, int o) {
    Mesh cubeMesh;
    std::vector<std::vector<std::vector<int>>> vertexIndices(m + 1, std::vector<std::vector<int>>(n + 1, std::vector<int>(o + 1, -1)));
    std::vector<std::vector<int>> faces; // New faces vector

    // Generate vertices
    int index = 0;
    for (int i = 0; i <= m; i++) {
        float x = -0.5f + i / float(m);
        for (int j = 0; j <= n; j++) {
            float y = -0.5f + j / float(n);
            for (int k = 0; k <= o; k++) {
                float z = -0.5f + k / float(o);
                cubeMesh.vertexPositions.emplace_back(x, y, z);
                vertexIndices[i][j][k] = index++;
            }
        }
    }

    // Generate faces and edges
    auto addQuad = [&](int v0, int v1, int v2, int v3) {
        faces.push_back({v0, v1, v2, v3});
    };

    // Generate faces
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            addQuad(vertexIndices[i][j][0], vertexIndices[i][j + 1][0], vertexIndices[i + 1][j + 1][0], vertexIndices[i + 1][j][0]); // +Z face
            addQuad(vertexIndices[i][j][o], vertexIndices[i + 1][j][o], vertexIndices[i + 1][j + 1][o], vertexIndices[i][j + 1][o]); // -Z face
        }
    }
    for (int i = 0; i < m; i++) {
        for (int k = 0; k < o; k++) {
            addQuad(vertexIndices[i][n][k], vertexIndices[i][n][k + 1], vertexIndices[i + 1][n][k + 1], vertexIndices[i + 1][n][k]); // +Y face
            addQuad(vertexIndices[i][0][k], vertexIndices[i + 1][0][k], vertexIndices[i + 1][0][k + 1], vertexIndices[i][0][k + 1]); // -Y face
        }
    }
    for (int j = 0; j < n; j++) {
        for (int k = 0; k < o; k++) {
            addQuad(vertexIndices[0][j][k], vertexIndices[0][j][k + 1], vertexIndices[0][j + 1][k + 1], vertexIndices[0][j + 1][k]); // +X face
            addQuad(vertexIndices[m][j][k], vertexIndices[m][j + 1][k], vertexIndices[m][j + 1][k + 1], vertexIndices[m][j][k + 1]); // -X face
        }
    }
    auto vertpos = cubeMesh.vertexPositions;
    getMeshFromVerts(cubeMesh, vertpos, faces);

    return cubeMesh;
}


Mesh loadOBJ(const std::string& filename) {
    Mesh mesh;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open OBJ file: " << filename << std::endl;
        return mesh;
    }

    std::vector<glm::vec3> tempNormals; //why temp?
    std::vector<glm::vec3> vertexPositions;
    std::vector<std::vector<int>> faces;

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string prefix;
        iss >> prefix;
        if (prefix == "v") {
            // Read vertex position
            glm::vec3 vertex;
            iss >> vertex.x >> vertex.y >> vertex.z;
            vertexPositions.push_back(vertex);

        } else if (prefix == "vn") {
            // Read vertex normal
            glm::vec3 normal;
            iss >> normal.x >> normal.y >> normal.z;
            tempNormals.push_back(normal);

        } else if (prefix == "f") {
            // Read face (triangles or n-gons)
            std::vector<int> faceIndices;
            std::string vertInfo;
            while (iss >> vertInfo) {
                std::istringstream vertStream(vertInfo);
                std::string vIndex;
                std::getline(vertStream, vIndex, '/'); // Extract vertex index (ignore textures/normals)
                int vertexIdx = std::stoi(vIndex) - 1; // Convert to 0-based indexing
                faceIndices.push_back(vertexIdx);
            }
            faces.push_back(faceIndices); 

            // // NO Triangulation in the mesh. Only for rendering.
            // // Triangulate n-gon into a triangle fan
            // for (size_t i = 1; i < faceIndices.size() - 1; i++) {
            //     mesh.triangles.emplace_back(faceIndices[0], faceIndices[i], faceIndices[i + 1]);
            // }

            // // Store original polygon edges for visualization
            // for (size_t i = 0; i < faceIndices.size(); i++) {
            //     mesh.edges.emplace_back(faceIndices[i], faceIndices[(i + 1) % faceIndices.size()]);
            // }
        }
    }

    file.close();

    //the below function also clears the mesh object first. So dont update anything in the mesh before this is called.
    getMeshFromVerts(mesh, vertexPositions, faces, tempNormals); //this will create the halfedges and faces as well.
    return mesh;
}

void Mesh::recomputeVertexNormals() 
{
    cout << "recomputing vertex normals" << endl;
    if(triangles.size() == 0)
    {
        cout << "triangulating mesh" << endl;
        this->triangulateMesh();
    }
    normals.assign(vertexPositions.size(), glm::vec3(0.0f));
    // Compute face normals and accumulate into vertex normals
    for (const auto& tri : triangles) {
        glm::vec3 v0 = vertexPositions[tri.x];
        glm::vec3 v1 = vertexPositions[tri.y];
        glm::vec3 v2 = vertexPositions[tri.z];

        // Compute face normal using cross product
        glm::vec3 normal = (glm::cross(v1 - v0, v2 - v0));
        // Accumulate normal weighted by face area
        float area = glm::length(glm::cross(v1 - v0, v2 - v0)) * 0.5f;
        normals[tri.x] += normal * area; 
        normals[tri.y] += normal * area; 
        normals[tri.z] += normal * area; 
    }

    // Normalize all vertex normals
    for (auto& normal : normals) {
        normal = glm::normalize(normal);
    }
}

// void getMeshFromVerts(Mesh &m, vector<vec3> &vertexPositions, vector<vector<int>> &faces, vector<vec3> &normals)
// {
//     getMeshFromVerts(m, vertexPositions, faces);
//     m.normals = normals;

// }

void getMeshFromVerts(Mesh &m, vector<vec3> &vertexPositions, vector<vector<int>> &faces, vector<vec3> normals)
{
    m.clear();
    m.vertexPositions = vertexPositions;
    m.normals = normals;
    for(int i = 0; i < vertexPositions.size(); i++)
    {
        m.verts.push_back(Vertex(i)); //creates a vertex object with the ith id.
    }  // no pointers to consider yet so can use emplaceback. 

    //now we will form an edge mapping. 
    //from this we will form the halfedges and faces. this map will store the index of the halfedge between the two vertices.
    map<pair<int,int>, int> halfEdgeMap;
    //the map stores the location of the half edge corresponding to these vertices in the direction of first to second.
    set<pair<int,int>> edges;
    m.faces = vectorMap<Face>(faces.size());
    int total_halfedges = 0;

    //we will iterate over the faces and create the halfedges and faces and only assign the face and vertex pointers.
    for(int i = 0; i < faces.size(); i++)
    {
        //now connect all the vertices in the face.
        vector<int> face = faces[i];        
        m.faces[i].num_sides = face.size(); 
        face.push_back(face[0]); //to close the loop, to make it easier in the following loops.
        
        //loop that initializes all the half edges for this face (and their pairs).
        for(int j = 0; j < face.size() - 1; j++)
        {
            int v1 = face[j], v2 = face[j+1];
            //we will make sure that we always go in the correct order 
            if(halfEdgeMap.count({v1,v2}) == 0) //if doesnt exist we create it.
            {
                //create a new halfedge here.
                //only apply the vert and faces pointers for now.
                m.halfEdges.push_back(HalfEdge(nullptr, nullptr, &m.verts[v2], &m.faces[i]));
                halfEdgeMap[{v1,v2}] = m.halfEdges.size() - 1;
            }
            m.halfEdges[halfEdgeMap[{v1,v2}]].face = &m.faces[i];
            m.halfEdges[halfEdgeMap[{v1,v2}]].head = &m.verts[v2];
            
            //we also create its twin if it doesn't exist yet. With no face initially.
            if(halfEdgeMap.count({v2,v1}) == 0)
            {
                //we create it if it doesn't exist. EZ
                m.halfEdges.push_back(HalfEdge(nullptr, nullptr, &m.verts[v1], nullptr));
                halfEdgeMap[{v2,v1}] = m.halfEdges.size() - 1;
            }  
            
            //finally update the set of edges.
            edges.insert({std::min(v1,v2), std::max(v1,v2)});
        }
    }

    //after this we have the halfEdges, verts, faces vector all ready without any further size changes.
    for(int i = 0; i < faces.size(); i++)
    {
        //now we can assign the next and pair pointers for all the halfEdges.
        vector<int> face = faces[i];
        face.push_back(face[0]); face.push_back(face[1]); 

        for(int j = 0; j < face.size() - 2; j++)
        {
            int v1 = face[j], v2 = face[j+1];
            int curHalfEdge = halfEdgeMap[{v1,v2}];
            int nextHalfEdge = halfEdgeMap[{v2,face[j+2]}];
            int pairedHalfEdge = halfEdgeMap[{v2, v1}];
            m.halfEdges[curHalfEdge].next = &m.halfEdges[nextHalfEdge];
            m.halfEdges[curHalfEdge].pair = &m.halfEdges[pairedHalfEdge];
            m.halfEdges[pairedHalfEdge].pair = &m.halfEdges[curHalfEdge];
            m.halfEdges[curHalfEdge].face = &m.faces[i]; //just for safety reassigning it.

            m.verts[v1].halfEdge = &m.halfEdges[curHalfEdge]; //v1 points to a half edge that points away from it.
        }

        //now we also need to assign an HalfEdge to each face. so we will do it for the first 2 vertex halfedge. EZ
        m.faces[i].halfEdge = &m.halfEdges[halfEdgeMap[{face[0], face[1]}]];
    }
    
    //finally set the edges vector as well.
    for(auto edge : edges)
    {
        m.edges.push_back(ivec2(edge.first, edge.second));
    }
    
}

HalfEdge *prev(HalfEdge *he)
{
    HalfEdge *prev = he;
    if(he->face == nullptr)
    {
        return nullptr; //if there is no face associated then the next pointer will also be null. But should it be?
    }
    while(prev->next != he)
    {
        prev = prev->next;
    }
    return prev;
}

vector<int> getVertexNeighbours(Mesh &m, int vid)
{
    vector<int> neighbours;
    HalfEdge *he = m.verts[vid].halfEdge;
    do
    {
        neighbours.push_back(he->head->id);
        he = he->pair->next;
    } while (he != nullptr && he != m.verts[vid].halfEdge);

    //now incase of vertices present at a boundary of a mesh. The above will not always capture all the neighbours.
    //so we need to go check the previous vertices as well.
    he = prev(m.verts[vid].halfEdge); 
    if(he == nullptr) return neighbours;
    //first we check if we have a boundary vertex by checking for the presence of this vertex in the neighbours list at the end.
    if(neighbours[neighbours.size() - 1] == he->pair->head->id)
    {
        //if this vertex is already added then we are done.
        return neighbours;
    }

    //otherwise we add all the vertices present in the opposite order frmo this halfedge as well.
    he = he->pair;
    do
    {
       neighbours.push_back(he->head->id);
       he = prev(he); if(he == nullptr) break;
       he = he->pair; 
    } while (he != nullptr);
    return neighbours;
}

void UmbrellaSmooth(Mesh &m, float lambda, int iterations)
{
    for(int iter = 0; iter < iterations; iter++)
    {
        vector<vec3> deltas(m.verts.size(), vec3(0.0f));
        for(int i = 0; i < m.verts.size(); i++)
        {
            Vertex v = m.verts[i];
            vector<int> neighbours = getVertexNeighbours(m, i);
            for(int j = 0; j < neighbours.size(); j++)
            {
                deltas[i] += (m.vertexPositions[neighbours[j]] - m.vertexPositions[v.id]);
            }
            deltas[i] /= neighbours.size();
        }

        //then we update the positions by lambda * delta
        for(int i = 0; i < m.verts.size(); i++)
        {
            m.vertexPositions[i] += lambda * deltas[i];
        }
    }
}

void extrude(Mesh &m, float offset, int faceid, vec3 direction, Face *f)
{
    if(f == nullptr)
    {
        f = &m.faces[faceid];
    }
    HalfEdge *he = f->halfEdge;

    //instead of doing it in a complicated way dependent on implementation. Lets bruteforce and use a map to remember all the new creations 
    //and then finally join them together using our halfedge knowledge.
    map<pair<int,int>, int> halfEdges;
    vector<int> newVerts;
    vector<int> newFaces;

    m.verts.push_back(Vertex(m.verts.size()));
    m.vertexPositions.push_back(m.vertexPositions[he->pair->head->id]);
    newVerts.push_back(m.verts.size() - 1);

    Vertex *back = he->pair->head;
    Vertex *backDup = &m.verts[m.verts.size() - 1]; //duplicated this beforehand.

    //and we also create a halfedge that points from back to backDup here.
    HalfEdge *backToBackDup = &m.halfEdges.push_back(HalfEdge());
    backToBackDup->head = backDup;

    Vertex *startBack = back, *startBackDup = backDup;
    HalfEdge *startBackToBackDup = backToBackDup;

    do
    {
        Vertex *vOg = he->head;
        HalfEdge *ogPair = he->pair; //the original pair of this halfedge. 
        Vertex *vDup; 
        if(he->next != f->halfEdge)
        {
            m.verts.push_back(Vertex(m.verts.size()));
            m.vertexPositions.push_back(m.vertexPositions[vOg->id]);
            newVerts.push_back(m.verts.size() - 1);
            vDup = &m.verts[m.verts.size() - 1];
        }
        else
        {
            //else we are on the last halfEdge. in this case no duplication is necessary as we already have vertices duplicated originally.
            vDup = startBackDup;
            vOg = startBack; //although he->head should give the same results imo
        }

        m.edges.push_back({std::min(vDup->id, backDup->id), std::max(vDup->id, backDup->id)});
        m.edges.push_back({std::min(vDup->id, vOg->id), std::max(vDup->id, vOg->id)});
        HalfEdge *h1 = &m.halfEdges.push_back(HalfEdge());
        HalfEdge *h2 = &m.halfEdges.push_back(HalfEdge());
        HalfEdge *h3 = &m.halfEdges.push_back(HalfEdge());
        HalfEdge *h4; 
        if(he->next != f->halfEdge)
        { h4 = &m.halfEdges.push_back(HalfEdge()); }
        else
        { h4 = startBackToBackDup; }
        
        //now we will make a new face that goes from backDup to vDup to vOg to backVert to backDup.
        Face *newFace = &m.faces.push_back(4); //a new quad face.
        newFace->halfEdge = h1; //pointing to the top halfedge.
        he->head = vDup; //pointing to the duplicate now.
        he->pair = h1; 
        
        h1->head = backDup; h1->face = newFace; h1->pair = he; h1->next = h2;

        h2->head = back; h2->face = newFace; h2->pair = backToBackDup; h2->next = h3; 
        backToBackDup->pair = h2; //gotta set this up. Can do it later too using a map, but this is possible so its fine.

        h3->head = vOg; h3->face = newFace; h3->pair = ogPair; h3->next = h4;

        h4->head = vDup; h4->face = newFace; h4->next = h1; //h4's pair is not yet set, or in the case of startBackToBackDup, already set beforehand.

        backToBackDup = h4; // to set its pair later in the next iteration.
        back = vOg;
        backDup = vDup;


        he = he->next;
    } while(he != f->halfEdge);
    

    //after this operation we can move it all a certain amount. Also gotta calculate its normal.
    if(direction == vec3(0.0f))
    {
        vec3 originpos = m.vertexPositions[newVerts[0]];
        for(int i = 1; i < newVerts.size() - 1; i++)
        {
            direction += (cross(m.vertexPositions[newVerts[i]] - originpos, m.vertexPositions[newVerts[i+1]] - originpos));
        }
    }
    direction = normalize(direction); //always normalize first.
    for(int i = 0; i < newVerts.size(); i++)
    {
        //update the position of all the vertices.
        m.vertexPositions[newVerts[i]] += offset * direction;
    }
}
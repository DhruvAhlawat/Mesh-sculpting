#include "mesh.hpp"


glm::vec3 randomVec3(float stddev) 
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::normal_distribution<float> dist(0.0f, stddev); // Mean = 0, Stddev = 1

    return glm::vec3(dist(gen), dist(gen), dist(gen));
}

void addNoise(Mesh &m, float threshold)
{
    for(int i = 0; i < m.vertexPositions.size(); i++)
    {
        m.vertexPositions[i] += randomVec3(threshold);
    }
}


// Clamps a value between min and max
float clampVal(float x, float minVal, float maxVal) {
    return std::max(minVal, std::min(maxVal, x));
}

// Computes squared distance from a point to a line segment
float pointToSegmentSq(const vec3& p, const vec3& a, const vec3& b) {
    vec3 ab = b - a;
    vec3 ap = p - a;
    float abLenSq = dot(ab, ab);
    
    if (abLenSq == 0.0f) return dot(ap, ap); // Degenerate segment (same point)

    float t = clampVal(dot(ap, ab) / abLenSq, 0.0f, 1.0f);
    vec3 closest = a + t * ab;
    return dot(p - closest, p - closest);
}

// Computes squared distance from a point to a triangle
float pointToTriSq(const vec3& p, const vec3& a, const vec3& b, const vec3& c) {
    vec3 ab = b - a, ac = c - a;
    vec3 normal = cross(ab, ac);
    float normLenSq = dot(normal, normal);

    if (normLenSq == 0.0f) { 
        // Degenerate triangle (collinear points), return min edge distance
        return std::min(
            pointToSegmentSq(p, a, b),
            std::min(pointToSegmentSq(p, a, c),
            pointToSegmentSq(p, b, c))
        );
    }

    normal = normalize(normal); // Normalize normal

    // Project point onto triangle plane
    float d = dot(normal, a);
    vec3 proj = p - normal * (dot(normal, p) - d);

    // Compute barycentric coordinates
    vec3 v0 = c - a, v1 = b - a, v2 = proj - a;
    float d00 = dot(v0, v0), d01 = dot(v0, v1);
    float d11 = dot(v1, v1), d20 = dot(v2, v0), d21 = dot(v2, v1);
    float denom = d00 * d11 - d01 * d01;

    if (denom == 0.0f) { // Just in case, dont want no errors
        // Degenerate triangle, return edge distances
        return std::min(
            pointToSegmentSq(p, a, b),
            std::min(pointToSegmentSq(p, a, c),
            pointToSegmentSq(p, b, c))
        );
    }

    float v = (d11 * d20 - d01 * d21) / denom;
    float w = (d00 * d21 - d01 * d20) / denom;
    float u = 1.0f - v - w;

    if (u >= 0 && v >= 0 && w >= 0) {
        return dot(p - proj, p - proj); // Inside the triangle
    }

    // Otherwise, return min distance to edges
    return std::min(
        pointToSegmentSq(p, a, b),
        std::min(pointToSegmentSq(p, a, c),
        pointToSegmentSq(p, b, c))
    );
}


vector<int> Mesh::nearestFacesToPoint(const vec3 point, int maxreturn)
{
    if(triangles.size() == 0)
    {
        cout << "triangulating mesh" << endl;
        this->triangulateMesh(); //gotta triangulate first.
    }
    map<int, float> faceDistances; 
    float minDist = std::numeric_limits<float>::max();

    for (size_t i = 0; i < triangles.size(); ++i) {
        ivec3 tri = triangles[i];

        float distSq = pointToTriSq(point, vertexPositions[tri.x], vertexPositions[tri.y], vertexPositions[tri.z]);
        faceDistances[triangle_to_face[i]] = distSq;
    }

    vector<pair<float, int>> sortedFaces;
    for(auto it = faceDistances.begin(); it != faceDistances.end(); it++)
    {
        sortedFaces.push_back(make_pair(it->second, it->first));
    }
    sort(sortedFaces.begin(), sortedFaces.end());
    vector<int> finalans;
    for(int i = 0; i < maxreturn; i++)
    {
        finalans.push_back(sortedFaces[i].second);
    }
    return finalans;
}

// Finds the nearest face ID to a given point
int Mesh::nearestFaceId(const vec3 point) {
    if(triangles.size() == 0)
    {
        cout << "triangulating mesh" << endl;
        this->triangulateMesh(); //gotta triangulate first.
    }
    int closestTri = -1;
    float minDist = std::numeric_limits<float>::max();
    // cout << "trying to "
    for (size_t i = 0; i < triangles.size(); ++i) {
        ivec3 tri = triangles[i];

        float distSq = pointToTriSq(point, vertexPositions[tri.x], vertexPositions[tri.y], vertexPositions[tri.z]);
        if (distSq < minDist) {
            minDist = distSq;
            closestTri = i;
        }
    }

    if (closestTri == -1) return -1; // No valid triangle found, although this should never occur. Butt just in case
    return triangle_to_face[closestTri];
}

void Mesh::triangulateMesh()
{
    this->triangles.clear();
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
            triangle_to_face[newTriangles.size() - 1] = i; //this triangle is a part of this face. 
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

vec3 spherePos(float phi, float theta)
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
            float z = cos(theta) * sin(phi);
            float x = sin(theta) * sin(phi);
            float y = cos(phi);
            sphereMesh.vertexPositions.emplace_back(x, y, z);
        }
    }

    // // Add poles
    // int northPoleIndex = sphereMesh.vertexPositions.size();
    // sphereMesh.vertexPositions.emplace_back(0.0f, 0.0f, 1.0f);
    // int southPoleIndex = sphereMesh.vertexPositions.size();
    // sphereMesh.vertexPositions.emplace_back(0.0f, 0.0f, -1.0f);

    int northPoleIndex = sphereMesh.vertexPositions.size();
    sphereMesh.vertexPositions.emplace_back(0.0f, 1.0f, 0.0f);
    int southPoleIndex = sphereMesh.vertexPositions.size();
    sphereMesh.vertexPositions.emplace_back(0.0f, -1.0f, 0.0f);

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
                if (i==0 || i==m || j==0 || j==n || k==0 || k==o){
                cubeMesh.vertexPositions.emplace_back(x, y, z);
                vertexIndices[i][j][k] = index++;
                }
            }
        }
    }

    // Generate faces and edges
    auto addQuad = [&](int v0, int v1, int v2, int v3) 
    {
        // cout << v0 << ", " << v1  << ", " << v2 << ", " << v3 << endl;
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
    int maxiterations = 1e2; // assuming that a a normal vertex probably doesn't have more than 100 neighbours.
    int iter = 0;
    cout << "   trying to find neighbours of " << vid << endl;
    HalfEdge *he = m.verts[vid].halfEdge;
    // cout << ""
    if(he == nullptr)
    {
        cout << "NULLPTR FOUND HalfEdge for vertex: " << vid << endl;
        return {}; //no neighbours incase of isolated vertex.
    }
    do
    {
        cout << he->head->id << " ";
        neighbours.push_back(he->head->id);
        he = he->pair->next;
        iter++;
        if(iter > maxiterations)
        {
            cout << "\n MAX ITERATIONS EXCEEDED IN FINDING NEIGHBOUR OF VERTEX: " << vid << endl;
            throw "MAX ITERATIONS EXCEEDED IN FINDING NEIGHBOUR OF VERTEX";
        }
        // cout << "going to next" << endl;
    } while (he != nullptr && he != m.verts[vid].halfEdge);
    cout << endl;
    //now incase of vertices present at a boundary of a mesh. The above will not always capture all the neighbours.
    //so we need to go check the previous vertices as well.
    // cout << "trying to get prev" << endl;
    he = prev(m.verts[vid].halfEdge); 
    // cout << "got prev" << endl;
    if(he == nullptr) return neighbours;
    //first we check if we have a boundary vertex by checking for the presence of this vertex in the neighbours list at the end.
    if(neighbours[neighbours.size() - 1] == he->pair->head->id)
    {
        //if this vertex is already added then we are done.
        return neighbours;
    }

    // cout << "prev didnt match. " << endl;
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

void umbrellaSmooth(Mesh &m, float lambda, int iterations)
{
    // getVertexNeighbours(m, 13);
    for(int iter = 0; iter < iterations; iter++)
    {
        vector<vec3> deltas(m.verts.size(), vec3(0.0f));
        for(int i = 0; i < m.verts.size(); i++)
        {
            // cout << "doing vert " << i << endl;
            Vertex v = m.verts[i];
            vector<int> neighbours = getVertexNeighbours(m, i);
            for(int j = 0; j < neighbours.size(); j++)
            {
                deltas.at(i) += (m.vertexPositions.at(neighbours.at(j)) - m.vertexPositions.at(v.id));
            }
            deltas.at(i) /= neighbours.size();
        }

        //then we update the positions by lambda * delta
        for(int i = 0; i < m.verts.size(); i++)
        {
            m.vertexPositions[i] += lambda * deltas[i];
        }
        // cout << "iteration done: " << iter << endl;
    }
}

vec3 faceNormal(Mesh &m, Face *f)
{
    HalfEdge *he = f->halfEdge;
    vector<vec3> vpos;
    vpos.push_back(m.vertexPositions[he->head->id]);

    he = he->next;
    do
    {
        vpos.push_back(m.vertexPositions[he->head->id]);
        he = he->next;
    } while (he != f->halfEdge);
    
    vec3 normal = vec3(0.0f);
    for(int i = 1; i < vpos.size() - 1; i++)
    {
        normal += cross(vpos[i] - vpos[0], vpos[i+1] - vpos[0]);
    }
    normal = normalize(normal);
    return normal;
}

void rotateFace(Mesh &m, int faceid, float theta, vec3 axis)
{
    Face *f = &m.faces[faceid];
    vec3 center = vec3(0.0f);

    HalfEdge *he = f->halfEdge;
    do
    {
        center += m.vertexPositions[he->head->id];
        he = he->next;
    } while (he != f->halfEdge);
    center /= f->num_sides;

    float angle = glm::radians(theta);  // Convert degrees to radians
    glm::mat4 rotationMatrix = glm::rotate(glm::mat4(1.0f), angle, glm::normalize(axis));

    he = f->halfEdge;
    do
    {
        vec3 ogpos = m.vertexPositions[he->head->id];
        m.vertexPositions[he->head->id] = glm::vec3(rotationMatrix * glm::vec4(ogpos - center, 1.0f)) + center;
        he = he->next;
    } while (he != f->halfEdge);

}
vec3 groupNormalFace(Mesh &m, vector<int> faceIds)
{
    vec3 normal = vec3(0.0f);
    for(int i = 0; i < faceIds.size(); i++)
    {
        normal += faceNormal(m, &m.faces[faceIds[i]]);
    }
    normal = normalize(normal);
    return normal;
}

void extrudeMultipleFaces(Mesh &m, float offset, vector<int> faceIds, vec3 direction)
{
    if (direction == vec3(0.0f))
    {
        for(int i =0;i < faceIds.size(); i++)
        {
            direction += faceNormal(m, &m.faces[faceIds[i]]);
        }
    }
    direction = normalize(direction);

    set<int> faceSet(faceIds.begin(), faceIds.end());
    if(faceIds.size() != faceSet.size())
    {
        cout << "WARNING: Duplicate faces in the input. Removing Duplicates" << endl;
        // return;
    }
    faceIds = vector<int>(faceSet.begin(), faceSet.end()); //to remove duplicates


    set<int> selectedFaces(faceIds.begin(), faceIds.end());
    set<int> verts;
    // set<int> ;

    map<pair<int, int>, int> edgeCount; // Track edge occurrences
    map<int, int> duplicatedVertices; // Original vertex -> duplicated vertex
    vector<int> newVerts;
    
    vector<HalfEdge*> boundaryHes;
    vector<Face*> boundaryFaces;


    // Step 1: Identify boundary edges and count occurrences
    for (int fid : faceIds)
    {
        Face *f = &m.faces[fid];
        HalfEdge *he = f->halfEdge;
        do
        {
            verts.insert(he->head->id);
            verts.insert(he->pair->head->id);

            pair<int,int> edgeKey = make_pair(std::min(he->head->id, he->pair->head->id),
                                     std::max(he->head->id, he->pair->head->id));
            edgeCount[edgeKey]++;
            he = he->next;
        } while (he != f->halfEdge);
    }

    for (int fid : faceIds)
    {
        Face *f = &m.faces[fid];
        HalfEdge *he = f->halfEdge;
        do
        {
            pair<int,int> edgeKey = make_pair(std::min(he->head->id, he->pair->head->id),
                                     std::max(he->head->id, he->pair->head->id));
            if (edgeCount[edgeKey] == 1) // Boundary edge. so we add this halfEdge to our boundary list, along with this face.
            {
                boundaryHes.push_back(he);
                boundaryFaces.push_back(f);
            }
            he = he->next;
        } while (he != f->halfEdge);
    }

    map<pair<int,int>, HalfEdge*> heMap; //for the new halfEdges created.
    set<HalfEdge*> boundaryHeSet(boundaryHes.begin(), boundaryHes.end()); //for quick lookup.
    set<pair<int,int>> edgesToDelete; //for the edges that we need to delete after the extrusion. (between boundary and inner faces.)
    vector<ivec2> newEdges;

    // Step 2: Duplicate boundary vertices, and create boundary halfEdges.
    for (const auto &entry : edgeCount)
    {
        if (entry.second == 1) // Boundary edge (appears once in selected faces)
        {
            int v1 = entry.first.first;
            int v2 = entry.first.second;
            if (duplicatedVertices.find(v1) == duplicatedVertices.end())
            {
                m.verts.push_back(Vertex(m.verts.size()));
                m.vertexPositions.push_back(m.vertexPositions[v1]);
                duplicatedVertices[v1] = m.verts.size() - 1;
                newVerts.push_back(m.verts.size() - 1);
                int dup = m.verts.size() - 1;
                // Create new half-edge
                HalfEdge *h = &m.halfEdges.push_back(HalfEdge());
                h->head = &m.verts[v1];
                heMap[make_pair(dup, v1)] = h;
                m.verts[dup].halfEdge = h;

                HalfEdge *h2 = &m.halfEdges.push_back(HalfEdge());
                h2->head = &m.verts[dup];
                heMap[make_pair(v1, dup)] = h2;
                m.verts[v1].halfEdge = h2;

                h->pair = h2; h2->pair = h;
                verts.insert(dup); verts.erase(v1);

                m.edges.push_back({std::min(dup,v1), std::max(dup, v1)});
            }
            if (duplicatedVertices.find(v2) == duplicatedVertices.end())
            {
                m.verts.push_back(Vertex(m.verts.size()));
                m.vertexPositions.push_back(m.vertexPositions[v2]);
                duplicatedVertices[v2] = m.verts.size() - 1;
                newVerts.push_back(m.verts.size() - 1);
                int dup = m.verts.size() - 1;

                HalfEdge *h = &m.halfEdges.push_back(HalfEdge());
                h->head = &m.verts[v2];
                heMap[make_pair(dup, v2)] = h;
                m.verts[dup].halfEdge = h;

                HalfEdge *h2 = &m.halfEdges.push_back(HalfEdge());
                h2->head = &m.verts[dup];
                heMap[make_pair(v2, dup)] = h2;
                m.verts[v2].halfEdge = h2;

                h->pair = h2; h2->pair = h;
                verts.insert(dup); verts.erase(v2); //v2 will remain at same position afterall.

                m.edges.push_back({std::min(dup,v2), std::max(dup, v2)});
            }
        }
    }
    for(int i = 0; i < boundaryHes.size(); i++)
    {
        HalfEdge *he = boundaryHes[i];
        HalfEdge *ogPair = he->pair;
        HalfEdge *ogNext = he->next;
        // Face *f = boundaryFaces[i];
        int v = he->head->id;
        int back = he->pair->head->id;
        int vDup = duplicatedVertices[v];
        int backDup = duplicatedVertices[back];

        // Create new half-edges and a quad face for the side
        HalfEdge *h1 = &m.halfEdges.push_back(HalfEdge());
        HalfEdge *h2 = heMap[make_pair(backDup, back)]; 
        HalfEdge *h3 = &m.halfEdges.push_back(HalfEdge());
        HalfEdge *h4 = heMap[make_pair(v, vDup)];

        Face *newFace = &m.faces.push_back(Face(4)); //create a new face for this face.
        newFace->halfEdge = h1;
        m.edges.push_back({std::min(vDup, backDup), std::max(vDup, backDup)});

        he->pair = h1; he->head = &m.verts[vDup];  
        
        // Connect half-edges to form a quad
        h1->head = &m.verts[backDup]; h1->next = h2; h1->face = newFace; h1->pair = he;
        
        h2->head = &m.verts[back]; h2->next = h3; h2->face = newFace; //pair already set up for this during creation.

        h3->head = &m.verts[v]; h3->next = h4; h3->face =newFace; h3->pair = ogPair;
        ogPair->pair = h3;

        h4->head = &m.verts[vDup]; h4->next = h1; h4->face = newFace; //pair already set up for this too.
     
        //after all this is done, we also should verify if the next and previous halfEdges are properly setup. 
        //since the next halfEdge might not be in the boundary. We will fix it and its pair in that case.
        if(boundaryHeSet.count(ogNext) == 0)
        {
            //ogNext is not in the boundary. So we need to create a new edge here and delete a previous edge.
            edgesToDelete.insert({std::min(ogNext->head->id, ogNext->pair->head->id), std::max(ogNext->head->id, ogNext->pair->head->id)});
            if(duplicatedVertices.count(ogNext->head->id))
            {
                ogNext->head = &m.verts[duplicatedVertices[ogNext->head->id]]; //now it will point to the duplicated vertex.
            }
            ogNext->pair->head = &m.verts[vDup];
            newEdges.push_back(ivec2(std::min(ogNext->head->id, vDup), std::max(ogNext->head->id, vDup)));
            // cout << "connected " << ogNext->head->id << " to " << vDup << endl;
        }
    }
 
    for(int i = 0; i < m.edges.size(); i++)
    {
        if(edgesToDelete.count({m.edges[i].x, m.edges[i].y}) == 0)
        {
            newEdges.push_back(m.edges[i]);
        }
    }
    m.edges = newEdges;
    // Step 4: Move duplicated boundary vertices
    for (int vid : verts) //verts contains the list of vertices we want to move.
    {
        // cout << "moving by " << offset << " in direction " << direction.x << " " << direction.y << " " << direction.z << endl;
        m.vertexPositions[vid] += offset * direction;
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
        ogPair->pair = h3; //also gotta update ogPair here.

        h4->head = vDup; h4->face = newFace; h4->next = h1; //h4's pair is not yet set, or in the case of startBackToBackDup, already set beforehand.

        //the new vertex needs to point to a halfEdge that points away from it. So for vDup, the best vertex to point to is:
        vDup->halfEdge = h1;

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

void revolveAndExtrude(Mesh &m, int faceid, float theta, float offset, int iterations, vec3 dir)
{
    //actually the face will remain the same.
    for(int iter = 0; iter < iterations; iter++)
    {
        extrude(m, offset, faceid, dir + faceNormal(m, &m.faces[faceid]));
        rotateFace(m, faceid, theta, vec3(0.0f, 1.0f, 0.0f));
    }
}


void catmullClarkSubdivision(Mesh& m, int iterations) {
    if (iterations>1){
        catmullClarkSubdivision(m, iterations-1);
    }
    std::vector<vec3> outVerts;
    std::vector<std::vector<int>> outFaces;
    outVerts.clear(); outFaces.clear();
    
    // Store original vertex positions
    const std::vector<vec3>& origVertPositions = m.vertexPositions;
    
    // Step 1: Create face points (one for each face)
    std::vector<vec3> facePoints;
    std::map<int, int> facePointIndices; // Maps original face index to new vertex index
    
    for (size_t faceIdx = 0; faceIdx < m.faces.size(); faceIdx++) {
        const Face& face = m.faces[faceIdx];
        vec3 facePoint(0, 0, 0);
        int vertCount = 0;
        
        // Get all vertices of this face
        HalfEdge* start = face.halfEdge;
        if (!start) continue;
        
        HalfEdge* current = start;
        do {
            if (current && current->head && current->head->id >= 0 && 
                current->head->id < origVertPositions.size()) {
                facePoint += origVertPositions[current->head->id];
                vertCount++;
            }
            current = current->next;
        } while (current && current != start);
        
        // Calculate average position
        if (vertCount > 0) {
            facePoint /= static_cast<float>(vertCount);
        }
        
        facePoints.push_back(facePoint);
        // The index of this face point will be: origVertCount + faceIdx
        facePointIndices[faceIdx] = origVertPositions.size() + facePoints.size() - 1;
    }
    
    // Step 2: Create edge points (one for each edge)
    std::map<std::pair<int, int>, int> edgePointIndices; // Maps edge (v1,v2) to new vertex index
    std::vector<vec3> edgePoints;
    
    for (size_t edgeIdx = 0; edgeIdx < m.halfEdges.size(); edgeIdx++) {
        const HalfEdge& edge = m.halfEdges[edgeIdx];
        if (!edge.head || !edge.pair || !edge.pair->head) continue;
        
        // Get vertex indices of this edge
        int v1 = edge.head->id;
        int v2 = edge.pair->head->id;
        
        // Create canonical edge key (smaller index first)
        std::pair<int, int> edgeKey = v1 < v2 ? std::make_pair(v1, v2) : std::make_pair(v2, v1);
        
        // Skip if we've already processed this edge
        if (edgePointIndices.find(edgeKey) != edgePointIndices.end()) continue;
        
        // Calculate edge point position (average of endpoints and adjoining face points)
        vec3 edgePoint(0, 0, 0);
        
        // Add the two endpoints
        if (v1 >= 0 && v1 < origVertPositions.size() && 
            v2 >= 0 && v2 < origVertPositions.size()) {
            edgePoint = origVertPositions[v1] + origVertPositions[v2];
        }
        
        // Find the face indices
        size_t face1Idx = SIZE_MAX, face2Idx = SIZE_MAX;
        
        // Find the faces connected to this edge
        if (edge.face) {
            for (size_t i = 0; i < m.faces.size(); i++) {
                if (&m.faces[i] == edge.face) {
                    face1Idx = i;
                    break;
                }
            }
        }
        
        if (edge.pair && edge.pair->face) {
            for (size_t i = 0; i < m.faces.size(); i++) {
                if (&m.faces[i] == edge.pair->face) {
                    face2Idx = i;
                    break;
                }
            }
        }
        
        // Add the face points if they exist
        if (face1Idx != SIZE_MAX && face1Idx < facePoints.size()) {
            edgePoint += facePoints[face1Idx];
        }
        
        if (face2Idx != SIZE_MAX && face2Idx < facePoints.size()) {
            edgePoint += facePoints[face2Idx];
        }
        
        // Calculate the average
        int divisor = 2; // Start with 2 for the endpoints
        if (face1Idx != SIZE_MAX) divisor++;
        if (face2Idx != SIZE_MAX) divisor++;
        
        edgePoint /= static_cast<float>(divisor);
        
        // Store the edge point
        edgePoints.push_back(edgePoint);
        // The index of this edge point will be: origVertCount + facePointCount + edgePointsSize - 1
        edgePointIndices[edgeKey] = origVertPositions.size() + facePoints.size() + edgePoints.size() - 1;
    }
    
    // Step 3: Update original vertex positions
    std::vector<vec3> updatedVertPositions(origVertPositions.size());
    
    for (size_t vertIdx = 0; vertIdx < m.verts.size(); vertIdx++) {
        const Vertex& vert = m.verts[vertIdx];
        if (vert.id < 0 || vert.id >= origVertPositions.size() || !vert.halfEdge) {
            if (vert.id >= 0 && vert.id < origVertPositions.size()) {
                updatedVertPositions[vert.id] = origVertPositions[vert.id];
            }
            continue;
        }
        
        // Calculate updated position using Catmull-Clark formula:
        // P' = (F + 2R + (n-3)P) / n
        // where F = average of face points, R = average of edge midpoints, P = original position, n = valence
        
        std::vector<vec3> adjacentFacePoints;
        std::vector<int> adjacentVertices;
        
        // Get all adjacent vertices and faces
        HalfEdge* start = vert.halfEdge;
        HalfEdge* current = start;
        
        do {
            if (current && current->face) {
                // Find the face index
                for (size_t i = 0; i < m.faces.size(); i++) {
                    if (&m.faces[i] == current->face) {
                        adjacentFacePoints.push_back(facePoints[i]);
                        break;
                    }
                }
            }
            
            if (current && current->pair && current->pair->head) {
                adjacentVertices.push_back(current->pair->head->id);
            }
            
            if (current && current->pair) {
                current = current->pair->next;
            } else {
                break;
            }
        } while (current && current != start);
        
        // If no adjacent faces or edges, keep original position
        if (adjacentFacePoints.empty() || adjacentVertices.empty()) {
            updatedVertPositions[vert.id] = origVertPositions[vert.id];
            continue;
        }
        
        // Calculate F (average of face points)
        vec3 F(0, 0, 0);
        for (const auto& fp : adjacentFacePoints) {
            F += fp;
        }
        F /= static_cast<float>(adjacentFacePoints.size());
        
        // Calculate R (average of edge midpoints)
        vec3 R(0, 0, 0);
        for (int adjVertId : adjacentVertices) {
            if (adjVertId >= 0 && adjVertId < origVertPositions.size()) {
                R += (origVertPositions[vert.id] + origVertPositions[adjVertId]) * 0.5f;
            }
        }
        R /= static_cast<float>(adjacentVertices.size());
        
        // Original position P
        const vec3& P = origVertPositions[vert.id];
        
        // Valence n
        int n = adjacentVertices.size();
        
        // Apply Catmull-Clark formula: P' = (F + 2R + (n-3)P) / n
        updatedVertPositions[vert.id] = (F + 2.0f * R + static_cast<float>(n - 3) * P) / static_cast<float>(n);
    }
    
    // Step 4: Build the output vertex list
    outVerts.clear();
    
    // Add updated original vertices
    outVerts.insert(outVerts.end(), updatedVertPositions.begin(), updatedVertPositions.end());
    
    // Add face points
    outVerts.insert(outVerts.end(), facePoints.begin(), facePoints.end());
    
    // Add edge points
    outVerts.insert(outVerts.end(), edgePoints.begin(), edgePoints.end());
    
    // Step 5: Create the new faces (quad for each original face edge)
    outFaces.clear();
    
    for (size_t faceIdx = 0; faceIdx < m.faces.size(); faceIdx++) {
        const Face& face = m.faces[faceIdx];
        if (!face.halfEdge) continue;
        
        // Find all half-edges that make up this face
        std::vector<HalfEdge*> faceEdges;
        HalfEdge* start = face.halfEdge;
        HalfEdge* current = start;
        
        do {
            faceEdges.push_back(current);
            current = current->next;
        } while (current && current != start);
        
        // Get the face point index
        int facePointIdx = facePointIndices[faceIdx];
        
        // For each edge, create a new quad face
        for (size_t i = 0; i < faceEdges.size(); i++) {
            HalfEdge* edge = faceEdges[i];
            if (!edge || !edge->head) continue;
            
            int vertIdx = edge->head->id;
            
            // Get next edge and next vertex
            HalfEdge* nextEdge = faceEdges[(i + 1) % faceEdges.size()];
            if (!nextEdge || !nextEdge->head) continue;
            
            int nextVertIdx = nextEdge->head->id;
            
            // Get edge point for current edge
            std::pair<int, int> edgeKey = vertIdx < nextVertIdx ? 
                std::make_pair(vertIdx, nextVertIdx) : 
                std::make_pair(nextVertIdx, vertIdx);
                
            auto edgeIt = edgePointIndices.find(edgeKey);
            if (edgeIt == edgePointIndices.end()) continue;
            int edgePointIdx = edgeIt->second;
            
            // Get edge point for previous edge
            int prevVertIdx = faceEdges[(i == 0 ? faceEdges.size() - 1 : i - 1)]->head->id;
            
            edgeKey = vertIdx < prevVertIdx ? 
                std::make_pair(vertIdx, prevVertIdx) : 
                std::make_pair(prevVertIdx, vertIdx);
                
            auto prevEdgeIt = edgePointIndices.find(edgeKey);
            if (prevEdgeIt == edgePointIndices.end()) continue;
            int prevEdgePointIdx = prevEdgeIt->second;
            
            // Create the quad face
            // Order vertices counterclockwise: vertex, current edge point, face point, previous edge point
            std::vector<int> newFace = {vertIdx, edgePointIdx, facePointIdx, prevEdgePointIdx};
            outFaces.push_back(newFace);
        }
    }
    
    // Finally, generate the new mesh from vertices and faces
    getMeshFromVerts(m, outVerts, outFaces);
}
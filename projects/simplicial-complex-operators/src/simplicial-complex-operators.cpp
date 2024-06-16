#include <iostream>
// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {

    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.
    
    size_t num_edges = mesh->nEdges();
    size_t num_vertices = mesh->nVertices();
    SparseMatrix<size_t> mat(num_edges, num_vertices);
    mat.reserve(2 * num_edges);
    for (Edge e : mesh->edges())
    {
        size_t i = e.getIndex();
        for (Vertex v : e.adjacentVertices())
        {
            size_t j = v.getIndex();
            mat.insert(i, j) = 1;
        }
    }
    return mat;
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {

    size_t num_faces = mesh->nFaces();
    size_t num_edges = mesh->nEdges();
    SparseMatrix<size_t> mat(num_faces, num_edges);
    mat.reserve(3 * num_faces);
    for (Face f : mesh->faces())
    {
        size_t i = f.getIndex();
        for (Edge e : f.adjacentEdges())
        {
            size_t j = e.getIndex();
            mat.insert(i, j) = 1;
        }
    }
    return mat;
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {

    Vector<size_t> vec = Vector<size_t>::Zero(mesh->nVertices());
    for (size_t v : subset.vertices)
    {
        vec[v] = 1;
    }
    return vec;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {

    Vector<size_t> vec = Vector<size_t>::Zero(mesh->nEdges());
    for (size_t e : subset.edges)
    {
        vec[e] = 1;
    }
    return vec;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {

    Vector<size_t> vec = Vector<size_t>::Zero(mesh->nFaces());
    for (size_t f : subset.faces)
    {
        vec[f] = 1;
    }
    return vec;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {

    MeshSubset str = subset.deepCopy();
    SparseMatrix<size_t> A10 = A1 * A0;
    
    for (size_t vertex : subset.vertices)
    {
        for (SparseMatrix<size_t>::InnerIterator it(A0, vertex); it; ++it)
        {
            str.addEdge(it.index());
        }
        
        for (SparseMatrix<size_t>::InnerIterator it(A10, vertex); it; ++it)
        {
            str.addFace(it.index());
        }
    }

    for (size_t edge : subset.edges)
    {
        for (SparseMatrix<size_t>::InnerIterator it(A1, edge); it; ++it)
        {
            str.addFace(it.index());
        }
    }
    return str;
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {

    MeshSubset cl = subset.deepCopy();
    SparseMatrix<size_t> A0t = A0.transpose();
    SparseMatrix<size_t> A1t = A1.transpose();
    SparseMatrix<size_t> A10t = (A1*A0).transpose();
    for (size_t edge : subset.edges)
    {
        for (SparseMatrix<size_t>::InnerIterator it(A0t, edge); it; ++it)
        {
            cl.addVertex(it.index());
        }
    }
    for (size_t face : subset.faces)
    {
        for (SparseMatrix<size_t>::InnerIterator it(A1t, face); it; ++it)
        {
            cl.addEdge(it.index());
        }
        for (SparseMatrix<size_t>::InnerIterator it(A10t, face); it; ++it)
        {
            cl.addVertex(it.index());
        }
    }
    return cl;
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {

    MeshSubset str = star(subset);
    MeshSubset cl = closure(str);
    MeshSubset lk = cl;
    for (size_t v : str.vertices)
    {
        lk.vertices.erase(v);
    }
    for (size_t e : str.edges)
    {
        lk.edges.erase(e);
    }
    for (size_t f : str.faces)
    {
        lk.faces.erase(f);
    }
    return lk;
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    MeshSubset cl = closure(subset);
    return cl.vertices == subset.vertices &&
            cl.edges == subset.edges;
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {

    if (!isComplex(subset))
    {
        return -1;
    }

    if (subset.edges.empty() && subset.faces.empty())
    {
        return 0;
    }
    if (!subset.edges.empty())
    {
        for (size_t v : subset.vertices)
        {
            int count = 0;
            for (SparseMatrix<size_t>::InnerIterator it(A0, v); it; ++it)
            {
                size_t incident_edge = it.index();
                if (subset.edges.find(incident_edge) != subset.edges.end())
                {
                    count += 1;
                }
            }
            if (count == 0)
            {
                return -1;
            }
        }
    }
    if (subset.faces.empty())
    {
        return 1;
    }
    for (size_t e : subset.edges)
    {
        int count = 0;
        for (SparseMatrix<size_t>::InnerIterator it(A1, e); it; ++it)
        {
            size_t incident_face = it.index();
            if (subset.faces.find(incident_face) != subset.faces.end())
            {
                count += 1;
            }
        }
        if (count == 0)
        {
            return -1;
        }
    }
    return 2;
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {

    MeshSubset bdry;
    if (isPureComplex(subset) == -1)
    {
        return bdry;
    }
    for (size_t edge : subset.edges)
    {
        int count = 0;
        for (SparseMatrix<size_t>::InnerIterator it(A1, edge); it; ++it)
        {
            if (subset.faces.find(it.index()) != subset.faces.end())
            {
                count += 1;
            }
        }
        if (count == 1)
        {
            bdry.edges.insert(edge);
        }
    }
    return closure(bdry);
}
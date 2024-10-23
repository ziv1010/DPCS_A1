// wireframe.cpp

#include "wireframe.h"
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>

// Constructor
Wireframe::Wireframe(const Projection2D& front, const Projection2D& top, const Projection2D& side)
    : frontView(front), topView(top), sideView(side) {
    generate3DVertices();
    generate3DEdges();
    removeRedundantEdges();
    removePathologicalElements();
}

// Helper function to find if a 3D vertex already exists
int Wireframe::findVertex3DIndex(float x, float y, float z) const {
    for (size_t i = 0; i < vertices3D.size(); ++i) {
        if (std::abs(vertices3D[i].getX() - x) < 1e-6 &&
            std::abs(vertices3D[i].getY() - y) < 1e-6 &&
            std::abs(vertices3D[i].getZ() - z) < 1e-6) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

// Generate 3D vertices by matching 2D vertices from projections
void Wireframe::generate3DVertices() {
    // Mapping from projection vertices to 3D vertices
    std::unordered_map<int, int> frontTo3D, topTo3D, sideTo3D;

    // For simplicity, assume that the projections are aligned and matching by indices
    // In practice, more sophisticated matching is required
    size_t numVertices = std::min({ frontView.vertices.size(), topView.vertices.size(), sideView.vertices.size() });

    for (size_t i = 0; i < numVertices; ++i) {
        const Vertex2D& vf = frontView.vertices[i];
        const Vertex2D& vt = topView.vertices[i];
        const Vertex2D& vs = sideView.vertices[i];

        float x = vf.getX();
        float y = vt.getY(); // From top view (Y axis)
        float z = vs.getY(); // From side view (Z axis)

        // Check if vertex already exists
        int existingIndex = findVertex3DIndex(x, y, z);
        if (existingIndex == -1) {
            vertices3D.emplace_back(x, y, z);
            existingIndex = static_cast<int>(vertices3D.size() - 1);
        }

        frontTo3D[i] = existingIndex;
        topTo3D[i] = existingIndex;
        sideTo3D[i] = existingIndex;
    }
}

// Generate 3D edges by combining 2D edges from projections
void Wireframe::generate3DEdges() {
    // A simple approach: for each edge in front view, check if corresponding edges exist in top and side views
    for (const auto& fe : frontView.edges) {
        int fv1_idx = fe.v1_index;
        int fv2_idx = fe.v2_index;

        // Get corresponding 3D vertex indices
        if (fv1_idx >= vertices3D.size() || fv2_idx >= vertices3D.size()) continue;

        int v1_idx = fv1_idx;
        int v2_idx = fv2_idx;

        // Check if this edge exists in top and side views
        bool edgeInTop = false;
        bool edgeInSide = false;

        for (const auto& te : topView.edges) {
            if ((te.v1_index == fv1_idx && te.v2_index == fv2_idx) ||
                (te.v1_index == fv2_idx && te.v2_index == fv1_idx)) {
                edgeInTop = true;
                break;
            }
        }

        for (const auto& se : sideView.edges) {
            if ((se.v1_index == fv1_idx && se.v2_index == fv2_idx) ||
                (se.v1_index == fv2_idx && se.v2_index == fv1_idx)) {
                edgeInSide = true;
                break;
            }
        }

        if (edgeInTop && edgeInSide) {
            // Add edge if it exists in all projections
            edges3D.emplace_back(v1_idx, v2_idx);
            // Update connectivity
            vertices3D[v1_idx].addConnectedEdge(static_cast<int>(edges3D.size() - 1));
            vertices3D[v2_idx].addConnectedEdge(static_cast<int>(edges3D.size() - 1));
        }
    }

    // Repeat for edges in top and side views to ensure all edges are captured
    // This can be enhanced based on specific requirements
}

// Remove redundant edges using the RER procedure
void Wireframe::removeRedundantEdges() {
    // Mark all edges as unexamined
    std::vector<bool> examined(edges3D.size(), false);

    for (size_t i = 0; i < edges3D.size(); ++i) {
        if (examined[i]) continue;

        const Edge3D& e1 = edges3D[i];
        examined[i] = true;

        std::vector<size_t> overlappingEdges;

        for (size_t j = i + 1; j < edges3D.size(); ++j) {
            if (examined[j]) continue;

            const Edge3D& e2 = edges3D[j];

            // Check if edges overlap (same vertices, regardless of direction)
            if ((e1.v1_index == e2.v1_index && e1.v2_index == e2.v2_index) ||
                (e1.v1_index == e2.v2_index && e1.v2_index == e2.v1_index)) {
                overlappingEdges.push_back(j);
            }
        }

        // Remove overlapping edges
        // Iterate in reverse to avoid invalidating indices
        for (auto it = overlappingEdges.rbegin(); it != overlappingEdges.rend(); ++it) {
            int idx = static_cast<int>(*it);
            // Remove edge from connected vertices
            vertices3D[edges3D[idx].v1_index].removeConnectedEdge(idx);
            vertices3D[edges3D[idx].v2_index].removeConnectedEdge(idx);
            // Remove edge from edges3D
            edges3D.erase(edges3D.begin() + idx);
            // Update connectivity indices in all vertices
            for (auto& vertex : vertices3D) {
                for (auto& edgeIdx : vertex.connectedEdges) {
                    if (edgeIdx > idx) {
                        edgeIdx--;
                    }
                }
            }
            // Adjust examined vector
            examined.erase(examined.begin() + idx);
        }
    }
}

// Helper function to check if two edges are collinear
bool Wireframe::areCollinear(int edge1, int edge2) const {
    if (edge1 >= edges3D.size() || edge2 >= edges3D.size()) return false;

    const Vertex3D& v1_start = vertices3D[edges3D[edge1].v1_index];
    const Vertex3D& v1_end = vertices3D[edges3D[edge1].v2_index];
    const Vertex3D& v2_start = vertices3D[edges3D[edge2].v1_index];
    const Vertex3D& v2_end = vertices3D[edges3D[edge2].v2_index];

    // Compute direction vectors
    float dx1 = v1_end.x - v1_start.x;
    float dy1 = v1_end.y - v1_start.y;
    float dz1 = v1_end.z - v1_start.z;

    float dx2 = v2_end.x - v2_start.x;
    float dy2 = v2_end.y - v2_start.y;
    float dz2 = v2_end.z - v2_start.z;

    // Compute cross product
    float cross_x = dy1 * dz2 - dz1 * dy2;
    float cross_y = dz1 * dx2 - dx1 * dz2;
    float cross_z = dx1 * dy2 - dy1 * dx2;

    // If cross product is (close to) zero vector, they are collinear
    return (std::abs(cross_x) < 1e-6 &&
            std::abs(cross_y) < 1e-6 &&
            std::abs(cross_z) < 1e-6);
}

// Helper function to check if a set of edges are coplanar
bool Wireframe::areCoplanar(const std::vector<int>& edgeIndices) const {
    if (edgeIndices.size() < 3) return true; // Any two edges are always coplanar

    // Get three distinct points to define a plane
    const Vertex3D& p1 = vertices3D[edges3D[edgeIndices[0]].v1_index];
    const Vertex3D& p2 = vertices3D[edges3D[edgeIndices[0]].v2_index];
    const Vertex3D& p3 = vertices3D[edges3D[edgeIndices[1]].v2_index];

    // Define vectors
    float v1x = p2.x - p1.x;
    float v1y = p2.y - p1.y;
    float v1z = p2.z - p1.z;

    float v2x = p3.x - p1.x;
    float v2y = p3.y - p1.y;
    float v2z = p3.z - p1.z; // Corrected from p3.z - p3.z

    // Compute normal vector of the plane
    float nx = v1y * v2z - v1z * v2y;
    float ny = v1z * v2x - v1x * v2z;
    float nz = v1x * v2y - v1y * v2x;

    // Normalize the normal vector
    float length = std::sqrt(nx * nx + ny * ny + nz * nz);
    if (length < 1e-6) return false; // Degenerate case, cannot define plane

    nx /= length;
    ny /= length;
    nz /= length;

    // Compute d using p1
    float d = nx * p1.x + ny * p1.y + nz * p1.z;

    // Check if all other edges lie on this plane
    for (size_t i = 0; i < edgeIndices.size(); ++i) {
        const Edge3D& edge = edges3D[edgeIndices[i]];
        const Vertex3D& va = vertices3D[edge.v1_index];
        const Vertex3D& vb = vertices3D[edge.v2_index];

        float dist_a = nx * va.x + ny * va.y + nz * va.z - d;
        float dist_b = nx * vb.x + ny * vb.y + nz * vb.z - d;

        if (std::abs(dist_a) > 1e-6 || std::abs(dist_b) > 1e-6) {
            return false;
        }
    }

    return true;
}

// Helper function to merge two collinear edges into one
void Wireframe::mergeEdges(int edge1, int edge2) {
    if (edge1 >= edges3D.size() || edge2 >= edges3D.size()) return;

    const Edge3D& e1 = edges3D[edge1];
    const Edge3D& e2 = edges3D[edge2];

    // Determine the common vertex
    int common = -1;
    if (e1.v1_index == e2.v1_index || e1.v1_index == e2.v2_index) {
        common = e1.v1_index;
    }
    else if (e1.v2_index == e2.v1_index || e1.v2_index == e2.v2_index) {
        common = e1.v2_index;
    }

    if (common == -1) return; // No common vertex, cannot merge

    // Determine the non-common vertices
    int non_common_e1 = (e1.v1_index == common) ? e1.v2_index : e1.v1_index;
    int non_common_e2 = (e2.v1_index == common) ? e2.v2_index : e2.v1_index;

    // Create a new edge connecting the two non-common vertices
    Edge3D newEdge(non_common_e1, non_common_e2); // Declaration
    edges3D.emplace_back(newEdge);
    int newEdgeIdx = static_cast<int>(edges3D.size() - 1);

    // Update connectivity
    vertices3D[non_common_e1].addConnectedEdge(newEdgeIdx);
    vertices3D[non_common_e2].addConnectedEdge(newEdgeIdx);

    // Remove the old edges from connected vertices
    vertices3D[e1.v1_index].removeConnectedEdge(edge1);
    vertices3D[e1.v2_index].removeConnectedEdge(edge1);
    vertices3D[e2.v1_index].removeConnectedEdge(edge2);
    vertices3D[e2.v2_index].removeConnectedEdge(edge2);

    // Remove the old edges from edges3D
    // Remove higher index first to avoid invalidating the lower index
    if (edge1 > edge2) {
        edges3D.erase(edges3D.begin() + edge1);
        edges3D.erase(edges3D.begin() + edge2);
    }
    else {
        edges3D.erase(edges3D.begin() + edge2);
        edges3D.erase(edges3D.begin() + edge1);
    }

    // Update connectivity indices in all vertices
    for (auto& vertex : vertices3D) {
        for (auto& eIdx : vertex.connectedEdges) {
            if (eIdx > edge1) eIdx--;
            if (eIdx > edge2) eIdx--;
            // New edge is already added at the end, no need to adjust
        }
    }
}

// Remove pathological edges and vertices using the PEVR procedure
void Wireframe::removePathologicalElements() {
    bool changesMade = true;

    while (changesMade) {
        changesMade = false;

        for (size_t i = 0; i < vertices3D.size(); ++i) {
            Vertex3D& vertex = vertices3D[i];
            int degree = vertex.degree();

            if (degree == 0) {
                // W3.1: Delete isolated vertex
                // Remove vertex and continue
                vertices3D.erase(vertices3D.begin() + i);
                // No edges connected, so no need to update edges
                // Update edge indices in all vertices
                for (auto& v : vertices3D) {
                    for (auto& edgeIdx : v.connectedEdges) {
                        if (edgeIdx > static_cast<int>(i)) {
                            edgeIdx--;
                        }
                    }
                }
                changesMade = true;
                i--; // Adjust index after erase
                continue;
            }
            else if (degree == 1) {
                // W3.2: Remove dangling edge and vertex
                int edgeIdx = vertex.connectedEdges[0];
                Edge3D& edge = edges3D[edgeIdx];
                // Find the other vertex
                int otherV = (edge.v1_index == static_cast<int>(i)) ? edge.v2_index : edge.v1_index;
                // Remove the edge from the other vertex's connectedEdges
                vertices3D[otherV].removeConnectedEdge(edgeIdx);
                // Remove the edge from edges3D
                edges3D.erase(edges3D.begin() + edgeIdx);
                // Remove the edge from current vertex's connectedEdges
                vertex.removeConnectedEdge(edgeIdx);
                // Remove vertex
                vertices3D.erase(vertices3D.begin() + i);
                // Update edge indices in all vertices
                for (auto& v : vertices3D) {
                    for (auto& eIdx : v.connectedEdges) {
                        if (eIdx > edgeIdx) {
                            eIdx--;
                        }
                    }
                }
                changesMade = true;
                i--; // Adjust index after erase
                continue;
            }
            else if (degree == 2) {
                // W3.3: Handle degree 2 vertices
                int edge1 = vertex.connectedEdges[0];
                int edge2 = vertex.connectedEdges[1];

                if (areCollinear(edge1, edge2)) {
                    // Merge the two collinear edges
                    mergeEdges(edge1, edge2);
                    changesMade = true;
                    i--; // Recheck the current index after modification
                    break;
                }
                else {
                    // Remove both edges and the vertex
                    // Remove edges from connected vertices
                    Edge3D& e1 = edges3D[edge1];
                    Edge3D& e2 = edges3D[edge2];

                    int otherV1 = (e1.v1_index == static_cast<int>(i)) ? e1.v2_index : e1.v1_index;
                    int otherV2 = (e2.v1_index == static_cast<int>(i)) ? e2.v2_index : e2.v1_index;

                    vertices3D[otherV1].removeConnectedEdge(edge1);
                    vertices3D[otherV2].removeConnectedEdge(edge2);

                    // Remove edges from edges3D
                    // Remove higher index first
                    if (edge1 > edge2) {
                        edges3D.erase(edges3D.begin() + edge1);
                        edges3D.erase(edges3D.begin() + edge2);
                    }
                    else {
                        edges3D.erase(edges3D.begin() + edge2);
                        edges3D.erase(edges3D.begin() + edge1);
                    }

                    // Remove vertex
                    vertices3D.erase(vertices3D.begin() + i);
                    // Update edge indices in all vertices
                    for (auto& v : vertices3D) {
                        for (auto& eIdx : v.connectedEdges) {
                            if (eIdx > edge1) eIdx--;
                            if (eIdx > edge2) eIdx--;
                        }
                    }
                    changesMade = true;
                    i--; // Adjust index after erase
                    break;
                }
            }
            else if (degree == 3) {
                // W3.4 and W3.5: Handle degree 3 vertices
                std::vector<int> connectedEdges = vertex.connectedEdges;
                bool twoCollinear = false;
                int collinearEdge1 = -1, collinearEdge2 = -1;
                for (size_t m = 0; m < connectedEdges.size(); ++m) {
                    for (size_t n = m + 1; n < connectedEdges.size(); ++n) {
                        if (areCollinear(connectedEdges[m], connectedEdges[n])) {
                            twoCollinear = true;
                            collinearEdge1 = connectedEdges[m];
                            collinearEdge2 = connectedEdges[n];
                            break;
                        }
                    }
                    if (twoCollinear) break;
                }

                if (twoCollinear) {
                    // W3.4: Two edges are collinear
                    // Merge the two collinear edges and remove the third edge
                    int thirdEdge = -1;
                    for (const auto& eIdx : connectedEdges) {
                        if (eIdx != collinearEdge1 && eIdx != collinearEdge2) {
                            thirdEdge = eIdx;
                            break;
                        }
                    }
                    if (thirdEdge != -1) {
                        // Remove the third edge
                        Edge3D& te = edges3D[thirdEdge];
                        int otherV = (te.v1_index == static_cast<int>(i)) ? te.v2_index : te.v1_index;
                        vertices3D[otherV].removeConnectedEdge(thirdEdge);
                        // Remove the edge from edges3D
                        edges3D.erase(edges3D.begin() + thirdEdge);
                        // Remove the edge from current vertex's connectedEdges
                        vertex.removeConnectedEdge(thirdEdge);
                        // Merge collinear edges
                        mergeEdges(collinearEdge1, collinearEdge2);
                        changesMade = true;
                        i--; // Recheck after modification
                        break;
                    }
                }
                else {
                    // Check coplanarity
                    bool coplanar = areCoplanar(connectedEdges);
                    if (coplanar) {
                        // W3.5: Three edges are coplanar but no two are collinear
                        // Remove the vertex and all three edges
                        for (const auto& eIdx : connectedEdges) {
                            Edge3D& e = edges3D[eIdx];
                            int otherV = (e.v1_index == static_cast<int>(i)) ? e.v2_index : e.v1_index;
                            vertices3D[otherV].removeConnectedEdge(eIdx);
                        }
                        // Remove edges in reverse order to maintain correct indices
                        std::sort(connectedEdges.begin(), connectedEdges.end(), std::greater<int>());
                        for (const auto& eIdx : connectedEdges) {
                            edges3D.erase(edges3D.begin() + eIdx);
                        }
                        // Remove vertex
                        vertices3D.erase(vertices3D.begin() + i);
                        // Update edge indices in all vertices
                        for (auto& v : vertices3D) {
                            for (auto& eIdx : v.connectedEdges) {
                                for (const auto& removedIdx : connectedEdges) {
                                    if (eIdx > removedIdx) {
                                        eIdx--;
                                    }
                                }
                            }
                        }
                        changesMade = true;
                        i--; // Adjust index after erase
                        break;
                    }
                }
            }
            else if (degree >= 4) {
                // W3.6: Handle degree >=4 vertices
                // For simplicity, handle only the first four edges
                std::vector<int> firstFourEdges;
                for (int eIdx : vertex.connectedEdges) {
                    firstFourEdges.push_back(eIdx);
                    if (firstFourEdges.size() == 4) break;
                }

                if (firstFourEdges.size() < 4) continue;

                bool coplanar = areCoplanar(firstFourEdges);
                if (coplanar) {
                    // Count number of collinear pairs
                    int collinearPairs = 0;
                    for (size_t m = 0; m < firstFourEdges.size(); ++m) {
                        for (size_t n = m + 1; n < firstFourEdges.size(); ++n) {
                            if (areCollinear(firstFourEdges[m], firstFourEdges[n])) {
                                collinearPairs++;
                            }
                        }
                    }

                    if (collinearPairs == 1) {
                        // Only two edges are collinear
                        // Merge the two collinear edges and remove two non-collinear edges
                        int collinearEdge1 = -1, collinearEdge2 = -1;
                        for (size_t m = 0; m < firstFourEdges.size(); ++m) {
                            for (size_t n = m + 1; n < firstFourEdges.size(); ++n) {
                                if (areCollinear(firstFourEdges[m], firstFourEdges[n])) {
                                    collinearEdge1 = firstFourEdges[m];
                                    collinearEdge2 = firstFourEdges[n];
                                    break;
                                }
                            }
                            if (collinearEdge1 != -1) break;
                        }

                        if (collinearEdge1 != -1 && collinearEdge2 != -1) {
                            // Identify two non-collinear edges
                            std::vector<int> nonCollinearEdges;
                            for (const auto& eIdx : firstFourEdges) {
                                if (eIdx != collinearEdge1 && eIdx != collinearEdge2) {
                                    nonCollinearEdges.push_back(eIdx);
                                }
                            }

                            if (nonCollinearEdges.size() >= 2) {
                                int nonCollinearEdge1 = nonCollinearEdges[0];
                                int nonCollinearEdge2 = nonCollinearEdges[1];

                                // Merge collinear edges
                                mergeEdges(collinearEdge1, collinearEdge2);

                                // Remove non-collinear edges
                                Edge3D& e1 = edges3D[nonCollinearEdge1];
                                Edge3D& e2 = edges3D[nonCollinearEdge2];
                                int otherV1 = (e1.v1_index == static_cast<int>(i)) ? e1.v2_index : e1.v1_index;
                                int otherV2 = (e2.v1_index == static_cast<int>(i)) ? e2.v2_index : e2.v1_index;

                                vertices3D[otherV1].removeConnectedEdge(nonCollinearEdge1);
                                vertices3D[otherV2].removeConnectedEdge(nonCollinearEdge2);

                                // Remove edges from edges3D
                                // Remove higher index first
                                if (nonCollinearEdge1 > nonCollinearEdge2) {
                                    edges3D.erase(edges3D.begin() + nonCollinearEdge1);
                                    edges3D.erase(edges3D.begin() + nonCollinearEdge2);
                                }
                                else {
                                    edges3D.erase(edges3D.begin() + nonCollinearEdge2);
                                    edges3D.erase(edges3D.begin() + nonCollinearEdge1);
                                }

                                // Remove vertex
                                vertices3D.erase(vertices3D.begin() + i);

                                // Update edge indices in all vertices
                                for (auto& v : vertices3D) {
                                    for (auto& eIdx : v.connectedEdges) {
                                        if (eIdx > nonCollinearEdge1) eIdx--;
                                        if (eIdx > nonCollinearEdge2) eIdx--;
                                    }
                                }

                                changesMade = true;
                                i--; // Adjust index after erase
                                break;
                            }
                        }
                    }
                    else if (collinearPairs >= 2) {
                        // Two pairs of collinear edges or all edges are non-collinear
                        // Remove vertex and all connected edges
                        for (const auto& eIdx : vertex.connectedEdges) {
                            Edge3D& e = edges3D[eIdx];
                            int otherV = (e.v1_index == static_cast<int>(i)) ? e.v2_index : e.v1_index;
                            vertices3D[otherV].removeConnectedEdge(eIdx);
                        }
                        // Remove edges in reverse order to maintain correct indices
                        std::vector<int> connectedEdgesCopy = vertex.connectedEdges;
                        std::sort(connectedEdgesCopy.begin(), connectedEdgesCopy.end(), std::greater<int>());
                        for (const auto& eIdx : connectedEdgesCopy) {
                            edges3D.erase(edges3D.begin() + eIdx);
                        }
                        // Remove vertex
                        vertices3D.erase(vertices3D.begin() + i);
                        // Update edge indices in all vertices
                        for (auto& v : vertices3D) {
                            for (auto& eIdx : v.connectedEdges) {
                                for (const auto& removedIdx : connectedEdgesCopy) {
                                    if (eIdx > removedIdx) {
                                        eIdx--;
                                    }
                                }
                            }
                        }
                        changesMade = true;
                        i--; // Adjust index after erase
                        break;
                    }
                }
            }
        }
    }
}
// Save wireframe to file
// Definition in wireframe.cpp
void Wireframe::saveToFile(const std::string& filename) const {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: Unable to write to file " << filename << std::endl;
        return;
    }

    // Output vertices
    int numVertices = static_cast<int>(vertices3D.size());
    outfile << numVertices << std::endl;
    for (const auto& vertex : vertices3D) {
        outfile << vertex.getX() << " " << vertex.getY() << " " << vertex.getZ() << std::endl;
    }

    // Output edges
    int numEdges = static_cast<int>(edges3D.size());
    outfile << numEdges << std::endl;
    for (const auto& edge : edges3D) {
        outfile << edge.v1_index << " " << edge.v2_index << std::endl;
    }

    outfile.close();
}

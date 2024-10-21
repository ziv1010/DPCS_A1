// 3Dto2DModel.cpp
#include <set>    // For std::set
#include <map>    // For std::map
#include <cmath>  // For mathematical functions like abs, sqrt
#include "3Dto2DModel.h"

// Default Constructor
Model_2D::Model_2D() {
    // Initialize direction to an undefined state
    direction = -1;

    // Initialize empty lists for vertices, edges, and surfaces
    v = std::vector<Vertex>();
    e = std::vector<Edge>();
    s = std::vector<Surface>();
}

// Add Multiple Vertices to the Model
void Model_2D::add_vertices(const std::vector<Vertex>& vt) {
    // Adds a list of Vertex objects to the Model_2D's vertex list.
    v.insert(v.end(), vt.begin(), vt.end());
}

// Add Multiple Edges to the Model
void Model_2D::add_edges(const std::vector<Edge>& et) {
    // Adds a list of Edge objects to the Model_2D's edge list.
    e.insert(e.end(), et.begin(), et.end());
}

// Add a Single Vertex to the Model
void Model_2D::add_vertex(const Vertex& vt) {
    // Adds a single Vertex object to the Model_2D's vertex list.
    v.push_back(vt);
}

// Retrieve All Vertices in the Model
const std::vector<Vertex>& Model_2D::ret_vertices() const {
    // Returns the list of Vertex objects in the Model_2D.
    return v;
}

// Retrieve All Edges in the Model
const std::vector<Edge>& Model_2D::ret_edges() const {
    // Returns the list of Edge objects in the Model_2D.
    return e;
}

// New method: Read model data from a file
void Model_2D::read_from_file(std::ifstream& inFile) {
    // Read number of vertices
    int nv;
    inFile >> nv;
    std::cout << "Number of vertices: " << nv << std::endl;

    // Read vertices
    for (int i = 0; i < nv; ++i) {
    int vNo;
    float coord1, coord2;
    inFile >> vNo >> coord1 >> coord2;
    Vertex vertex;

    // Assign coordinates based on direction
    if (direction == TOP_VIEW) {
        vertex = Vertex(coord1, coord2, 0.0f, vNo);
    } else if (direction == FRONT_VIEW) {
        vertex = Vertex(coord1, 0.0f, coord2, vNo);
    } else if (direction == SIDE_VIEW) {
        vertex = Vertex(0.0f, coord1, coord2, vNo);
    } else {
        std::cerr << "Invalid direction: " << direction << std::endl;
        return;
    }

    // Assign the vertex number
    vertex.vNo = vNo;

    v.push_back(vertex);
}

    // Read number of edges
    int ne;
    inFile >> ne;
    std::cout << "Number of edges: " << ne << std::endl;

    // Read edges
    for (int i = 0; i < ne; ++i) {
    int eno, vNoA, vNoB, hiddenFlag;
    inFile >> eno >> vNoA >> vNoB >> hiddenFlag;

    // Find vertices in v
    Vertex* a = nullptr;
    Vertex* b = nullptr;
    for (auto& vert : v) {
        if (vert.vNo == vNoA) {
            a = &vert;
        }
        if (vert.vNo == vNoB) {
            b = &vert;
        }
        if (a && b) {
            break;
        }
    }

    if (a && b) {
        Edge edge(*a, *b, eno);
        edge.hidden = (hiddenFlag == 1);
        e.push_back(edge);
    } else {
        std::cerr << "Error: Edge " << eno << " references non-existent vertices." << std::endl;
    }
}
// After reading vertices
for (const auto& vert : v) {
    std::cout << "Vertex vNo: " << vert.vNo << ", Coordinates: (" << vert.x << ", " << vert.y << ", " << vert.z << ")" << std::endl;
}
}


// Check if Two Edges Intersect and Return the Intersection Vertex
Vertex* Model_2D::is_intersect(int p, int q) {
    // Ensure indices are within bounds
    if (p >= e.size() || q >= e.size()) {
        return nullptr;
    }

    // Coordinates for the first edge
    float px1, py1, px2, py2;
    // Coordinates for the second edge
    float qx1, qy1, qx2, qy2;

    // Determine the projection direction and assign coordinates accordingly
    switch (direction) {
        case TOP_VIEW:    // Top view (x, y)
            px1 = e[p].a.x; py1 = e[p].a.y;
            px2 = e[p].b.x; py2 = e[p].b.y;

            qx1 = e[q].a.x; qy1 = e[q].a.y;
            qx2 = e[q].b.x; qy2 = e[q].b.y;

            // Set Z-coordinates to 0 as per top view projection
            e[p].a.z = e[p].b.z = e[q].a.z = e[q].b.z = 0.0f;
            break;

        case FRONT_VIEW:    // Front view (x, z)
            px1 = e[p].a.x; py1 = e[p].a.z;
            px2 = e[p].b.x; py2 = e[p].b.z;

            qx1 = e[q].a.x; qy1 = e[q].a.z;
            qx2 = e[q].b.x; qy2 = e[q].b.z;

            // Set Y-coordinates to 0 as per front view projection
            e[p].a.y = e[p].b.y = e[q].a.y = e[q].b.y = 0.0f;
            break;

        case SIDE_VIEW:    // Side view (y, z)
            px1 = e[p].a.y; py1 = e[p].a.z;
            px2 = e[p].b.y; py2 = e[p].b.z;

            qx1 = e[q].a.y; qy1 = e[q].a.z;
            qx2 = e[q].b.y; qy2 = e[q].b.z;

            // Set X-coordinates to 0 as per side view projection
            e[p].a.x = e[p].b.x = e[q].a.x = e[q].b.x = 0.0f;
            break;

        default:
            // Undefined direction
            return nullptr;
    }

    // Calculate denominators
    float den_s = (qx2 - qx1) * (py2 - py1) - (qy2 - qy1) * (px2 - px1);
    float den_t = den_s;

    // Avoid division by zero
    if (std::abs(den_s) < EPS || std::abs(den_t) < EPS) {
        return nullptr;
    }

    // Calculate numerators
    float num_s = (qx2 - qx1) * (py1 - qy1) - (qy2 - qy1) * (px1 - qx1);
    float num_t = (px2 - px1) * (py1 - qy1) - (py2 - py1) * (px1 - qx1);

    // Calculate parameters s and t
    float s_j = num_s / den_s;
    float t_j = num_t / den_t;

    // Check for valid intersection conditions
    if (s_j > 0.0f && s_j < 1.0f && t_j > 0.0f && t_j < 1.0f) {
        // Calculate the exact intersection point coordinates
        float x_intersection = px1 + (px2 - px1) * s_j;
        float y_intersection = py1 + (py2 - py1) * s_j;

        // Create a new Vertex object for the intersection
        Vertex v_intersect;

        // Assign coordinates based on the projection direction
        switch (direction) {
            case TOP_VIEW:
                v_intersect.x = x_intersection;
                v_intersect.y = y_intersection;
                v_intersect.z = 0.0f;
                break;

            case FRONT_VIEW:
                v_intersect.x = x_intersection;
                v_intersect.y = 0.0f;
                v_intersect.z = y_intersection;
                break;

            case SIDE_VIEW:
                v_intersect.x = 0.0f;
                v_intersect.y = x_intersection;
                v_intersect.z = y_intersection;
                break;

            default:
                // Undefined direction
                return nullptr;
        }

        // Assign a new vertex number based on the current vertex list size
        v_intersect.vNo = v.size() + 1;
        v_intersect.isTrue = true;

        // Add the intersection vertex to the model
        add_vertex(v_intersect);

        // Split the original edges at the intersection point

        // Create new edges replacing the first original edge
        Edge pe1(e[p].a, v_intersect, e.size() + 1);
        Edge pe2(v_intersect, e[p].b, e.size() + 2);

        // Create new edges replacing the second original edge
        Edge qe1(e[q].a, v_intersect, e.size() + 3);
        Edge qe2(v_intersect, e[q].b, e.size() + 4);

        // Add the new edges to the edge list
        e.push_back(pe1);
        e.push_back(pe2);
        e.push_back(qe1);
        e.push_back(qe2);

        // Mark the original edges as inactive
        e[p].isTrue = false;
        e[q].isTrue = false;

        // Return the intersection vertex
        // Return pointer to the newly added vertex
        return &v.back();
    } else {
        // No valid intersection detected
        return nullptr;
    }
}

// Check if a Point Lies on a Line Segment
bool Model_2D::onSegment(const glm::vec2& p, const glm::vec2& q, const glm::vec2& r) const {
    // Check if Q's coordinates are within the bounding rectangle of P and R
    return (q.x <= std::max(p.x, r.x) + EPS && q.x >= std::min(p.x, r.x) - EPS &&
            q.y <= std::max(p.y, r.y) + EPS && q.y >= std::min(p.y, r.y) - EPS);
}

// Determine the Orientation of Three Points (p, q, r)
int Model_2D::orientation(const glm::vec2& p, const glm::vec2& q, const glm::vec2& r) const {
    // Calculate the orientation value
    float val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);

    if (std::abs(val) < EPS) {
        return 0; // Colinear
    }
    return (val > 0) ? 1 : 2; // Clockwise or Counterclockwise
}

// Check if Two Line Segments (p1, q1) and (p2, q2) Intersect
bool Model_2D::doIntersect(const glm::vec2& p1, const glm::vec2& q1,
                           const glm::vec2& p2, const glm::vec2& q2) const {
    // Find the four orientations needed for the general and special cases
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    // General case
    if (o1 != o2 && o3 != o4)
        return true;

    // Special Cases
    // p1, q1 and p2 are colinear and p2 lies on segment p1q1
    if (o1 == 0 && onSegment(p1, p2, q1)) return true;

    // p1, q1 and q2 are colinear and q2 lies on segment p1q1
    if (o2 == 0 && onSegment(p1, q2, q1)) return true;

    // p2, q2 and p1 are colinear and p1 lies on segment p2q2
    if (o3 == 0 && onSegment(p2, p1, q2)) return true;

    // p2, q2 and q1 are colinear and q1 lies on segment p2q2
    if (o4 == 0 && onSegment(p2, q1, q2)) return true;

    // Doesn't fall in any of the above cases
    return false;
}

// Determine if a Vertex is Inside a Defined Boundary (Polygon) in the Model
bool Model_2D::is_inside(const Vertex& t, int q) const {
    if (q >= s.size()) {
        std::cout << "Invalid surface index!" << std::endl;
        return false;
    }

    const Surface& surface = s[q];
    int n = surface.boundary.size();

    // Check if the number of edges is sufficient to form a polygon
    if (n < 3) {
        std::cout << "Wrong surface!!" << std::endl;
        return false;    // Invalid surface, cannot contain the point
    }

    // Extract unique vertices from the surface boundary
    std::set<int> vertexNumbers;
    for (const auto& edge : surface.boundary) {
        vertexNumbers.insert(edge.a.vNo);
        vertexNumbers.insert(edge.b.vNo);
    }

    // Map vertex numbers to their coordinates
    std::map<int, glm::vec2> vertexCoords;
    for (const auto& vert : v) {
        if (vertexNumbers.count(vert.vNo)) {
            glm::vec2 coord;
            switch (direction) {
                case FRONT_VIEW:
                    coord = glm::vec2(vert.x, vert.z);
                    break;
                case TOP_VIEW:
                    coord = glm::vec2(vert.x, vert.y);
                    break;
                case SIDE_VIEW:
                    coord = glm::vec2(vert.y, vert.z);
                    break;
                default:
                    return false;
            }
            vertexCoords[vert.vNo] = coord;
        }
    }

    // Build the polygon by ordering the vertices (assumes convex polygon)
    std::vector<glm::vec2> poly;
    for (const auto& vno : vertexNumbers) {
        poly.push_back(vertexCoords[vno]);
    }

    // Use the centroid to sort the vertices in order around the center
    glm::vec2 centroid(0.0f, 0.0f);
    for (const auto& point : poly) {
        centroid += point;
    }
    centroid /= static_cast<float>(poly.size());

    std::sort(poly.begin(), poly.end(), [&centroid](const glm::vec2& a, const glm::vec2& b) {
        float angleA = atan2(a.y - centroid.y, a.x - centroid.x);
        float angleB = atan2(b.y - centroid.y, b.x - centroid.x);
        return angleA < angleB;
    });

    // Now perform the point-in-polygon test
    glm::vec2 p;
    switch (direction) {
        case FRONT_VIEW:
            p = glm::vec2(t.x, t.z);
            break;
        case TOP_VIEW:
            p = glm::vec2(t.x, t.y);
            break;
        case SIDE_VIEW:
            p = glm::vec2(t.y, t.z);
            break;
        default:
            return false;
    }

    // Ray Casting Algorithm
    int count = 0;
    int nVertices = poly.size();
    glm::vec2 extreme(INF, p.y);

    for (int i = 0; i < nVertices; ++i) {
        int next = (i + 1) % nVertices;
        if (doIntersect(poly[i], poly[next], p, extreme)) {
            if (orientation(poly[i], p, poly[next]) == 0) {
                return onSegment(poly[i], p, poly[next]);
            }
            count++;
        }
    }

    return (count % 2 == 1);
}

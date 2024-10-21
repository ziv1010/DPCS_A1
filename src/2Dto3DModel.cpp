// 2Dto3DModel.cpp

#include "2Dto3DModel.h"
#include <map>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <set>
#include <limits>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp> // For glm::to_string (debugging purposes)

// Constructor
Model_3D::Model_3D() {
    // Initialize empty vectors for vertices, edges, and surfaces
    v = std::vector<Vertex>();
    e = std::vector<Edge>();
    s = std::vector<Surface>();
}

// Function to generate the model from input file
void Model_3D::generate_model(const std::string& inputFilePath) {
    std::ifstream inFile(inputFilePath);
    if (!inFile) {
        std::cerr << "Error opening input file: " << inputFilePath << std::endl;
        return;
    }

    int n;
    inFile >> n;
    if (n == 2) {
        std::cout << "File is 2D" << std::endl;

        // Read Top View
        inFile >> topView.direction;
        std::cout << "Top View Direction: " << topView.direction << std::endl;
        topView.read_from_file(inFile);

        // Read Front View
        inFile >> frontView.direction;
        std::cout << "Front View Direction: " << frontView.direction << std::endl;
        frontView.read_from_file(inFile);

        // Read Side View
        inFile >> sideView.direction;
        std::cout << "Side View Direction: " << sideView.direction << std::endl;
        sideView.read_from_file(inFile);

        // Form the wireframe
        form_wireframe();
    } else if (n == 3) {
        std::cout << "File is 3D" << std::endl;
        // Handle 3D input (not covered here)
    } else {
        std::cerr << "Invalid value of n: " << n << std::endl;
    }

    inFile.close();
}

// Function: basic_proj(direction)
// Develop a basic 2D projection of the given 3D solid without any hidden lines.
Model_2D Model_3D::basic_proj(int direction) {
    Model_2D t_view;
    t_view.direction = direction;

    // Project vertices
    for (const auto& vert : v) {
        if (vert.isTrue) {
            Vertex proj_vertex;
            switch (direction) {
                case 0: // Top view (x, y)
                    proj_vertex = Vertex(vert.x, vert.y, 0.0f, vert.vNo);
                    break;
                case 1: // Front view (x, z)
                    proj_vertex = Vertex(vert.x, 0.0f, vert.z, vert.vNo);
                    break;
                case 2: // Side view (y, z)
                    proj_vertex = Vertex(0.0f, vert.y, vert.z, vert.vNo);
                    break;
                default:
                    // Do nothing for invalid direction
                    break;
            }
            t_view.v.push_back(proj_vertex);
        }
    }

    // Project edges
    for (const auto& edge : e) {
        Edge e_temp(
            t_view.v[edge.a.vNo - 1],
            t_view.v[edge.b.vNo - 1],
            edge.eno
        );
        t_view.e.push_back(e_temp);
    }

    // Project surfaces
    for (const auto& surface : s) {
        Surface s_temp;
        for (const auto& boundary_edge : surface.boundary) {
            s_temp.boundary.push_back(t_view.e[boundary_edge.eno - 1]);
        }
        s_temp.sno = surface.sno;
        t_view.s.push_back(s_temp);
    }

    return t_view;
}

// Function: calc_2d(direction)
// Perform segmentation of edges in 2D views if two edges intersect in 2D view.
void Model_3D::calc_2d(int direction) {
    Model_2D* view = nullptr;
    if (direction == 0) {
        view = &topView;
    } else if (direction == 1) {
        view = &frontView;
    } else if (direction == 2) {
        view = &sideView;
    } else {
        std::cout << "Wrong direction input!!" << std::endl;
        return;
    }

    for (size_t i = 0; i < view->e.size(); ++i) {
        if (!view->e[i].isTrue) continue;

        for (size_t j = i + 1; j < view->e.size(); ++j) {
            if (!view->e[j].isTrue) continue;

            Vertex* v_intersect_2D = view->is_intersect(i, j);
            if (v_intersect_2D != nullptr) {
                // Create new vertices v_int1 and v_int2 as copies of v_intersect_2D
                Vertex v_int1 = *v_intersect_2D;
                Vertex v_int2 = *v_intersect_2D;

                // Calculate missing coordinates based on the edges
                if (direction == 0) {
                    // Top view
                    if (view->e[i].a.y == view->e[i].b.y) {
                        v_int1.z = (v[view->e[i].a.vNo - 1].z + v[view->e[i].b.vNo - 1].z) / 2.0f;
                    } else {
                        v_int1.z = (v[view->e[i].a.vNo - 1].z + v[view->e[i].b.vNo - 1].z) / 2.0f;
                    }

                    if (view->e[j].a.y == view->e[j].b.y) {
                        v_int2.z = (v[view->e[j].a.vNo - 1].z + v[view->e[j].b.vNo - 1].z) / 2.0f;
                    } else {
                        v_int2.z = (v[view->e[j].a.vNo - 1].z + v[view->e[j].b.vNo - 1].z) / 2.0f;
                    }
                } else if (direction == 1) {
                    // Front view
                    if (view->e[i].a.z == view->e[i].b.z) {
                        v_int1.y = (v[view->e[i].a.vNo - 1].y + v[view->e[i].b.vNo - 1].y) / 2.0f;
                    } else {
                        v_int1.y = (v[view->e[i].a.vNo - 1].y + v[view->e[i].b.vNo - 1].y) / 2.0f;
                    }

                    if (view->e[j].a.z == view->e[j].b.z) {
                        v_int2.y = (v[view->e[j].a.vNo - 1].y + v[view->e[j].b.vNo - 1].y) / 2.0f;
                    } else {
                        v_int2.y = (v[view->e[j].a.vNo - 1].y + v[view->e[j].b.vNo - 1].y) / 2.0f;
                    }
                } else if (direction == 2) {
                    // Side view
                    if (view->e[i].a.z == view->e[i].b.z) {
                        v_int1.x = (v[view->e[i].a.vNo - 1].x + v[view->e[i].b.vNo - 1].x) / 2.0f;
                    } else {
                        v_int1.x = (v[view->e[i].a.vNo - 1].x + v[view->e[i].b.vNo - 1].x) / 2.0f;
                    }

                    if (view->e[j].a.z == view->e[j].b.z) {
                        v_int2.x = (v[view->e[j].a.vNo - 1].x + v[view->e[j].b.vNo - 1].x) / 2.0f;
                    } else {
                        v_int2.x = (v[view->e[j].a.vNo - 1].x + v[view->e[j].b.vNo - 1].x) / 2.0f;
                    }
                }

                // Assign new vertex numbers and mark as false (temporary)
                v_int1.vNo = v.size() + 1;
                v_int1.isTrue = false;
                v.push_back(v_int1);

                v_int2.vNo = v.size() + 1;
                v_int2.isTrue = false;
                v.push_back(v_int2);

                // Update view vertices
                view->v[v_int1.vNo - 1] = v_int1;
                view->v[v_int2.vNo - 1] = v_int2;

                // Create new edges replacing the original edges
                Edge pe1(view->e[i].a, v_int1, view->e.size() + 1);
                pe1.isTrue = false;
                view->e.push_back(pe1);

                Edge pe2(v_int1, view->e[i].b, view->e.size() + 1);
                pe2.isTrue = false;
                view->e.push_back(pe2);

                Edge qe1(view->e[j].a, v_int2, view->e.size() + 1);
                qe1.isTrue = false;
                view->e.push_back(qe1);

                Edge qe2(v_int2, view->e[j].b, view->e.size() + 1);
                qe2.isTrue = false;
                view->e.push_back(qe2);

                // Mark the original edges as inactive
                view->e[i].isTrue = false;
                view->e[j].isTrue = false;
            }
        }
    }

    calc_2d_1(direction);
}

// Function: calc_2d_1(direction)
// Perform the point-in-polygon algorithm on all segmented edges over the surfaces to determine hidden lines.
void Model_3D::calc_2d_1(int direction) {
    Model_2D* view = nullptr;
    if (direction == 0) {
        view = &topView;
    } else if (direction == 1) {
        view = &frontView;
    } else if (direction == 2) {
        view = &sideView;
    } else {
        return;
    }

    for (auto& edge : view->e) {
        if (!edge.isTrue) continue;

        for (const auto& surface : s) {
            float coefficient = 0.0f;
            switch (direction) {
                case 0:
                    coefficient = surface.coeff[2];
                    break;
                case 1:
                    coefficient = surface.coeff[1];
                    break;
                case 2:
                    coefficient = surface.coeff[0];
                    break;
                default:
                    break;
            }

            if (coefficient != 0) {
                Vertex v_mid = edge.midpoint();

                // Create vertex R as a copy of one of the surface's vertices and modify it
                Vertex R = surface.boundary[0].a;
                switch (direction) {
                    case 0:
                        R.z += 10.0f;
                        break;
                    case 1:
                        R.y += 10.0f;
                        break;
                    case 2:
                        R.x += 10.0f;
                        break;
                    default:
                        break;
                }

                float result = surface.calcProj(v_mid) * surface.calcProj(R);

                if (result < 0) {
                    bool res_inside = view->is_inside(v_mid, surface.sno - 1);
                    if (res_inside) {
                        edge.hidden = true;
                        break; // Break out of the surfaces loop
                    }
                }
            }
        }
    }

     for (auto& edge : view->e) {
        if (!edge.isTrue) continue;

        for (size_t j = 0; j < s.size(); ++j) {
            const auto& surface = s[j];
            float coefficient = 0.0f;
            switch (direction) {
                case 0:
                    coefficient = surface.coeff[2];
                    break;
                case 1:
                    coefficient = surface.coeff[1];
                    break;
                case 2:
                    coefficient = surface.coeff[0];
                    break;
                default:
                    break;
            }
            const float EPSILON = 1e-6f;
            if (std::abs(coefficient) > EPSILON) {
                Vertex v_mid = edge.midpoint();
                Vertex R = surface.boundary[0].a;

                switch (direction) {
                    case 0:
                        R.z += 10.0f;
                        break;
                    case 1:
                        R.y += 10.0f;
                        break;
                    case 2:
                        R.x += 10.0f;
                        break;
                    default:
                        break;
                }

                float result = surface.calcProj(v_mid) * surface.calcProj(R);

                if (result < 0) {
                    bool res_inside = view->is_inside(v_mid, surface.sno - 1);

                    // Debug statements
                    std::cout << "Edge " << edge.eno << " is inside surface " << surface.sno << std::endl;
                    std::cout << "Midpoint: (" << v_mid.x << ", " << v_mid.y << ", " << v_mid.z << ")" << std::endl;

                    if (res_inside) {
                        edge.hidden = true;
                        std::cout << "Edge " << edge.eno << " is hidden." << std::endl;
                        break;
                    } else {
                        std::cout << "Edge " << edge.eno << " is visible." << std::endl;
                    }
                }
            }
        }
    }
}



void Model_3D::form_wireframe() {
    // Build maps from vNo to Vertex for each view
    std::map<int, Vertex> topViewVertices;
    std::map<int, Vertex> frontViewVertices;
    std::map<int, Vertex> sideViewVertices;

    for (const auto& vert : topView.v) {
        topViewVertices[vert.vNo] = vert;
    }

    for (const auto& vert : frontView.v) {
        frontViewVertices[vert.vNo] = vert;
    }

    for (const auto& vert : sideView.v) {
        sideViewVertices[vert.vNo] = vert;
    }

    // Debug: Print the number of vertices in each map
    std::cout << "Top View Vertices: " << topViewVertices.size() << std::endl;
    std::cout << "Front View Vertices: " << frontViewVertices.size() << std::endl;
    std::cout << "Side View Vertices: " << sideViewVertices.size() << std::endl;

    // Build 3D vertices by combining coordinates from different views
    for (const auto& [vno, topVert] : topViewVertices) {
        float x_top = topVert.x; // x from Top View
        float y_top = topVert.y; // y from Top View

        // Debug: Print current top view vertex details
        std::cout << "Top View Vertex vNo: " << vno << " | x: " << x_top << " | y: " << y_top << std::endl;

        float x_front = 0.0f, z_front = 0.0f;
        bool has_front = false;
        if (frontViewVertices.find(vno) != frontViewVertices.end()) {
            has_front = true;
            x_front = frontViewVertices[vno].x;
            z_front = frontViewVertices[vno].z;
            // Debug: Print matching front view vertex details
            std::cout << "Front View Vertex vNo: " << vno << " | x: " << x_front << " | z: " << z_front << std::endl;
        } else {
            std::cerr << "Warning: vNo " << vno << " not found in Front View." << std::endl;
        }
        
        float y_side = 0.0f, z_side = 0.0f;
        bool has_side = false;
        if (sideViewVertices.find(vno) != sideViewVertices.end()) {
            has_side = true;
            y_side = sideViewVertices[vno].y; // y from Side View
            z_side = sideViewVertices[vno].z; // z from Side View
            // Debug: Print matching side view vertex details
            std::cout << "Side View Vertex vNo: " << vno << " | y: " << y_side << " | z: " << z_side << std::endl;
        } else {
            std::cerr << "Warning: vNo " << vno << " not found in Side View." << std::endl;
        }

        // Initialize final coordinates with Top View values
        float x = x_top;
        float y = y_top;
        float z = 0.0f;

        // Check consistency and update x and z from Front View
        if (has_front) {
            if (std::abs(x - x_front) > Model_2D::EPS) {
                std::cerr << "Inconsistent x coordinate for vNo " << vno << " between Top and Front views." << std::endl;
            }
            x = x_front; // Update x from Front View if different
            z = z_front; // z from Front View
        }

        // Check consistency and update y and z from Side View
        if (has_side) {
            if (std::abs(y - y_side) > Model_2D::EPS) {
                std::cerr << "Inconsistent y coordinate for vNo " << vno << " between Top and Side views." << std::endl;
                y = y_side; // Update y from Side View
            }
            // Check z consistency
            if (has_front && std::abs(z - z_side) > Model_2D::EPS) {
                std::cerr << "Inconsistent z coordinate for vNo " << vno << " between Front and Side views." << std::endl;
            } else {
                z = z_side; // Update z from Side View if consistent
            }
        }

        // Debug: Print final computed 3D vertex coordinates
        std::cout << "Final 3D Vertex vNo: " << vno << " | x: " << x << " | y: " << y << " | z: " << z << std::endl;

        // Create 3D Vertex with the computed coordinates
        Vertex v_3d(x, y, z, vno);
        v.push_back(v_3d);
    }

    // Debug: Print the 3D vertices created
    std::cout << "3D Vertices Created:" << std::endl;
    for (const auto& vert : v) {
        std::cout << "vNo: " << vert.vNo << " | x: " << vert.x << " | y: " << vert.y << " | z: " << vert.z << std::endl;
    }

    // Build map from vNo to index in v
    std::map<int, int> vNoToIndex;
    for (size_t i = 0; i < v.size(); ++i) {
        vNoToIndex[v[i].vNo] = static_cast<int>(i);
    }

    // Debug: Print mapping of vNo to index
    std::cout << "vNo to Index Mapping:" << std::endl;
    for (const auto& [vno, idx] : vNoToIndex) {
        std::cout << "vNo: " << vno << " -> Index: " << idx << std::endl;
    }

    // Build 3D edges based on topView edges
    for (const auto& edge : topView.e) {
        int vNoA = edge.a.vNo;
        int vNoB = edge.b.vNo;

        // Debug: Print edge information before checking
        std::cout << "Processing Edge " << edge.eno << ": vNo " << vNoA << " to vNo " << vNoB << std::endl;

        // Check if both vertices exist in the mapping
        if (vNoToIndex.find(vNoA) == vNoToIndex.end() || vNoToIndex.find(vNoB) == vNoToIndex.end()) {
            std::cerr << "Error: Vertex numbers " << vNoA << " or " << vNoB << " not found in 3D vertices." << std::endl;
            continue;
        }

        int idxA = vNoToIndex[vNoA];
        int idxB = vNoToIndex[vNoB];

        Vertex a = v[idxA];
        Vertex b = v[idxB];

        // Avoid adding edges between the same vertex
        if (a.vNo == b.vNo) {
            std::cerr << "Same vertices " << a.vNo << " and " << b.vNo << ". No edge added!!" << std::endl;
            continue;
        }

        // Create Edge
        Edge e_3d(a, b, edge.eno);

        // Use epsilon comparison to check if vertices are actually different
        const float EPSILON = 1e-6f;
        if (std::abs(a.x - b.x) < EPSILON &&
            std::abs(a.y - b.y) < EPSILON &&
            std::abs(a.z - b.z) < EPSILON) {
            std::cerr << "Same vertices " << a.vNo << " and " << b.vNo << " after epsilon check. No edge added!!" << std::endl;
            continue;
        }

        // Debug: Print the newly created edge information
        std::cout << "Edge Created: " << e_3d.eno << " | vNoA: " << e_3d.a.vNo << " | vNoB: " << e_3d.b.vNo << std::endl;

        e_3d.hidden = edge.hidden;  // Preserve hidden flag from 2D view
        e.push_back(e_3d);
    }

    // Debug: Print the number of edges created
    std::cout << "3D Edges Created: " << e.size() << std::endl;

    // Proceed to form 3D faces
    //form_3dfaces();
}


// Function: form_3dfaces()
// Form 3D faces on the wireframe model using the Minimum Surface Angle approach.
void Model_3D::form_3dfaces() {
    int nos = 0; // Number of surfaces
    size_t edge_count = e.size();
    std::vector<std::vector<int>> matrix(edge_count + 2, std::vector<int>(edge_count + 2, 0));
    int ibc = 0;

    for (const auto& vertex : v) {
        // Get neighbours of the vertex
        std::vector<int> neighbours = get_neighbours(vertex);

        for (size_t j = 0; j < neighbours.size(); ++j) {
            Edge* e1 = is_edge(vertex.vNo, neighbours[j]);
            if (!e1) continue;

            int checkno = neighbours[j];

            for (size_t k = j + 1; k < neighbours.size(); ++k) {
                Edge* e2 = is_edge(vertex.vNo, neighbours[k]);
                if (!e2) continue;

                if (matrix[e1->eno][e2->eno] == 1) continue;

                Surface s_temp;
                s_temp.boundary.push_back(*e1);
                s_temp.boundary.push_back(*e2);
                s_temp.calc_coeff();

                // Check for NaN coefficients
                if (std::isnan(s_temp.coeff[0]) || std::isnan(s_temp.coeff[1]) ||
                    std::isnan(s_temp.coeff[2]) || std::isnan(s_temp.coeff[3])) {
                    continue;
                }

                int vno2 = neighbours[k];
                int prev_vno = vertex.vNo;

                bool found = false;
                int loop_count = 0;

                do {
                    std::vector<int> n1 = get_neighbours(v[vno2 - 1]);

                    bool dummy = false;
                    for (size_t l = 0; l < n1.size(); ++l) {
                        if (n1[l] == prev_vno) continue;

                        Edge* e3 = is_edge(vno2, n1[l]);
                        if (!e3) continue;

                        bool var = true;
                        for (const auto& edge_in_surface : s_temp.boundary) {
                            if (e3->overlap(edge_in_surface)) {
                                var = false;
                                break;
                            }
                        }

                        if (s_temp.dotProduct(*e3) == 0 && var) {
                            ++ibc;
                            if (ibc == 1000) {
                                std::cout << "Infinite loop detected. Exiting." << std::endl;
                                return;
                            }

                            s_temp.boundary.push_back(*e3);
                            prev_vno = vno2;
                            vno2 = n1[l];
                            dummy = true;

                            if (s_temp.boundary.size() >= 8) {
                                found = true;
                                break;
                            }

                            if (checkno == n1[l]) {
                                for (const auto& edge_m : s_temp.boundary) {
                                    for (const auto& edge_m1 : s_temp.boundary) {
                                        matrix[edge_m.eno][edge_m1.eno] = 1;
                                        matrix[edge_m1.eno][edge_m.eno] = 1;
                                    }
                                }
                                // After adding a surface to s
                                s_temp.sno = nos + 1;
                                ++nos;

                                // Print surface details
                                std::cout << "Formed Surface " << s_temp.sno << " with edges: ";
                                for (const auto& edge : s_temp.boundary) {
                                    std::cout << edge.eno << " ";
                                }
                                std::cout << std::endl;

                                s.push_back(s_temp);
                                
                                found = true;
                                break;
                            }

                            break; // Break from inner loop
                        }
                    }
                    if (!dummy || found) {
                        break; // Break from do-while loop
                    }
                    ++loop_count;
                    if (loop_count > 1000) {
                        std::cout << "Loop count exceeded. Exiting." << std::endl;
                        return;
                    }
                } while (true);

                if (found) {
                    break; // Break from k loop
                }
            }
        }
    }
}

// Function: get_neighbours(vn)
// Returns vertex numbers of all the neighbours of a given vertex.
std::vector<int> Model_3D::get_neighbours(const Vertex& vn) {
    std::set<int> arr_set;
    for (const auto& edge : e) {
        if (edge.a.vNo == vn.vNo) {
            arr_set.insert(edge.b.vNo);
        } else if (edge.b.vNo == vn.vNo) {
            arr_set.insert(edge.a.vNo);
        }
    }
    return std::vector<int>(arr_set.begin(), arr_set.end());
}

// Function: is_edge(a, b)
// Return the edge corresponding to vertices whose vertex numbers are given as integers.
Edge* Model_3D::is_edge(int a, int b) {
    for (auto& edge : e) {
        if ((edge.a.vNo == a && edge.b.vNo == b) || (edge.a.vNo == b && edge.b.vNo == a)) {
            return &edge;
        }
    }
    return nullptr; // If no such edge is found
}

// Function: model_3d_verification()
// Verify vertices, edges, and surfaces using a set of decision rules and remove false ones.
void Model_3D::model_3d_verification() {
    // Function implementation is empty; intended for future development.
}

// Function: save2Dmodels(out)
// Save the 2D views of the solid in a .txt file.
void Model_3D::save2Dmodels(std::ofstream& out) {
    out << "2" << std::endl;

    for (int direction = 0; direction <= 2; ++direction) {
        out << direction << std::endl;

        Model_2D* view = nullptr;
        if (direction == 0) {
            view = &topView;
        } else if (direction == 1) {
            view = &frontView;
        } else if (direction == 2) {
            view = &sideView;
        }

        out << view->v.size() << std::endl;
        for (const auto& vert : view->v) {
            switch (direction) {
                case 0:
                    out << vert.vNo << " " << vert.x << " " << vert.y << std::endl;
                    break;
                case 1:
                    out << vert.vNo << " " << vert.x << " " << vert.z << std::endl;
                    break;
                case 2:
                    out << vert.vNo << " " << vert.y << " " << vert.z << std::endl;
                    break;
                default:
                    break;
            }
        }

        out << view->e.size() << std::endl;
        for (const auto& edge : view->e) {
            out << edge.eno << " " << edge.a.vNo << " " << edge.b.vNo << " " << (edge.hidden ? 1 : 0) << std::endl;
        }
    }
}

// Function: save3Dmodel(out)
// Save the 3D view of the solid in a .txt file.
void Model_3D::save3Dmodel(std::ofstream& out) {
    out << "3" << std::endl;
    out << v.size() << std::endl;
    for (const auto& vert : v) {
        if (vert.isTrue) {
            out << vert.vNo << " " << vert.x << " " << vert.y << " " << vert.z << std::endl;
        }
    }

    out << e.size() << std::endl;
    for (const auto& edge : e) {
        if (edge.isTrue) {
            out << edge.eno << " " << edge.a.vNo << " " << edge.b.vNo << std::endl;
        }
    }
    out << s.size() << std::endl;
    for (const auto& surface : s) {
        out << surface.sno << " " << surface.boundary.size() << " ";
        for (const auto& edge : surface.boundary) {
            out << edge.eno << " ";
        }
        out << std::endl;
    }
}
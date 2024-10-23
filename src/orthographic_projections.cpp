#include "orthographic_projections.h"
#include "stb_image_write.h"
#include <vector>
#include <cmath>
#include <cfloat>
#include <fstream>
#include <iostream>


void projectToTopView(const Object3D& object, Projection2D& projection) {
    projection.projectedVertices.clear();
    projection.projectedEdges.clear();

    for (const auto& vertex : object.vertices) {
        projection.projectedVertices.push_back({vertex.getX(), vertex.getY()});
    }

    for (const auto& edge : object.edges) {
        projection.projectedEdges.push_back({edge.v1_index, edge.v2_index});
    }
}

void projectToFrontView(const Object3D& object, Projection2D& projection) {
    projection.projectedVertices.clear();
    projection.projectedEdges.clear();

    for (const auto& vertex : object.vertices) {
        projection.projectedVertices.push_back({-vertex.getZ(), vertex.getY()});  // x' = -z, y' = y
    }

    for (const auto& edge : object.edges) {
        projection.projectedEdges.push_back({edge.v1_index, edge.v2_index});
    }
}

void projectToSideView(const Object3D& object, Projection2D& projection) {
    projection.projectedVertices.clear();
    projection.projectedEdges.clear();

    for (const auto& vertex : object.vertices) {
        projection.projectedVertices.push_back({vertex.getX(), -vertex.getZ()});  // x' = x, y' = -z
    }

    for (const auto& edge : object.edges) {
        projection.projectedEdges.push_back({edge.v1_index, edge.v2_index});
    }
}

void saveProjectionsToTextFile(const std::string& filename,
                               const Projection2D& topView,
                               const Projection2D& frontView,
                               const Projection2D& sideView) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return;
    }

    int N = topView.projectedVertices.size();
    outfile << N << std::endl;

    // Output the projected vertices
    for (int i = 0; i < N; ++i) {
        outfile << topView.projectedVertices[i].first << " " << topView.projectedVertices[i].second << std::endl;
        outfile << frontView.projectedVertices[i].first << " " << frontView.projectedVertices[i].second << std::endl;
        outfile << sideView.projectedVertices[i].first << " " << sideView.projectedVertices[i].second << std::endl;
    }

    // Output number of edges in top view
    int Te = topView.projectedEdges.size();
    outfile << Te << std::endl;

    // Output edges in top view
    for (const auto& edge : topView.projectedEdges) {
        int v1 = edge.first;
        int v2 = edge.second;
        float x0 = topView.projectedVertices[v1].first;
        float y0 = topView.projectedVertices[v1].second;
        float x1 = topView.projectedVertices[v2].first;
        float y1 = topView.projectedVertices[v2].second;
        outfile << x0 << " " << y0 << " " << x1 << " " << y1 << std::endl;
    }

    // Repeat for front and side views
    int Fe = frontView.projectedEdges.size();
    outfile << Fe << std::endl;
    for (const auto& edge : frontView.projectedEdges) {
        int v1 = edge.first;
        int v2 = edge.second;
        float x0 = frontView.projectedVertices[v1].first;
        float y0 = frontView.projectedVertices[v1].second;
        float x1 = frontView.projectedVertices[v2].first;
        float y1 = frontView.projectedVertices[v2].second;
        outfile << x0 << " " << y0 << " " << x1 << " " << y1 << std::endl;
    }

    int Se = sideView.projectedEdges.size();
    outfile << Se << std::endl;
    for (const auto& edge : sideView.projectedEdges) {
        int v1 = edge.first;
        int v2 = edge.second;
        float x0 = sideView.projectedVertices[v1].first;
        float y0 = sideView.projectedVertices[v1].second;
        float x1 = sideView.projectedVertices[v2].first;
        float y1 = sideView.projectedVertices[v2].second;
        outfile << x0 << " " << y0 << " " << x1 << " " << y1 << std::endl;
    }

    outfile.close();
}

// Function to scale and draw projections onto one combined PNG
void saveCombinedProjectionAsImage(const std::string& filename,
                                   const Projection2D& topView,
                                   const Projection2D& frontView,
                                   const Projection2D& sideView) {
    const int singleViewWidth = 800;
    const int singleViewHeight = 800;
    const int canvasWidth = singleViewWidth * 3;  // One section for each view
    const int canvasHeight = singleViewHeight;

    // Create a blank white image
    std::vector<unsigned char> image(canvasWidth * canvasHeight * 3, 255); // White background

    auto drawProjection = [&](const Projection2D& projection, int offsetX, int offsetY) {
        // Find the min and max coordinates for scaling
        float minX = FLT_MAX, minY = FLT_MAX, maxX = -FLT_MAX, maxY = -FLT_MAX;

        for (const auto& vertex : projection.projectedVertices) {
            minX = std::min(minX, vertex.first);
            minY = std::min(minY, vertex.second);
            maxX = std::max(maxX, vertex.first);
            maxY = std::max(maxY, vertex.second);
        }

        // Add a small margin for padding
        float margin = 0.05f * std::max(maxX - minX, maxY - minY);
        minX -= margin;
        minY -= margin;
        maxX += margin;
        maxY += margin;

        // Scale and map each vertex into image space
        std::vector<std::pair<int, int>> imageCoords;
        for (const auto& vertex : projection.projectedVertices) {
            int x = static_cast<int>(((vertex.first - minX) / (maxX - minX)) * (singleViewWidth - 1)) + offsetX;
            int y = static_cast<int>(((vertex.second - minY) / (maxY - minY)) * (singleViewHeight - 1)) + offsetY;

            // Invert y to match image coordinate system
            y = singleViewHeight - 1 - y + offsetY;

            // Store the mapped coordinates
            imageCoords.push_back({x, y});

            // Draw the vertex as a small dot
            if (x >= 0 && x < canvasWidth && y >= 0 && y < canvasHeight) {
                int index = (y * canvasWidth + x) * 3;
                image[index] = 0;   // Red
                image[index + 1] = 0; // Green
                image[index + 2] = 0; // Blue
            }
        }

        // Draw the edges
        for (const auto& edge : projection.projectedEdges) {
            int v1 = edge.first;
            int v2 = edge.second;
            int x0 = imageCoords[v1].first;
            int y0 = imageCoords[v1].second;
            int x1 = imageCoords[v2].first;
            int y1 = imageCoords[v2].second;

            // Draw a line from (x0, y0) to (x1, y1)
            int dx = abs(x1 - x0);
            int dy = abs(y1 - y0);
            int sx = x0 < x1 ? 1 : -1;
            int sy = y0 < y1 ? 1 : -1;
            int err = dx - dy;

            int x = x0;
            int y = y0;

            while (true) {
                // Set pixel
                if (x >= 0 && x < canvasWidth && y >= 0 && y < canvasHeight) {
                    int index = (y * canvasWidth + x) * 3;
                    image[index] = 0;   // Red
                    image[index + 1] = 0; // Green
                    image[index + 2] = 0; // Blue
                }

                if (x == x1 && y == y1) break;
                int e2 = 2 * err;
                if (e2 > -dy) { err -= dy; x += sx; }
                if (e2 < dx) { err += dx; y += sy; }
            }
        }
    };

    // Draw top view in the first section
    drawProjection(topView, 0, 0);

    // Draw front view in the second section
    drawProjection(frontView, singleViewWidth, 0);

    // Draw side view in the third section
    drawProjection(sideView, singleViewWidth * 2, 0);

    // Save the image using stb_image_write
    if (stbi_write_png(filename.c_str(), canvasWidth, canvasHeight, 3, image.data(), canvasWidth * 3)) {
        std::cout << "Saved combined image to " << filename << std::endl;
    } else {
        std::cerr << "Error: Failed to save combined image to " << filename << std::endl;
    }
}
#ifndef MARCHINGCUBEGRID_H
#define MARCHINGCUBEGRID_H

#include <thread>
#include <vector>
#include <glm/glm.hpp>

#include "mesh.h"
#include "spatialGrid.h"

class MarchingCubeGrid
{
public:
    struct MarchingCubeVertex
    {
        MarchingCubeVertex(): value(1.0e20)
        {
            points[0] = -1;
            points[1] = -1;
            points[2] = -1;
        }

        double value;
        int points[3];
        int gridIndex;
        glm::vec3 normal;
    };

public:
    MarchingCubeGrid(const double cubeSize, const glm::vec3 minVolume, const glm::vec3 maxVolume);
    ~MarchingCubeGrid();

    void initializeGrid(const double cubeSize, const glm::vec3 minVolume, const glm::vec3 maxVolume);
    void computeIsoValues(std::vector<unsigned int> surfaceVertices, double influenceRadius, SpatialGridPoints spatialGrid);
    void triangulate(Mesh& mesh, std::vector<glm::vec3> pointNormals, bool computeNormals);

    int getNbVertices() {return _resX*_resY*_resZ;}

private:
    unsigned int _nbVertices;
    std::vector<int> _vertices;
    std::vector<MarchingCubeVertex> _verticesData;

    unsigned int _resX;
    unsigned int _resY;
    unsigned int _resZ;
    double _cubeSize;
    glm::vec3 _volMin;
    glm::vec3 _dimensions;

    void updateNormals();
    unsigned int getGridIndex(unsigned int xIndex, unsigned int yIndex, unsigned int zIndex);
    unsigned int getEdgePoint(std::vector<MarchingCubeVertex>, int edgeNo, std::vector<glm::vec3>& points);
    unsigned int getIndex(unsigned int gridIndex, int component);
    glm::vec3 getVertexPosition(unsigned int xIndex, unsigned int yIndex, unsigned int zIndex);
    void setScalarValue(unsigned int xIndex, unsigned int yIndex, unsigned int zIndex, double value);
    bool hasVertexIndexes(std::vector<int> vertexIndexes);
};

#endif // MARCHINGCUBEGRID_H

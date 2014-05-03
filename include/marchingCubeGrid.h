#ifndef MARCHINGCUBEGRID_H
#define MARCHINGCUBEGRID_H

#include <thread>
#include <vector>
#include <glm/glm.hpp>

#include "mesh.h"
#include "spatialGrid.h"
//#include "computeIsoValues_CUDA.h"

typedef std::shared_ptr<SpatialGridPoints> SpatialGridPtr;

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

//cloudVolume.resolution, cloudVolume.minimum, cloudVolume.maximum

public:
    MarchingCubeGrid();
    MarchingCubeGrid(SpatialGridPtr spatialGrid, CloudVolume volume);
    ~MarchingCubeGrid();

    void initializeGrid(const float cubeSize, const glm::vec3 minVolume, const glm::vec3 maxVolume);
    void computeIsoValues(const std::vector<glm::vec3> points, const std::vector<glm::vec3> normals, float resolution);
    void triangulate(Mesh& mesh);

    int getNbVertices() {return _resX*_resY*_resZ;}

private:
    SpatialGridPtr _spatialGrid;

    unsigned int _nbVertices;
    std::vector<int> _vertices;
    std::vector<MarchingCubeVertex> _verticesData;

    unsigned int _resX;
    unsigned int _resY;
    unsigned int _resZ;
    float _cubeSize;
    glm::vec3 _volMin;
    glm::vec3 _dimensions;

    int _count;

    void updateNormals();
    unsigned int getGridIndex(unsigned int xIndex, unsigned int yIndex, unsigned int zIndex);
    unsigned int getEdgePoint(MarchingCubeVertex& v1, MarchingCubeVertex& v2, MarchingCubeVertex& v3, MarchingCubeVertex& v4,
                              MarchingCubeVertex& v5, MarchingCubeVertex& v6, MarchingCubeVertex& v7, MarchingCubeVertex& v8,
                              int edgeNo, std::vector<glm::vec3>& points, std::vector<glm::vec3>& normals);
    unsigned int getIndex(unsigned int gridIndex, int component);
    glm::vec3 getVertexPosition(unsigned int xIndex, unsigned int yIndex, unsigned int zIndex);
    CloudVolume getCellsInRadius(const glm::vec3 position, double radius);
    void setScalarValue(unsigned int xIndex, unsigned int yIndex, unsigned int zIndex, double value);
    double getScalarValue(unsigned int xIndex, unsigned int yIndex, unsigned int zIndex);
    bool hasVertexIndexes(std::vector<int> vertexIndexes);
    glm::vec3& getNormal(uint xIndex, uint yIndex, uint zIndex);
};

#endif // MARCHINGCUBEGRID_H

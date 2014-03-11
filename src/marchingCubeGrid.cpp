#include "marchingCubeGrid.h"
#include "marchingCubeLookupTable.h"

MarchingCubeGrid::MarchingCubeGrid(const double cubeSize, const glm::vec3 minVolume, const glm::vec3 maxVolume)
{
    initializeGrid(cubeSize, minVolume, maxVolume);
}

MarchingCubeGrid::~MarchingCubeGrid()
{
}

void MarchingCubeGrid::initializeGrid(const double cubeSize, const glm::vec3 minVolume, const glm::vec3 maxVolume)
{
    _resX = static_cast<int>(std::ceil((maxVolume.x-minVolume.x)/cubeSize));
    _resY = static_cast<int>(std::ceil((maxVolume.y-minVolume.y)/cubeSize));
    _resZ = static_cast<int>(std::ceil((maxVolume.z-minVolume.z)/cubeSize));

    _cubeSize = cubeSize;
    _volMin = minVolume;

    _dimensions = glm::vec3(_resX, _resY, _resZ);
    _dimensions *= _cubeSize;

    _nbVertices = _resX*_resY*_resZ;
    _vertices.resize(_nbVertices);

    for (unsigned int i = 0; i < _nbVertices; ++i)
        _vertices[i] = -1;
}

void MarchingCubeGrid::updateNormals()
{
}

unsigned int MarchingCubeGrid::getGridIndex(unsigned int xIndex, unsigned int yIndex, unsigned int zIndex)
{
    return xIndex + (yIndex*_resX) + (zIndex*_resX*_resY);
}

unsigned int MarchingCubeGrid::getIndex(unsigned int gridIndex, int component)
{
    unsigned int index;

    switch (component)
    {
        case 0:
            index = gridIndex % _resX;
        break;
        case 1:
            index = gridIndex % (_resX*_resY) / _resX;
        break;
        case 2:
            index = gridIndex / (_resX*_resY);
        break;
    }

    return index;
}

unsigned int MarchingCubeGrid::getEdgePoint(std::vector<MarchingCubeVertex> vertices, int edgeIndex, std::vector<glm::vec3>& points)
{
    MarchingCubeVertex* v1 = 0x0;
    MarchingCubeVertex* v2 = 0x0;
    int axis = 0;

    if (edgeIndex == 1)
    {
        v1 = &vertices.at(0);
        v2 = &vertices.at(1);
        axis = 0;
    }
    else if (edgeIndex == 2)
    {
        v1 = &vertices.at(1);
        v2 = &vertices.at(2);
        axis = 1;
    }
    else if (edgeIndex == 3)
    {
        v1 = &vertices.at(3);
        v2 = &vertices.at(2);
        axis = 0;
    }
    else if (edgeIndex == 4)
    {
        v1 = &vertices.at(0);
        v2 = &vertices.at(3);
        axis = 1;
    }
    else if (edgeIndex == 5)
    {
        v1 = &vertices.at(4);
        v2 = &vertices.at(5);
        axis = 0;
    }
    else if (edgeIndex == 6)
    {
        v1 = &vertices.at(5);
        v2 = &vertices.at(6);
        axis = 1;
    }
    else if (edgeIndex == 7)
    {
        v1 = &vertices.at(7);
        v2 = &vertices.at(6);
        axis = 0;
    }
    else if (edgeIndex == 8)
    {
        v1 = &vertices.at(4);
        v2 = &vertices.at(7);
        axis = 1;
    }
    else if (edgeIndex == 9)
    {
        v1 = &vertices.at(0);
        v2 = &vertices.at(4);
        axis = 2;
    }
    else if (edgeIndex == 10)
    {
        v1 = &vertices.at(1);
        v2 = &vertices.at(5);
        axis = 2;
    }
    else if (edgeIndex == 11)
    {
        v1 = &vertices.at(3);
        v2 = &vertices.at(7);
        axis = 2;
    }
    else
    {
        v1 = &vertices.at(2);
        v2 = &vertices.at(6);
        axis = 2;
    }

    int pointID = v1->points[axis];
    if (pointID == -1)
    {
        glm::vec3 pos1, pos2, pos;
        pos1 = getVertexPosition(getIndex(v1->gridIndex, 0),
                          getIndex(v1->gridIndex, 1),
                          getIndex(v1->gridIndex, 2));
        pos2 = getVertexPosition(getIndex(v2->gridIndex, 0),
                          getIndex(v2->gridIndex, 1),
                          getIndex(v2->gridIndex, 2));

        double t = (0.0 - v1->value) / (v2->value - v1->value);
        pos1 *= 1.0-t;
        pos2 *= t;
        pos = pos1 + pos2;

        pointID = points.size();
        points.push_back(pos);
        v1->points[axis] = pointID;
    }

    return pointID;
}

bool MarchingCubeGrid::hasVertexIndexes(std::vector<int> vertexIndexes)
{
    for (unsigned int i = 0; i < vertexIndexes.size(); ++i)
        if (vertexIndexes.at(i) < 0)
            return false;

    return true;
}

glm::vec3 MarchingCubeGrid::getVertexPosition(unsigned int xIndex, unsigned int yIndex, unsigned int zIndex)
{
    glm::vec3 position;

    position = glm::vec3(xIndex, yIndex, zIndex);
    position *= _cubeSize;
    position += _volMin;

    return position;
}

// Lorensen1987
void MarchingCubeGrid::triangulate(Mesh& mesh, std::vector<glm::vec3>& normals, bool computeNormals)
{
    std::vector<Mesh::Triangle>& triangles = mesh.triangles();
    std::vector<glm::vec3>&points = mesh.points();

    if (computeNormals)
        updateNormals();

    int nbVerticesData = _verticesData.size();
    for (int i = 0; i < nbVerticesData; ++i)
    {
        MarchingCubeVertex& vertex = _verticesData[i];

        unsigned int xIndex = getIndex(vertex.gridIndex, 0);
        unsigned int yIndex = getIndex(vertex.gridIndex, 1);
        unsigned int zIndex = getIndex(vertex.gridIndex, 2);
        if ((xIndex < (_resX-1)) && (yIndex < (_resY-1)) && (zIndex < (_resZ-1)))
        {
            std::vector<int> vertexIndexes;
            vertexIndexes.push_back(_vertices[getGridIndex(xIndex, yIndex, zIndex)]);
            vertexIndexes.push_back(_vertices[getGridIndex(xIndex+1, yIndex, zIndex)]);
            vertexIndexes.push_back(_vertices[getGridIndex(xIndex+1, yIndex+1, zIndex)]);
            vertexIndexes.push_back(_vertices[getGridIndex(xIndex, yIndex+1, zIndex)]);
            vertexIndexes.push_back(_vertices[getGridIndex(xIndex, yIndex, zIndex+1)]);
            vertexIndexes.push_back(_vertices[getGridIndex(xIndex+1, yIndex, zIndex+1)]);
            vertexIndexes.push_back(_vertices[getGridIndex(xIndex+1, yIndex+1, zIndex+1)]);
            vertexIndexes.push_back(_vertices[getGridIndex(xIndex, yIndex+1, zIndex+1)]);

            if (hasVertexIndexes(vertexIndexes))
            {
                std::vector<MarchingCubeVertex> vertices;
                for (unsigned int i = 0; i < vertexIndexes.size(); ++i)
                    vertices.push_back(_verticesData[vertexIndexes.at(i)]);

                unsigned int cubeIndex = 0;
                for (unsigned int i = 0; vertices.size(); ++i)
                    if (vertices.at(i).value < 0.0)
                        cubeIndex |= (int)(pow(2.0, (double)i));

                if ((cubeIndex != 0) && (cubeIndex != 255))
                {
                    int i = 0;
                    while (MarchingCubeLookupTable::triangleList[cubeIndex][i] != -1)
                    {
                        int p1, p2, p3;
                        p1 = getEdgePoint(vertices,
                                          MarchingCubeLookupTable::triangleList[cubeIndex][i+0]+1,
                                          points);
                        p2 = getEdgePoint(vertices,
                                          MarchingCubeLookupTable::triangleList[cubeIndex][i+1]+1,
                                          points);
                        p3 = getEdgePoint(vertices,
                                          MarchingCubeLookupTable::triangleList[cubeIndex][i+2]+1,
                                          points);

                        triangles.push_back(Mesh::Triangle(p1, p2, p3));

                        i += 3;
                    }
                }
            }
        }
    }
}

#include "marchingCubeGrid.h"
#include "marchingCubeLookupTable.h"

#include <iostream>
#include <array>
#include <math.h>

MarchingCubeGrid::MarchingCubeGrid():
    _resX(0),
    _resY(0),
    _resZ(0),
    _cubeSize(0.0f)
{

}

MarchingCubeGrid::MarchingCubeGrid(SpatialGridPtr spatialGrid, CloudVolume volume):
    _resX(0),
    _resY(0),
    _resZ(0),
    _cubeSize(0.0f)
{
    _spatialGrid = spatialGrid;
    initializeGrid(volume.resolution, volume.minimum, volume.maximum);
}

MarchingCubeGrid::~MarchingCubeGrid()
{
}

void MarchingCubeGrid::initializeGrid(const float cubeSize, const glm::vec3 minVolume, const glm::vec3 maxVolume)
{
    _count = 0;

    _resX = static_cast<int>(ceil((maxVolume.x-minVolume.x)/cubeSize));
    _resY = static_cast<int>(ceil((maxVolume.y-minVolume.y)/cubeSize));
    _resZ = static_cast<int>(ceil((maxVolume.z-minVolume.z)/cubeSize));

    _cubeSize = cubeSize;
    _volMin = minVolume;

    _dimensions = glm::vec3(_resX, _resY, _resZ);
    _dimensions *= _cubeSize;

    _nbVertices = _resX*_resY*_resZ;
    _vertices.resize(_nbVertices);

    for (unsigned int i = 0; i < _nbVertices; ++i)
        _vertices[i] = -1;
}

double MarchingCubeGrid::getScalarValue(unsigned int xIndex, unsigned int yIndex, unsigned int zIndex)
{
    int dataIndex = _vertices[getGridIndex(xIndex, yIndex, zIndex)];
    return (dataIndex > -1) ? _verticesData[dataIndex].value : 1.0e20;
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

unsigned int MarchingCubeGrid::getEdgePoint(MarchingCubeVertex& v1,
                                            MarchingCubeVertex& v2,
                                            MarchingCubeVertex& v3,
                                            MarchingCubeVertex& v4,
                                            MarchingCubeVertex& v5,
                                            MarchingCubeVertex& v6,
                                            MarchingCubeVertex& v7,
                                            MarchingCubeVertex& v8,
                                            int edgeNo,
                                            std::vector<glm::vec3>& points,
                                            std::vector<glm::vec3>& normals)
{
    MarchingCubeVertex* va = 0x0;
    MarchingCubeVertex* vb = 0x0;
    int axis = 0;

    if (edgeNo == 1)
    {
        va = &v1;
        vb = &v2;
        axis = 0;
    }
    else if (edgeNo == 2)
    {
        va = &v2;
        vb = &v3;
        axis = 1;
    }
    else if (edgeNo == 3)
    {
        va = &v4;
        vb = &v3;
        axis = 0;
    }
    else if (edgeNo == 4)
    {
        va = &v1;
        vb = &v4;
        axis = 1;
    }
    else if (edgeNo == 5)
    {
        va = &v5;
        vb = &v6;
        axis = 0;
    }
    else if (edgeNo == 6)
    {
        va = &v6;
        vb = &v7;
        axis = 1;
    }
    else if (edgeNo == 7)
    {
        va = &v8;
        vb = &v7;
        axis = 0;
    }
    else if (edgeNo == 8)
    {
        va = &v5;
        vb = &v8;
        axis = 1;
    }
    else if (edgeNo == 9)
    {
        va = &v1;
        vb = &v5;
        axis = 2;
    }
    else if (edgeNo == 10)
    {
        va = &v2;
        vb = &v6;
        axis = 2;
    }
    else if (edgeNo == 11)
    {
        va = &v4;
        vb = &v8;
        axis = 2;
    }
    else //if (edgeNo == 12)
    {
        va = &v3;
        vb = &v7;
        axis = 2;
    }

    int pointID = va->points[axis];
    if (pointID == -1)
    {
        glm::vec3 posA, posB, pos;
        posA = getVertexPosition(getIndex(va->gridIndex, 0),
                          getIndex(va->gridIndex, 1),
                          getIndex(va->gridIndex, 2));
        posB = getVertexPosition(getIndex(vb->gridIndex, 0),
                          getIndex(vb->gridIndex, 1),
                          getIndex(vb->gridIndex, 2));

        double t = (0.0 - va->value) / (vb->value - va->value);

        posA *= 1.0-t;
        posB *= t;

        pos = posA;
        pos += posB;

        glm::vec3 normalA = getNormal(getIndex(va->gridIndex, 0),
                                      getIndex(va->gridIndex, 1),
                                      getIndex(va->gridIndex, 2));
        normalA *= 1.0-t;

        glm::vec3 normal = getNormal(getIndex(vb->gridIndex, 0),
                                    getIndex(vb->gridIndex, 1),
                                    getIndex(vb->gridIndex, 2));

        normal *= t;
        normal += normalA;
        normal = glm::normalize(normal);

        /*float direction = glm::dot(pos, normal);
        if (direction > 0.0f)
        {*/
            normals.push_back(normal);
            pointID = points.size();
            points.push_back(pos);
            va->points[axis] = pointID;
        /*}
        else
            pointID = -1;*/
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

CloudVolume MarchingCubeGrid::getCellsInRadius(const glm::vec3 position, double radius)
{
    CloudVolume volume;

    int xmin, ymin, zmin;

    xmin = static_cast<int>(floor( ((position.x-_volMin.x)-radius)/_cubeSize  ));
    ymin = static_cast<int>(floor( ((position.y-_volMin.y)-radius)/_cubeSize  ));
    zmin = static_cast<int>(floor( ((position.z-_volMin.z)-radius)/_cubeSize  ));

    volume.maximum.x = static_cast<int>(floor( ((position.x-_volMin.x)+radius)/_cubeSize  ));
    volume.maximum.y = static_cast<int>(floor( ((position.y-_volMin.y)+radius)/_cubeSize  ));
    volume.maximum.z = static_cast<int>(floor( ((position.z-_volMin.z)+radius)/_cubeSize  ));

    volume.minimum.x = (xmin<0) ? 0 : static_cast<unsigned int>(xmin);
    if (volume.maximum.x >= _resX) volume.maximum.x = _resX-1;

    volume.minimum.y = (ymin<0) ? 0 : static_cast<unsigned int>(ymin);
    if (volume.maximum.y >= _resY) volume.maximum.y = _resY-1;

    volume.minimum.z = (zmin<0) ? 0 : static_cast<unsigned int>(zmin);
    if (volume.maximum.z >= _resZ) volume.maximum.z = _resZ-1;

    return volume;
}

glm::vec3 MarchingCubeGrid::getVertexPosition(unsigned int xIndex, unsigned int yIndex, unsigned int zIndex)
{
    glm::vec3 position;

    position = glm::vec3(xIndex, yIndex, zIndex);
    position *= _cubeSize;
    position += _volMin;

    return position;
}

void MarchingCubeGrid::updateNormals()
{
    for (long v = 0; v< getNbVertices(); ++v)
    {
        // Only process vertices with data
        if (_vertices[v] == -1)
        {
            continue;
        }

        unsigned int ix = getIndex(v, 0);
        unsigned int iy = getIndex(v, 1);
        unsigned int iz = getIndex(v, 2);

        double value = getScalarValue(ix, iy, iz);

        glm::vec3& normal = getNormal(ix, iy, iz);

        int previousDataIndex = (ix>0) ? _vertices[getGridIndex(ix-1, iy, iz)] : -1;
        int nextDataIndex = (ix<_resX-1) ? _vertices[getGridIndex(ix+1, iy, iz)] : -1;
        if (previousDataIndex == -1)
        {
            if (nextDataIndex == -1)
                normal.x = 0.0;
            else
                normal.x = 2.0 * (_verticesData[nextDataIndex].value - value);
        }
        else if (nextDataIndex == -1)
            normal.x = 2.0 * (value - _verticesData[previousDataIndex].value);
        else
            normal.x = _verticesData[nextDataIndex].value - _verticesData[previousDataIndex].value;

        previousDataIndex = (iy>0) ? _vertices[getGridIndex(ix, iy-1, iz)] : -1;
        nextDataIndex = (iy<_resY-1) ? _vertices[getGridIndex(ix, iy+1, iz)] : -1;
        if (previousDataIndex == -1)
        {
            if (nextDataIndex == -1)
                normal.y = 0.0;
            else
                normal.y = 2.0 * (_verticesData[nextDataIndex].value - value);
        }
        else if (nextDataIndex == -1)
            normal.y = 2.0 * (value - _verticesData[previousDataIndex].value);
        else
            normal.y = _verticesData[nextDataIndex].value -
                _verticesData[previousDataIndex].value;

        previousDataIndex = (iz>0) ? _vertices[getGridIndex(ix, iy, iz-1)] : -1;
        nextDataIndex = (iz<_resZ-1) ? _vertices[getGridIndex(ix, iy, iz+1)] : -1;
        if (previousDataIndex == -1)
        {
            if (nextDataIndex == -1)
                normal.z = 0.0;
            else
                normal.z = 2.0 * (_verticesData[nextDataIndex].value - value);
        }
        else if (nextDataIndex == -1)
            normal.z = 2.0 * (value - _verticesData[previousDataIndex].value);
        else
            normal.z = _verticesData[nextDataIndex].value - _verticesData[previousDataIndex].value;

        normal = glm::normalize(normal);
    }
}

glm::vec3& MarchingCubeGrid::getNormal(uint xIndex, uint yIndex, uint zIndex)
{
    int gridIndex = getGridIndex(xIndex, yIndex, zIndex);
    int dataIndex = _vertices[gridIndex];

    if (dataIndex < 0)
    {
        dataIndex = _verticesData.size();
        _verticesData.push_back(MarchingCubeVertex());
        _vertices[gridIndex] = dataIndex;
        _verticesData[dataIndex].gridIndex = gridIndex;
    }

    return _verticesData[dataIndex].normal;
}

void MarchingCubeGrid::setScalarValue(unsigned int xIndex, unsigned int yIndex, unsigned int zIndex, double value)
{
    int gridIndex = getGridIndex(xIndex, yIndex, zIndex);
    int dataIndex = _vertices[gridIndex];

    if (dataIndex < 0)
    {
        dataIndex = _verticesData.size();
        _verticesData.push_back(MarchingCubeVertex());
        _vertices[gridIndex] = dataIndex;
        _verticesData[dataIndex].gridIndex = gridIndex;
    }

    _verticesData[dataIndex].value = value;
}

// Lorensen1987
void MarchingCubeGrid::computeIsoValues(const std::vector<glm::vec3> points, const std::vector<glm::vec3> normals, float resolution)
{
    double influenceRadius = resolution * 4.0;
    double influenceRadius2 = influenceRadius*influenceRadius;
    double influenceRadius6 = pow(influenceRadius, 6);

    std::vector<double> sumWj;
    std::vector<glm::vec3> sumRjWj;

    int nbGridVertices = getNbVertices();
    sumWj.resize(nbGridVertices, 0.0);
    sumRjWj.resize(nbGridVertices, glm::vec3(0.0,0.0,0.0));

    int nbPoints = points.size();
    for (int p = 0; p < nbPoints; ++p)
    {
        CloudVolume volume;
        volume = getCellsInRadius(points[p], (resolution * 4.0));

        glm::vec3 vertexPos;
        for (int iz=volume.minimum.z; iz<=volume.maximum.z; ++iz)
        {
            for (int iy=volume.minimum.y; iy<=volume.maximum.y; ++iy)
            {
                for (int ix=volume.minimum.x; ix<=volume.maximum.x; ++ix)
                {
                    unsigned int cellIndex = getGridIndex(ix, iy, iz);
                    /*if (!_spatialGrid->isCellEmpty(cellIndex))
                    {*/
                        vertexPos = getVertexPosition(ix, iy, iz);

                        glm::vec3 delta(vertexPos);
                        delta -= points[p];

                        double dist2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                        if (dist2 < influenceRadius2)
                        {
                            /*glm::vec3 normal = normals.at(p);
                            float distance = glm::dot((vertexPos - points[p]), normal);
                            float epsilon = -0.065f;

                            if (distance > epsilon)
                            {*/
                                double dist = std::sqrt(dist2);
                                double Wj = std::pow((1.0 - std::pow(dist/influenceRadius,2)), 3);

                                glm::vec3 gradWj(delta);
                                gradWj *= -6.0 * std::pow(influenceRadius2-dist2, 2) / influenceRadius6;

                                sumWj[cellIndex] += Wj;

                                sumRjWj[cellIndex].x += points[p].x*Wj;
                                sumRjWj[cellIndex].y += points[p].y*Wj;
                                sumRjWj[cellIndex].z += points[p].z*Wj;
                            //}
                        }
                    //}
                }
            }
        }
    }

    glm::vec3 vertexPos;
    for (int c = 0; c < nbGridVertices; ++c)
    {
        /*if (sumWj[c] > 0.0f)
        {*/
            double isoValue = 1.0;
            unsigned int ix = getIndex(c, 0);
            unsigned int iy = getIndex(c, 1);
            unsigned int iz = getIndex(c, 2);

            vertexPos = getVertexPosition(ix, iy, iz);

            glm::vec3 averagePosition(sumRjWj[c]);
            averagePosition /= sumWj[c];

            glm::vec3 deltaToAverage(vertexPos);
            deltaToAverage -= averagePosition;

            isoValue = sqrt(deltaToAverage.x*deltaToAverage.x +
                                   deltaToAverage.y*deltaToAverage.y +
                                   deltaToAverage.z*deltaToAverage.z);
            isoValue -= resolution;

            setScalarValue(ix, iy, iz, isoValue);
        //}
    }
}

void MarchingCubeGrid::triangulate(Mesh& mesh)
{
    std::vector<Mesh::Triangle>& triangles = mesh.triangles();
    std::vector<glm::vec3>& points = mesh.points();
    std::vector<glm::vec3>& normals = mesh.normals();

    updateNormals();

    int nbVerticesData = _verticesData.size();
    for (int v = 0; v < nbVerticesData; ++v)
    {
        MarchingCubeVertex& vertex = _verticesData[v];

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
                MarchingCubeVertex& vertex1 = _verticesData[vertexIndexes.at(0)];
                MarchingCubeVertex& vertex2 = _verticesData[vertexIndexes.at(1)];
                MarchingCubeVertex& vertex3 = _verticesData[vertexIndexes.at(2)];
                MarchingCubeVertex& vertex4 = _verticesData[vertexIndexes.at(3)];
                MarchingCubeVertex& vertex5 = _verticesData[vertexIndexes.at(4)];
                MarchingCubeVertex& vertex6 = _verticesData[vertexIndexes.at(5)];
                MarchingCubeVertex& vertex7 = _verticesData[vertexIndexes.at(6)];
                MarchingCubeVertex& vertex8 = _verticesData[vertexIndexes.at(7)];

                unsigned int cubeIndex = 0;
                if (vertex1.value<0.0) cubeIndex |= 1;
                if (vertex2.value<0.0) cubeIndex |= 2;
                if (vertex3.value<0.0) cubeIndex |= 4;
                if (vertex4.value<0.0) cubeIndex |= 8;
                if (vertex5.value<0.0) cubeIndex |= 16;
                if (vertex6.value<0.0) cubeIndex |= 32;
                if (vertex7.value<0.0) cubeIndex |= 64;
                if (vertex8.value<0.0) cubeIndex |= 128;

                if ((cubeIndex != 0) && (cubeIndex != 255))
                {
                    int i = 0;
                    while (MarchingCubeLookupTable::triangleList[cubeIndex][i] != -1)
                    {
                        int p1, p2, p3;
                        p1 = getEdgePoint(vertex1, vertex2, vertex3, vertex4, vertex5, vertex6, vertex7, vertex8,
                                          MarchingCubeLookupTable::triangleList[cubeIndex][i+0]+1,
                                          points, normals);
                        p2 = getEdgePoint(vertex1, vertex2, vertex3, vertex4, vertex5, vertex6, vertex7, vertex8,
                                          MarchingCubeLookupTable::triangleList[cubeIndex][i+1]+1,
                                          points, normals);
                        p3 = getEdgePoint(vertex1, vertex2, vertex3, vertex4, vertex5, vertex6, vertex7, vertex8,
                                          MarchingCubeLookupTable::triangleList[cubeIndex][i+2]+1,
                                          points, normals);

                        /*if (p1 != -1 && p2 != -1 && p3 != -1)
                        {
                            glm::vec3 vertex1 = points.at(p1);
                            glm::vec3 vertex2 = points.at(p2);
                            glm::vec3 vertex3 = points.at(p3);

                            glm::vec3 normal = glm::cross((vertex1 - vertex2), (vertex2 - vertex3));
                            normal = glm::normalize(normal);

                            float direction = glm::dot(normal, vertex1);
                            //inverse for winding
                            if (direction > 0.0f)*/
                                triangles.push_back(Mesh::Triangle(p3, p2, p1));
                        //}

                        i += 3;
                    }
                }

            }
        }
    }
}

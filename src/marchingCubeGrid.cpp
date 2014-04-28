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

unsigned int MarchingCubeGrid::getEdgePoint(MarchingCubeVertex& v1,
                                            MarchingCubeVertex& v2,
                                            MarchingCubeVertex& v3,
                                            MarchingCubeVertex& v4,
                                            MarchingCubeVertex& v5,
                                            MarchingCubeVertex& v6,
                                            MarchingCubeVertex& v7,
                                            MarchingCubeVertex& v8,
                                            int edgeNo,
                                            std::vector<glm::vec3>& points)
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

        //std::clog << "posA: [" << posA.x << "," << posA.y << "," << posA.z << "]" << std::endl;
        //std::clog << "posB: [" << posB.x << "," << posB.y << "," << posB.z << "]" << std::endl;

        /*std::vector<SpatialGridPoint*> elements;
        _spatialGrid->getElements(getIndex(va->gridIndex, 0), getIndex(va->gridIndex, 1), getIndex(va->gridIndex, 2), elements);

        if (elements.size() > 0)
        {
            glm::vec3 normal = elements.at(0)->normal;*/

            double t = (0.0 - va->value) / (vb->value - va->value);

            posA *= 1.0-t;
            posB *= t;

            pos = posA;
            pos += posB;

            //float distance = glm::dot(pos, normal);

            //std::clog << "pos: [" << pos.x << "," << pos.y << "," << pos.z << "]" << std::endl;

            /*if (distance > 0.0f)
            {*/
                pointID = points.size();
                points.push_back(pos);
                va->points[axis] = pointID;
            /*}
            else
                pointID = -1;*/
        //}
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
    std::clog << "points: " << points.size() << std::endl;
    for (int p = 0; p < nbPoints; ++p)
    {
        CloudVolume volume;
        volume = getCellsInRadius(points[p], (resolution * 2.0));

        glm::vec3 vertexPos;
        for (int iz=volume.minimum.z; iz<=volume.maximum.z; ++iz)
        {
            for (int iy=volume.minimum.y; iy<=volume.maximum.y; ++iy)
            {
                for (int ix=volume.minimum.x; ix<=volume.maximum.x; ++ix)
                {
                    unsigned int cellIndex = getGridIndex(ix, iy, iz);
                    //vertexPos = getVertexPosition(ix, iy, iz);

                    //std::clog << "vertex: [" << vertexPos.x << "," << vertexPos.y << "," << vertexPos.z << "]" << std::endl;

                    //vertexPos = glm::vec3((ix*_cubeSize)+_volMin.x, (iy*_cubeSize)+_volMin.y, (iz*_cubeSize)+_volMin.z);
                    //vertexPos.x = (ix*_cubeSize)+_volMin.x;

                    //std::cout << "resolution: " << volume.resolution << std::endl;
                    //std::cout << "cubeSize: " << _cubeSize << std::endl;

                    //unsigned int cellIndex = 0;

                    /*float resolutionFactor = 20.0f;
                    float ixf = static_cast<float>(ix) / resolutionFactor;
                    float iyf = static_cast<float>(iy) / resolutionFactor;
                    float izf = static_cast<float>(iz) / resolutionFactor;*/

                    //float ixf = static_cast<float>(ix) * resolution;
                    //float iyf = static_cast<float>(iy) * resolution;
                    //float izf = static_cast<float>(iz) * resolution;

                    //float ixf = static_cast<float>(99); float iyf = static_cast<float>(99); float izf = static_cast<float>(99);

                    //vertexPos = glm::vec3(ixf+_volMin.x, iyf+_volMin.y, izf+_volMin.z);

                    vertexPos = getVertexPosition(ix, iy, iz);

                    //std::clog << "vertexFloat: [" << vertexPos.x << "," << vertexPos.y << "," << vertexPos.z << "]" << std::endl;

                    //std::vector<SpatialGridPoint*> elements;
                    //_spatialGrid->getElements(ix, iy, iz, elements);

                    /*if (elements.size() > 0)
                    {*/
                        //glm::vec3 pos = elements.at(0)->pos;


                            //double distance = sqrt(pow(delta.x, 2.0)+pow(delta.y, 2.0)+pow(delta.z, 2.0));
                            glm::vec3 delta(vertexPos);
                            delta -= points[p];

                            double dist2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                            if (dist2 < influenceRadius2)
                            {
                                //glm::vec3 origin = glm::vec3(0.0f, 0.0f, 0.0f);
                                glm::vec3 normal = normals.at(p);
                                //glm::vec3 delta = points[p] + (normal * resolution);
                                //float distanceDelta = glm::distance(delta, origin);
                                //float distancePoint = glm::distance(vertexPos, origin);
                                //float distance = distanceDelta - distancePoint;

                                float distance = glm::dot((vertexPos - points[p]), normal);

                                //std::clog << "length1: " << vertexPos.length() << std::endl;
                                //std::clog << "length2: " << delta.length() << std::endl;
                                //std::clog << "distance: " << deltaDistance << std::endl;

                                //std::clog << "pos: [" << vertexPos.x << "," << vertexPos.y << "," << vertexPos.z << "]" << std::endl;
                                //std::clog << "normal: [" << normal.x << "," << normal.y << "," << normal.z << "]" << std::endl;
                                //std::clog << "delta: [" << delta.x << "," << delta.y << "," << delta.z << "]" << std::endl;
                                //glm::vec3 pointDirection = vertexPos - points.at(p);

                                //std::clog << "vertexPos: [" << vertexPos.x << "," << vertexPos.y << "," << vertexPos.z << "]" << std::endl;
                                //std::clog << "points[" << p << "]: [" << points[p].x << "," << points[p].y << "," << points[p].z << "]" << std::endl;

                                //float direction = glm::dot(pointDirection, normal);
                                //float epsilon = influenceRadius6;//-6.0*pow(influenceRadius2-dist2, 2) / influenceRadius6;

                                //float epsilon = -pow(influenceRadius, influenceRadius);
                                //float epsilon = 0.0f;
                                float epsilon = -influenceRadius6;//0.0f;//-0.02f;

                                /*if (deltaDistance < 0.0f)
                                    std::clog << "test!" << std::endl;*/

                                /*if (direction > -epsilon)
                                {*/
                                if (distance > epsilon)
                                {
                                    double dist = sqrt(dist2);
                                    double Wj = pow((1.0 - pow(dist/influenceRadius,2)), 3);

                                    glm::vec3 gradWj(delta);
                                    gradWj *= -6.0*pow(influenceRadius2-dist2, 2) / influenceRadius6;

                                    sumWj[cellIndex] += Wj;

                                    sumRjWj[cellIndex].x += points[p].x*Wj;
                                    sumRjWj[cellIndex].y += points[p].y*Wj;
                                    sumRjWj[cellIndex].z += points[p].z*Wj;
                                }
                            }
                        }
                    //}
                //}
            }
        }
    }

    glm::vec3 vertexPos;
    for (int c = 0; c < nbGridVertices; ++c)
    {

        if (sumWj[c] > 0.0f)
        {
            double isoValue = 1.0;
            unsigned int ix = getIndex(c, 0);
            unsigned int iy = getIndex(c, 1);
            unsigned int iz = getIndex(c, 2);


            vertexPos = getVertexPosition(ix, iy, iz);

            /*std::vector<SpatialGridPoint*> points;
            _spatialGrid->getElements(ix, iy, iz, points);

            if (points.size() > 0)
            {*/
                //std::clog << "points: " << points.size() << std::endl;

                /*for (int i = 0; i < points.size(); ++i)
                {
                    glm::vec3 point = points.at(i)->pos;
                    glm::vec3 normal = points.at(i)->normal;

    //                std::clog << "point: [" << point.x << "," << point.y << "," << point.z << "]" << std::endl;
    //                std::clog << "normal: [" << normal.x << "," << normal.y << "," << normal.z << "]" << std::endl;
                    //double distance = sqrt(pow(vertexPos.x-point.x,2.0) + pow(vertexPos.y-point.y,2.0) + pow(vertexPos.z-point.z,2.0));

                    glm::vec3 pointVector = vertexPos - point;
                    float direction = glm::dot(pointVector, normal);

                    //std::clog << "direction: " << direction << std::endl;
                    if (direction > -1.0f)
                    {*/
                        glm::vec3 averagePosition(sumRjWj[c]);
                        averagePosition /= sumWj[c];

                        glm::vec3 deltaToAverage(vertexPos);
                        deltaToAverage -= averagePosition;

                        isoValue = sqrt(deltaToAverage.x*deltaToAverage.x +
                                               deltaToAverage.y*deltaToAverage.y +
                                               deltaToAverage.z*deltaToAverage.z);
                        isoValue -= resolution;
                    /*}

                }
            }
            else
                isoValue = 1.0;*/
            setScalarValue(ix, iy, iz, isoValue);

        }

    }
}

void MarchingCubeGrid::triangulate(Mesh& mesh)
{
    std::vector<Mesh::Triangle>& triangles = mesh.triangles();
    std::vector<glm::vec3>& points = mesh.points();

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
                                          points);
                        p2 = getEdgePoint(vertex1, vertex2, vertex3, vertex4, vertex5, vertex6, vertex7, vertex8,
                                          MarchingCubeLookupTable::triangleList[cubeIndex][i+1]+1,
                                          points);
                        p3 = getEdgePoint(vertex1, vertex2, vertex3, vertex4, vertex5, vertex6, vertex7, vertex8,
                                          MarchingCubeLookupTable::triangleList[cubeIndex][i+2]+1,
                                          points);

                        //inverse for winding
                        triangles.push_back(Mesh::Triangle(p3, p2, p1));

                        i += 3;
                    }
                }

            }
        }
    }
}

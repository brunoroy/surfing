#include "marchingCubeGrid.h"
#include "marchingCubeLookupTable.h"

#include <iostream>
#include <array>

MarchingCubeGrid::MarchingCubeGrid(const double cubeSize, const glm::dvec3 minVolume, const glm::dvec3 maxVolume)
{
    initializeGrid(cubeSize, minVolume, maxVolume);
}

MarchingCubeGrid::~MarchingCubeGrid()
{
}

void MarchingCubeGrid::initializeGrid(const double cubeSize, const glm::dvec3 minVolume, const glm::dvec3 maxVolume)
{
    _count = 0;

    _resX = static_cast<int>(std::ceil((maxVolume.x-minVolume.x)/cubeSize));
    _resY = static_cast<int>(std::ceil((maxVolume.y-minVolume.y)/cubeSize));
    _resZ = static_cast<int>(std::ceil((maxVolume.z-minVolume.z)/cubeSize));

    _cubeSize = cubeSize;
    _volMin = minVolume;

    _dimensions = glm::dvec3(_resX, _resY, _resZ);
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
                                            std::vector<glm::dvec3>& points)
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
        glm::dvec3 posA, posB, pos;
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

        pointID = points.size();
        points.push_back(pos);
        va->points[axis] = pointID;
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

CloudVolume MarchingCubeGrid::getCellsInRadius(const glm::dvec3 position, double radius)
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

glm::dvec3 MarchingCubeGrid::getVertexPosition(unsigned int xIndex, unsigned int yIndex, unsigned int zIndex)
{
    glm::dvec3 position;

    position = glm::dvec3(xIndex, yIndex, zIndex);
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
void MarchingCubeGrid::computeIsoValues(const std::vector<glm::dvec3> points, double influenceRadius)
{
    double influenceRadius2 = influenceRadius*influenceRadius;
    double influenceRadius6 = pow(influenceRadius, 6);

    std::vector<double> sumWj;
    std::vector<glm::dmat3x3> sumRjGradWjT;
    std::vector<glm::dvec3> sumGradWj;
    std::vector<glm::dvec3> sumRjWj;

    int nbGridVertices = getNbVertices();
    sumWj.resize(nbGridVertices, 0.0);
    sumRjGradWjT.resize(nbGridVertices, glm::dmat3x3(0.0));
    sumGradWj.resize(nbGridVertices, glm::dvec3(0.0,0.0,0.0));
    sumRjWj.resize(nbGridVertices, glm::dvec3(0.0,0.0,0.0));

    int nbPoints = points.size();
    for (int p = 0; p < nbPoints; ++p)
    {
        CloudVolume volume;
        volume = getCellsInRadius(points[p], influenceRadius);

        glm::dvec3 vertexPos;
        for (int iz=volume.minimum.z; iz<=volume.maximum.z; ++iz)
        {
            for (int iy=volume.minimum.y; iy<=volume.maximum.y; ++iy)
            {
                for (int ix=volume.minimum.x; ix<=volume.maximum.x; ++ix)
                {
                    unsigned int cellIndex = getGridIndex(ix, iy, iz);
                    vertexPos = getVertexPosition(ix, iy, iz);

                    glm::dvec3 delta(vertexPos);
                    delta -= points[p];

                    double dist2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                    if (dist2 < influenceRadius2)
                    {
                        double dist = sqrt(dist2);
                        double Wj = pow((1.0 - pow(dist/influenceRadius,2)), 3);

                        glm::dvec3 gradWj(delta);
                        gradWj *= -6.0*pow(influenceRadius2-dist2, 2) / influenceRadius6;

                        sumWj[cellIndex] += Wj;

                        sumRjGradWjT[cellIndex][0][0] += points[p].x*gradWj.x;
                        sumRjGradWjT[cellIndex][1][0] += points[p].x*gradWj.y;
                        sumRjGradWjT[cellIndex][2][0] += points[p].x*gradWj.z;
                        sumRjGradWjT[cellIndex][0][1] += points[p].y*gradWj.x;
                        sumRjGradWjT[cellIndex][1][1] += points[p].y*gradWj.y;
                        sumRjGradWjT[cellIndex][2][1] += points[p].y*gradWj.z;
                        sumRjGradWjT[cellIndex][0][2] += points[p].z*gradWj.x;
                        sumRjGradWjT[cellIndex][1][2] += points[p].z*gradWj.y;
                        sumRjGradWjT[cellIndex][2][2] += points[p].z*gradWj.z;

                        sumGradWj[cellIndex] += gradWj;

                        sumRjWj[cellIndex].x += points[p].x*Wj;
                        sumRjWj[cellIndex].y += points[p].y*Wj;
                        sumRjWj[cellIndex].z += points[p].z*Wj;
                    }
                }
            }
        }
    }

    glm::dvec3 vertexPos;
    for (int c=0; c<nbGridVertices; ++c)
    {
        unsigned int ix = getIndex(c, 0);
        unsigned int iy = getIndex(c, 1);
        unsigned int iz = getIndex(c, 2);

        double isoValue = 1.0;
        if (sumWj[c] > 0.0)
        {
            vertexPos = getVertexPosition(ix, iy, iz);

            glm::dvec3 averagePosition(sumRjWj[c]);
            averagePosition /= sumWj[c];

            glm::dmat3x3 sumGradWjSumRjWjT;
            sumGradWjSumRjWjT[0][0] = sumGradWj[c].x*sumRjWj[c].x;
            sumGradWjSumRjWjT[0][1] = sumGradWj[c].x*sumRjWj[c].y;
            sumGradWjSumRjWjT[0][2] = sumGradWj[c].x*sumRjWj[c].z;
            sumGradWjSumRjWjT[1][0] = sumGradWj[c].y*sumRjWj[c].x;
            sumGradWjSumRjWjT[1][1] = sumGradWj[c].y*sumRjWj[c].y;
            sumGradWjSumRjWjT[1][2] = sumGradWj[c].y*sumRjWj[c].z;
            sumGradWjSumRjWjT[2][0] = sumGradWj[c].z*sumRjWj[c].x;
            sumGradWjSumRjWjT[2][1] = sumGradWj[c].z*sumRjWj[c].y;
            sumGradWjSumRjWjT[2][2] = sumGradWj[c].z*sumRjWj[c].z;

            double apTerm1 = 1.0/sumWj[c];
            glm::mat3x3 apTerm2;
            for (int i = 0; i <= 2; ++i)
                for (int j = 0; j <= 2; ++j)
                    apTerm2[i][j] = apTerm1 * sumRjGradWjT[c][i][j];

            double apTerm3 = 1.0/(sumWj[c]*sumWj[c]);
            glm::mat3x3 gradAvgPosition;
            for (int i = 0; i <= 2; ++i)
                for (int j = 0; j <= 2; ++j)
                {
                    double apTerm2IJ = apTerm2[i][j];
                    double sumIJ = sumGradWjSumRjWjT[i][j];
                    gradAvgPosition[i][j] = apTerm2IJ - (apTerm3 * sumIJ);
                }

            double x[3] = { 1.0, 1.0, 1.0 };
            double newX[3];
            double maxValue, oldMaxValue = std::numeric_limits<double>::max();
            double threshold = 0.00001;
            double error = std::numeric_limits<double>::max();
            for (int i=0; (error > threshold) && i<500; ++i)
            {
                newX[0] = gradAvgPosition[0][0]*x[0] + gradAvgPosition[0][1]*x[1] + gradAvgPosition[0][2]*x[2];
                newX[1] = gradAvgPosition[1][0]*x[0] + gradAvgPosition[1][1]*x[1] + gradAvgPosition[1][2]*x[2];
                newX[2] = gradAvgPosition[2][0]*x[0] + gradAvgPosition[2][1]*x[1] + gradAvgPosition[2][2]*x[2];

                double absNewX0 = fabs(newX[0]);
                double absNewX1 = fabs(newX[1]);
                double absNewX2 = fabs(newX[2]);
                if ( (absNewX0 >= absNewX1) && (absNewX0 >= absNewX2) )
                {
                    maxValue = newX[0];
                }
                else if (absNewX1 >= absNewX2)
                {
                    maxValue = newX[1];
                }
                else
                {
                    maxValue = newX[2];
                }

                if (maxValue==0.0)
                {
                    break;
                }

                x[0] = newX[0] / maxValue;
                x[1] = newX[1] / maxValue;
                x[2] = newX[2] / maxValue;

                if (i>0)
                {
                    error = fabs(maxValue-oldMaxValue);
                    oldMaxValue = maxValue;
                }
                else
                {
                    oldMaxValue = maxValue;
                }
            }

            double EVmax = fabs(maxValue);

            double f = 1.0;
            const double tHigh = 2.0;
            const double tLow = 0.4;
            if (EVmax > tLow)
            {
                if (EVmax > tHigh) EVmax = tHigh;
                double gamma = (tHigh-EVmax) / (tHigh-tLow);
                f = gamma*gamma*gamma - 3.0*gamma*gamma + 3.0*gamma;
            }

            glm::dvec3 deltaToAverage(vertexPos);
            deltaToAverage -= averagePosition;

            isoValue = sqrt(deltaToAverage.x*deltaToAverage.x +
                                   deltaToAverage.y*deltaToAverage.y +
                                   deltaToAverage.z*deltaToAverage.z);
            isoValue -= (influenceRadius/4.0)*f;
        }

        setScalarValue(ix, iy, iz, isoValue);
    }
}

void MarchingCubeGrid::triangulate(Mesh& mesh, std::vector<glm::dvec3> pointNormals, bool computeNormals)
{
    std::vector<Mesh::Triangle>& triangles = mesh.triangles();
    std::vector<glm::dvec3>& points = mesh.points();
    std::vector<glm::dvec3>& normals = mesh.normals();

    if (computeNormals)
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
                                          points);
                        p2 = getEdgePoint(vertex1, vertex2, vertex3, vertex4, vertex5, vertex6, vertex7, vertex8,
                                          MarchingCubeLookupTable::triangleList[cubeIndex][i+1]+1,
                                          points);
                        p3 = getEdgePoint(vertex1, vertex2, vertex3, vertex4, vertex5, vertex6, vertex7, vertex8,
                                          MarchingCubeLookupTable::triangleList[cubeIndex][i+2]+1,
                                          points);

                        triangles.push_back(Mesh::Triangle(p1, p2, p3));

                        normals.push_back(pointNormals.at(i/3));
                        normals.push_back(pointNormals.at(i/3+1));
                        normals.push_back(pointNormals.at(i/3+2));

                        i += 3;
                    }
                }                
            }
        }
    }
}

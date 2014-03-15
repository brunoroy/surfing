#include "marchingCubeGrid.h"
#include "marchingCubeLookupTable.h"

#include <iostream>

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
    //std::cout << "resolution = [" << _resX <<  ", " << _resY << ", " << _resZ << "]" << std::endl;
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
void MarchingCubeGrid::computeIsoValues(std::vector<unsigned int> surfaceVertices, double influenceRadius, SpatialGridPoints spatialGrid)
{
    double influenceRadius2 = influenceRadius*influenceRadius;
    double influenceRadius6 = pow(influenceRadius, 6);
    std::vector<SpatialGridPoint*>	nearbyPoints;

    // Init cells properties used to compute iso value
    std::vector<double> sumWj;			// SUM(Wj)
    std::vector<glm::mat3> sumRjGradWjT;	// SUM(rj*Gradient(Wj)')
    std::vector<glm::vec3> sumGradWj;		// SUM(Gradient(Wj))
    std::vector<glm::vec3> sumRjWj;		// SUM(rj*Wj)

    int nbGridVertices = getNbVertices();
    sumWj.resize(nbGridVertices, 0.0);
    sumRjGradWjT.resize(nbGridVertices, glm::mat3(0));
    sumGradWj.resize(nbGridVertices, glm::vec3(0.0,0.0,0.0));
    sumRjWj.resize(nbGridVertices, glm::vec3(0.0,0.0,0.0));

    // Compute cells isoValues
    glm::vec3 vertexPos;
    for (int c=0; c<surfaceVertices.size(); ++c)
    {
        unsigned int cellIndex = surfaceVertices[c];
        unsigned int xIndex = getIndex(cellIndex, 0);
        unsigned int yIndex = getIndex(cellIndex, 1);
        unsigned int zIndex = getIndex(cellIndex, 2);
        vertexPos = getVertexPosition(xIndex, yIndex, zIndex);

        // Compute contribution of nearby particles
        spatialGrid.getElements(vertexPos, influenceRadius, nearbyPoints);
        std::vector<SpatialGridPoint*>::iterator it = nearbyPoints.begin();
        std::vector<SpatialGridPoint*>::iterator itEnd = nearbyPoints.end();
        for ( ; it!=itEnd; ++it)
        {
            glm::vec3 point = (*it)->pos;

            // Is cell inside influence radius?
            glm::vec3 delta(vertexPos);
            delta -= point;

            double dist2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
            if (dist2 < influenceRadius2)
            {
                // Compute kernel and it's gradient
                double dist = sqrt(dist2);
                double Wj = pow((1.0 - pow(dist/influenceRadius,2)), 3);
                glm::vec3 gradWj(delta);
                gradWj *= -6.0*pow(influenceRadius2-dist2, 2) / influenceRadius6;

                // Update summation terms of cell
                sumWj[cellIndex] += Wj;

                sumRjGradWjT[cellIndex][0][0] += point.x*gradWj.x;
                sumRjGradWjT[cellIndex][1][0] += point.x*gradWj.y;
                sumRjGradWjT[cellIndex][2][0] += point.x*gradWj.z;
                sumRjGradWjT[cellIndex][0][1] += point.y*gradWj.x;
                sumRjGradWjT[cellIndex][1][1] += point.y*gradWj.y;
                sumRjGradWjT[cellIndex][2][1] += point.y*gradWj.z;
                sumRjGradWjT[cellIndex][0][2] += point.z*gradWj.x;
                sumRjGradWjT[cellIndex][1][2] += point.z*gradWj.y;
                sumRjGradWjT[cellIndex][2][2] += point.z*gradWj.z;

                sumGradWj[cellIndex] += gradWj;

                sumRjWj[cellIndex].x += point.x*Wj;
                sumRjWj[cellIndex].y += point.y*Wj;
                sumRjWj[cellIndex].z += point.z*Wj;
            }
        }

        // Compute isoValue based on particles contributions
        double isoValue = 1.0;
        if (sumWj[cellIndex] > 0.0)
        {
            vertexPos = getVertexPosition(xIndex, yIndex, zIndex);

            // Compute average position ( SUM(rj*Wj)/SUM(Wj) )
            glm::vec3 averagePosition(sumRjWj[cellIndex]);
            averagePosition /= sumWj[cellIndex];

            // Compute the gradient of the average position
            // (1/SUM(Wj)) * SUM(rj*gradWj') - (1/SUM(Wj)^2) * SUM(gradWj) * SUM(rj*Wj)'
            glm::mat3 sumGradWjSumRjWjT;	// SUM(gradWj) * SUM(rj*Wj)'
            sumGradWjSumRjWjT[0][0] = sumGradWj[cellIndex].x*sumRjWj[cellIndex].x;
            sumGradWjSumRjWjT[0][1] = sumGradWj[cellIndex].x*sumRjWj[cellIndex].y;
            sumGradWjSumRjWjT[0][2] = sumGradWj[cellIndex].x*sumRjWj[cellIndex].z;
            sumGradWjSumRjWjT[1][0] = sumGradWj[cellIndex].y*sumRjWj[cellIndex].x;
            sumGradWjSumRjWjT[1][1] = sumGradWj[cellIndex].y*sumRjWj[cellIndex].y;
            sumGradWjSumRjWjT[1][2] = sumGradWj[cellIndex].y*sumRjWj[cellIndex].z;
            sumGradWjSumRjWjT[2][0] = sumGradWj[cellIndex].z*sumRjWj[cellIndex].x;
            sumGradWjSumRjWjT[2][1] = sumGradWj[cellIndex].z*sumRjWj[cellIndex].y;
            sumGradWjSumRjWjT[2][2] = sumGradWj[cellIndex].z*sumRjWj[cellIndex].z;

            glm::mat3 gradAvgPosition =
                ((1.0/sumWj[cellIndex]) * sumRjGradWjT[cellIndex])
                - ((1.0/(sumWj[cellIndex]*sumWj[cellIndex])) * sumGradWjSumRjWjT);

            // Find maximum eigenvalue of the gradient using the
            // Power method
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
                // TODO: We could check if maxValue moves away from range [tlow, thigh] and
                // terminate earlier the algorithm! (Smarter, faster!)
            }

            double EVmax = fabs(maxValue);

            // Compute Radius correction based on EVmax
            double f = 1.0;
            const double tHigh = 2.0;
            const double tLow = 0.4;
            if (EVmax > tLow)
            {
                if (EVmax > tHigh) EVmax = tHigh;
                double gamma = (tHigh-EVmax) / (tHigh-tLow);
                f = gamma*gamma*gamma - 3.0*gamma*gamma + 3.0*gamma;
            }

            // Compute isoValue!!! (Finally...)
            glm::vec3 deltaToAverage(vertexPos);
            deltaToAverage -= averagePosition;

            isoValue = sqrt(deltaToAverage.x*deltaToAverage.x +
                                   deltaToAverage.y*deltaToAverage.y +
                                   deltaToAverage.z*deltaToAverage.z);
            //isoValue -= particleRadius*f;
        }

        setScalarValue(xIndex, yIndex, zIndex, isoValue);
    }
}

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
                for (unsigned int i = 0; i < vertices.size(); ++i)
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

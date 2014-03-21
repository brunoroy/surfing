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
    _resX = static_cast<int>(std::ceil((maxVolume.x-minVolume.x)/cubeSize));
    _resY = static_cast<int>(std::ceil((maxVolume.y-minVolume.y)/cubeSize));
    _resZ = static_cast<int>(std::ceil((maxVolume.z-minVolume.z)/cubeSize));

    _cubeSize = cubeSize;
    _volMin = minVolume;

    _dimensions = glm::dvec3(_resX, _resY, _resZ);
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

unsigned int MarchingCubeGrid::getEdgePoint(std::vector<MarchingCubeVertex> vertices, int edgeIndex, std::vector<glm::dvec3>& points)
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
        glm::dvec3 pos1, pos2, pos;
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

        std::clog << points.size() << std::endl;
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
//void computeIsoValues(std::vector<unsigned int> surfaceVertices, double influenceRadius, SpatialGridPoints spatialGrid)
void MarchingCubeGrid::computeIsoValues(const std::vector<glm::dvec3> points, double influenceRadius)
{
    double influenceRadius2 = influenceRadius*influenceRadius;
    double influenceRadius6 = pow(influenceRadius, 6);

    // Init cells properties used to compute iso value
    std::vector<double> sumWj;			// SUM(Wj)
    std::vector<glm::dmat3x3> sumRjGradWjT;	// SUM(rj*Gradient(Wj)')
    std::vector<glm::dvec3> sumGradWj;		// SUM(Gradient(Wj))
    std::vector<glm::dvec3> sumRjWj;		// SUM(rj*Wj)

    int nbGridVertices = getNbVertices();
    sumWj.resize(nbGridVertices, 0.0);
    sumRjGradWjT.resize(nbGridVertices, glm::dmat3x3(0.0));
    sumGradWj.resize(nbGridVertices, glm::dvec3(0.0,0.0,0.0));
    sumRjWj.resize(nbGridVertices, glm::dvec3(0.0,0.0,0.0));

    //std::cout << "influenceRadius2: " << influenceRadius2 << std::endl;

    /*std::clog << "test" << std::endl;
    glm::dmat3x3 mat3x3 = glm::dmat3x3(0.0);
    std::cout << "[" << mat3x3[0][0] << " " << mat3x3[0][1] << " " << mat3x3[0][2] << "]" << std::endl;
    std::cout << "[" << mat3x3[1][0] << " " << mat3x3[1][1] << " " << mat3x3[1][2] << "]" << std::endl;
    std::cout << "[" << mat3x3[2][0] << " " << mat3x3[2][1] << " " << mat3x3[2][2] << "]" << std::endl;*/

    // Traverse points and add their contribution to nearby cells
    int nbPoints = points.size();
    //std::cout << "nbPoints: " << nbPoints << std::endl;
    for (int p = 0; p < nbPoints; ++p)
    {
        // Get Nearby cells
        //unsigned int minX, maxX, minY, maxY, minZ, maxZ;
        CloudVolume volume;
        volume = getCellsInRadius(points[p], influenceRadius);

        //int c = 0;

        // Process nearby cells
        glm::dvec3 vertexPos;
        for (int iz=volume.minimum.z; iz<=volume.maximum.z; ++iz)
        {
            for (int iy=volume.minimum.y; iy<=volume.maximum.y; ++iy)
            {
                for (int ix=volume.minimum.x; ix<=volume.maximum.x; ++ix)
                {
                    unsigned int cellIndex = getGridIndex(ix, iy, iz);

                    /*if (p == 0)
                        std::cout << "p = [" << points[p].x  << "," << points[p].y << "," << points[p].z << "]" << std::endl;

                    if (p == 0)
                        std::cout << "cellIndex: " << cellIndex << std::endl;*/

                    vertexPos = getVertexPosition(ix, iy, iz);

                    /*if (p == 0)
                        std::cout << "vertexPos: [" << vertexPos.x << "," << vertexPos.y << "," << vertexPos.z << "]" << std::endl;*/

                    // Is cell inside influence radius?
                    glm::dvec3 delta(vertexPos);
                    delta -= points[p];

                    double dist2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                    //if (p == 0)
                    //    std::cout << "delta: [" << delta.x << "," << delta.y << "," << delta.z << "]" << std::endl;

                    /*if (p == 0)
                        std::clog << "dist2: " << dist2 << std::endl;*/

                    if (dist2 < influenceRadius2)
                    {
                        // Compute kernel and it's gradient
                        double dist = sqrt(dist2);
                        double Wj = pow((1.0 - pow(dist/influenceRadius,2)), 3);
                        /*if (p == 0)
                            std::cout << "Wj = " << Wj << std::endl;*/
                        glm::dvec3 gradWj(delta);
                        gradWj *= -6.0*pow(influenceRadius2-dist2, 2) / influenceRadius6;

                        // Update summation terms of cell
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

                        /*if (p == 0)
                        {
                            c++;
                            std::cout << "sumRjGradWjT" << std::endl;
                            std::cout << std::scientific;
                            std::cout << "[" << sumRjGradWjT[cellIndex][0][0] << " " << sumRjGradWjT[cellIndex][0][1] << " " << sumRjGradWjT[cellIndex][0][2] << "]" << std::endl;
                            std::cout << "[" << sumRjGradWjT[cellIndex][1][0] << " " << sumRjGradWjT[cellIndex][1][1] << " " << sumRjGradWjT[cellIndex][1][2] << "]" << std::endl;
                            std::cout << "[" << sumRjGradWjT[cellIndex][2][0] << " " << sumRjGradWjT[cellIndex][2][1] << " " << sumRjGradWjT[cellIndex][2][2] << "]" << std::endl;
                        }*/

                        sumGradWj[cellIndex] += gradWj;

                        sumRjWj[cellIndex].x += points[p].x*Wj;
                        sumRjWj[cellIndex].y += points[p].y*Wj;
                        sumRjWj[cellIndex].z += points[p].z*Wj;


                    }
                }
            }
        }

        /*if (p == 0)
            std::cout << "test count: " << c << std::endl;*/
    }

    //std::cout << "out!" << std::endl;

    // Compute cells isoValues
    glm::dvec3 vertexPos;
    for (int c=0; c<nbGridVertices; ++c)
    {
        unsigned int ix = getIndex(c, 0);
        unsigned int iy = getIndex(c, 1);
        unsigned int iz = getIndex(c, 2);

        // Make sure there was contribution from at least one particle
        double isoValue = 1.0;
        if (sumWj[c] > 0.0)
        {


            /*if (c < 10000)
                std::cout << "cell = " << c << std::endl;*/

            vertexPos = getVertexPosition(ix, iy, iz);

            // Compute average position ( SUM(rj*Wj)/SUM(Wj) )
            glm::dvec3 averagePosition(sumRjWj[c]);
            averagePosition /= sumWj[c];

            // Compute the gradient of the average position
            // (1/SUM(Wj)) * SUM(rj*gradWj') - (1/SUM(Wj)^2) * SUM(gradWj) * SUM(rj*Wj)'
            glm::dmat3x3 sumGradWjSumRjWjT;	// SUM(gradWj) * SUM(rj*Wj)'
            sumGradWjSumRjWjT[0][0] = sumGradWj[c].x*sumRjWj[c].x;
            sumGradWjSumRjWjT[0][1] = sumGradWj[c].x*sumRjWj[c].y;
            sumGradWjSumRjWjT[0][2] = sumGradWj[c].x*sumRjWj[c].z;
            sumGradWjSumRjWjT[1][0] = sumGradWj[c].y*sumRjWj[c].x;
            sumGradWjSumRjWjT[1][1] = sumGradWj[c].y*sumRjWj[c].y;
            sumGradWjSumRjWjT[1][2] = sumGradWj[c].y*sumRjWj[c].z;
            sumGradWjSumRjWjT[2][0] = sumGradWj[c].z*sumRjWj[c].x;
            sumGradWjSumRjWjT[2][1] = sumGradWj[c].z*sumRjWj[c].y;
            sumGradWjSumRjWjT[2][2] = sumGradWj[c].z*sumRjWj[c].z;

            /*if (c == 4370)
            {
                std::cout << "sumGradWjSumRjWjT" << std::endl;
                std::cout << std::scientific;
                std::cout << "[" << sumGradWjSumRjWjT[0][0] << " " << sumGradWjSumRjWjT[0][1] << " " << sumGradWjSumRjWjT[0][2] << "]" << std::endl;
                std::cout << "[" << sumGradWjSumRjWjT[1][0] << " " << sumGradWjSumRjWjT[1][1] << " " << sumGradWjSumRjWjT[1][2] << "]" << std::endl;
                std::cout << "[" << sumGradWjSumRjWjT[2][0] << " " << sumGradWjSumRjWjT[2][1] << " " << sumGradWjSumRjWjT[2][2] << "]" << std::endl;
            }*/

            double apTerm1 = 1.0/sumWj[c];
            //glm::dmat3x3 apTerm2 = apTerm1 * sumRjGradWjT[c];
            glm::mat3x3 apTerm2;
            for (int i = 0; i <= 2; ++i)
                for (int j = 0; j <= 2; ++j)
                    apTerm2[i][j] = apTerm1 * sumRjGradWjT[c][i][j];

            double apTerm3 = 1.0/(sumWj[c]*sumWj[c]);
            glm::mat3x3 gradAvgPosition;// = apTerm2 - (apTerm3 * sumGradWjSumRjWjT);
            for (int i = 0; i <= 2; ++i)
                for (int j = 0; j <= 2; ++j)
                {
                    double apTerm2IJ = apTerm2[i][j];
                    double sumIJ = sumGradWjSumRjWjT[i][j];
                    gradAvgPosition[i][j] = apTerm2IJ - (apTerm3 * sumIJ);
                    /*if (c == 4370)
                    {
                        std::clog << "apTerm2IJ: " << apTerm2IJ << std::endl;
                        std::clog << "sumIJ: " << sumIJ << std::endl;
                        std::clog << "arrayGradAvgPosition[" << i << "][" << j << "] = " << gradAvgPosition[i][j] << std::endl;
                    }*/
                }

            /*if (c == 4370)
            {
                std::cout << "apTerm1 = " << apTerm1 << std::endl;
                std::cout << "apTerm2" << std::endl;
                std::cout << std::scientific;
                std::cout << "[" << apTerm2[0][0] << " " << apTerm2[0][1] << " " << apTerm2[0][2] << "]" << std::endl;
                std::cout << "[" << apTerm2[1][0] << " " << apTerm2[1][1] << " " << apTerm2[1][2] << "]" << std::endl;
                std::cout << "[" << apTerm2[2][0] << " " << apTerm2[2][1] << " " << apTerm2[2][2] << "]" << std::endl;
                std::cout << "apTerm3 = " << apTerm1 << std::endl;

                std::cout << "sumWj[c] = " << sumWj[c] << std::endl;

                std::cout << "sumRjGradWjT" << std::endl;
                std::cout << std::scientific;
                std::cout << "[" << sumRjGradWjT[c][0][0] << " " << sumRjGradWjT[c][0][1] << " " << sumRjGradWjT[c][0][2] << "]" << std::endl;
                std::cout << "[" << sumRjGradWjT[c][1][0] << " " << sumRjGradWjT[c][1][1] << " " << sumRjGradWjT[c][1][2] << "]" << std::endl;
                std::cout << "[" << sumRjGradWjT[c][2][0] << " " << sumRjGradWjT[c][2][1] << " " << sumRjGradWjT[c][2][2] << "]" << std::endl;

                std::cout << "sumGradWjSumRjWjT" << std::endl;
                std::cout << std::scientific;
                std::cout << "[" << sumGradWjSumRjWjT[0][0] << " " << sumGradWjSumRjWjT[0][1] << " " << sumGradWjSumRjWjT[0][2] << "]" << std::endl;
                std::cout << "[" << sumGradWjSumRjWjT[1][0] << " " << sumGradWjSumRjWjT[1][1] << " " << sumGradWjSumRjWjT[1][2] << "]" << std::endl;
                std::cout << "[" << sumGradWjSumRjWjT[2][0] << " " << sumGradWjSumRjWjT[2][1] << " " << sumGradWjSumRjWjT[2][2] << "]" << std::endl;


                std::cout << "gradAvgPosition" << std::endl;
                std::cout << std::scientific;
                std::cout << "[" << gradAvgPosition[0][0] << " " << gradAvgPosition[0][1] << " " << gradAvgPosition[0][2] << "]" << std::endl;
                std::cout << "[" << gradAvgPosition[1][0] << " " << gradAvgPosition[1][1] << " " << gradAvgPosition[1][2] << "]" << std::endl;
                std::cout << "[" << gradAvgPosition[2][0] << " " << gradAvgPosition[2][1] << " " << gradAvgPosition[2][2] << "]" << std::endl;
            }*/

            // Find maximum eigenvalue of the gradient using the
            // Power method
            double x[3] = { 1.0, 1.0, 1.0 };
            double newX[3];
            double maxValue, oldMaxValue = std::numeric_limits<double>::max();
            double threshold = 0.00001;
            double error = std::numeric_limits<double>::max();
            for (int i=0; (error > threshold) && i<500; ++i)
            {
                //newX[0] = gradAvgPosition[0][0]*x[0] + gradAvgPosition[0][1]*x[1] + gradAvgPosition[0][2]*x[2];
                //newX[1] = gradAvgPosition[1][0]*x[0] + gradAvgPosition[1][1]*x[1] + gradAvgPosition[1][2]*x[2];
                //newX[2] = gradAvgPosition[2][0]*x[0] + gradAvgPosition[2][1]*x[1] + gradAvgPosition[2][2]*x[2];

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

                /*if (c == 4370 && i == 499)
                {
                    std::cout << "x = [" << x[0] << "," << x[1] << "," << x[2] << "]" << std::endl;
                    std::cout << "i = " << i << std::endl;
                    std::cout << "error = " << error << std::endl;
                    std::cout << "threshold = " << threshold << std::endl;
                }*/

                //std::cout << "I'm in!" << std::endl;
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
            glm::dvec3 deltaToAverage(vertexPos);
            deltaToAverage -= averagePosition;

            isoValue = sqrt(deltaToAverage.x*deltaToAverage.x +
                                   deltaToAverage.y*deltaToAverage.y +
                                   deltaToAverage.z*deltaToAverage.z);
            isoValue -= (influenceRadius/4.0)*f;

            /*if (c < 100)
                std::clog << "isovalue: " << isoValue << std::endl;*/
        }

        // Set value
        setScalarValue(ix, iy, iz, isoValue);
    }
}
/*void MarchingCubeGrid::computeIsoValues(std::vector<unsigned int> surfaceVertices, double influenceRadius, SpatialGridPoints spatialGrid)
{
    double pointRadius = influenceRadius / 5.0;
    double influenceRadius2 = influenceRadius*influenceRadius;
    double influenceRadius6 = pow(influenceRadius, 6);
    std::vector<SpatialGridPoint*>	nearbyPoints;

    // Init cells properties used to compute iso value
    std::vector<double> sumWj;			// SUM(Wj)
    std::vector<glm::mat3> sumRjGradWjT;	// SUM(rj*Gradient(Wj)')
    std::vector<glm::dvec3> sumGradWj;		// SUM(Gradient(Wj))
    std::vector<glm::dvec3> sumRjWj;		// SUM(rj*Wj)

    int nbGridVertices = getNbVertices();
    sumWj.resize(nbGridVertices, 0.0);
    sumRjGradWjT.resize(nbGridVertices, glm::mat3(0));
    sumGradWj.resize(nbGridVertices, glm::dvec3(0.0,0.0,0.0));
    sumRjWj.resize(nbGridVertices, glm::dvec3(0.0,0.0,0.0));

    // Compute cells isoValues
    glm::dvec3 vertexPos;
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
            glm::dvec3 point = (*it)->pos;

            // Is cell inside influence radius?
            glm::dvec3 delta(vertexPos);
            delta -= point;

            double dist2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
            if (dist2 < influenceRadius2)
            {
                // Compute kernel and it's gradient
                double dist = sqrt(dist2);
                double Wj = pow((1.0 - pow(dist/influenceRadius,2)), 3);
                glm::dvec3 gradWj(delta);
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
            glm::dvec3 averagePosition(sumRjWj[cellIndex]);
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
            glm::dvec3 deltaToAverage(vertexPos);
            deltaToAverage -= averagePosition;

            isoValue = sqrt(deltaToAverage.x*deltaToAverage.x +
                                   deltaToAverage.y*deltaToAverage.y +
                                   deltaToAverage.z*deltaToAverage.z);
            isoValue -= pointRadius*f;
        }

        setScalarValue(xIndex, yIndex, zIndex, isoValue);
    }
}*/

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

                        normals.push_back(pointNormals.at(i/3));
                        normals.push_back(pointNormals.at(i/3+1));
                        normals.push_back(pointNormals.at(i/3+2));

                        i += 3;
                    }
                }                
            }
        }
        //normals.push_back(pointNormals.at(v));
    }

    std::cout << "nbPoints: " << points.size() << std::endl;
}

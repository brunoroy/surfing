#include "marchingCubeGrid.h"
#include "computeIsoValues_CUDA.h"

void MarchingCubeGrid::computeIsoValues(const std::vector<glm::vec3> points, double resolution)
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
        volume = getCellsInRadius(points[p], influenceRadius);

        glm::vec3 vertexPos;
        for (int iz=volume.minimum.z; iz<=volume.maximum.z; ++iz)
        {
            for (int iy=volume.minimum.y; iy<=volume.maximum.y; ++iy)
            {
                for (int ix=volume.minimum.x; ix<=volume.maximum.x; ++ix)
                {
                    unsigned int cellIndex = getGridIndex(ix, iy, iz);
                    vertexPos = getVertexPosition(ix, iy, iz);
                    //unsigned int cellIndex = 0;//getGridIndex(ix, iy, iz);
                    //vertexPos = glm::vec3(0.0f,0.0f,0.0f);//getVertexPosition(ix, iy, iz);

                    glm::vec3 delta(vertexPos);
                    delta -= points[p];

                    double dist2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                    if (dist2 < influenceRadius2)
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
        }
    }

    glm::vec3 vertexPos;
    for (int c = 0; c < nbGridVertices; ++c)
    {
        unsigned int ix = getIndex(c, 0);
        unsigned int iy = getIndex(c, 1);
        unsigned int iz = getIndex(c, 2);

        double isoValue = 1.0;
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
    }
}

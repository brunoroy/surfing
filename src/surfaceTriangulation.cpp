#include "surfaceTriangulation.h"

SurfaceTriangulation::SurfaceTriangulation(const CloudVolume cloudVolume, const std::vector<glm::vec3> points, const std::vector<glm::vec3> normals)
{
    _spatialGrid.reset(new SpatialGridPoints(cloudVolume));
    buildSpatialGrid(points, normals);
    _marchingCubeGrid.reset(new MarchingCubeGrid(_spatialGrid, cloudVolume));
}

SurfaceTriangulation::~SurfaceTriangulation()
{
}

void SurfaceTriangulation::buildSpatialGrid(const std::vector<glm::vec3> points, const std::vector<glm::vec3> normals)
{
    int nbPoints = points.size();
    for (int i = 0; i < nbPoints; ++i)
        _spatialGrid->insert(SpatialGridPoint(points.at(i), normals.at(i), i), points.at(i));
}

std::shared_ptr<MarchingCubeGrid> SurfaceTriangulation::getMarchingCubeGrid()
{
    return _marchingCubeGrid;
}

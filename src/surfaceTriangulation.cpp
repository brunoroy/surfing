#include "surfaceTriangulation.h"

SurfaceTriangulation::SurfaceTriangulation(const CloudVolume cloudVolume)
{
    _grid.reset(new MarchingCubeGrid(cloudVolume.resolution, cloudVolume.minimum, cloudVolume.maximum));
}

SurfaceTriangulation::~SurfaceTriangulation()
{
}

void SurfaceTriangulation::triangulate(Mesh& mesh, std::vector<glm::vec3>& normals, bool computeNormals)
{
    _grid->triangulate(mesh, normals, computeNormals);
}

std::shared_ptr<MarchingCubeGrid> SurfaceTriangulation::getMarchingCubeGrid()
{
    return _grid;
}

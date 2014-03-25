#include "surfaceTriangulation.h"

SurfaceTriangulation::SurfaceTriangulation(const CloudVolume cloudVolume)
{
    _grid.reset(new MarchingCubeGrid(cloudVolume.resolution, cloudVolume.minimum, cloudVolume.maximum));
}

SurfaceTriangulation::~SurfaceTriangulation()
{
}

std::shared_ptr<MarchingCubeGrid> SurfaceTriangulation::getMarchingCubeGrid()
{
    return _grid;
}

#ifndef SURFACETRIANGULATION_H
#define SURFACETRIANGULATION_H

#include "cloudVolume.h"
#include "spatialGrid.h"
#include "marchingCubeGrid.h"

class SurfaceTriangulation
{
public:
    SurfaceTriangulation(const CloudVolume cloudVolume, const std::vector<glm::vec3> points, const std::vector<glm::vec3> normals);
    ~SurfaceTriangulation();

    void buildSpatialGrid(const std::vector<glm::vec3> points, const std::vector<glm::vec3> normals);
    std::shared_ptr<MarchingCubeGrid> getMarchingCubeGrid();

private:
    std::shared_ptr<SpatialGridPoints> _spatialGrid;
    std::shared_ptr<MarchingCubeGrid> _marchingCubeGrid;
};

#endif // SURFACETRIANGULATION_H

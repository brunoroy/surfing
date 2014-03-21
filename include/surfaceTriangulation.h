#ifndef SURFACETRIANGULATION_H
#define SURFACETRIANGULATION_H

#include "marchingCubeGrid.h"
#include "cloudVolume.h"

class SurfaceTriangulation
{
public:
    SurfaceTriangulation(const CloudVolume cloudVolume);
    ~SurfaceTriangulation();

    std::shared_ptr<MarchingCubeGrid> getMarchingCubeGrid();
    //void computeIsoValues(std::vector<unsigned int>& surfaceVertices, double influenceRadius, SpatialGrid<SpatialGridPoint> spatialGrid);
    //void triangulate(Mesh& mesh, std::vector<glm::dvec3>& normals, bool computeNormals);

private:
    std::shared_ptr<MarchingCubeGrid> _grid;
};

#endif // SURFACETRIANGULATION_H

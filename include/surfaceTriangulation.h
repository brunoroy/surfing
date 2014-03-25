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

private:
    std::shared_ptr<MarchingCubeGrid> _grid;
};

#endif // SURFACETRIANGULATION_H

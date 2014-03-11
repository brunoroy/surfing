#ifndef SURFACETRIANGULATION_H
#define SURFACETRIANGULATION_H

#include "marchingCubeGrid.h"
#include "cloudVolume.h"

class SurfaceTriangulation
{
public:
    SurfaceTriangulation(const CloudVolume cloudVolume);
    ~SurfaceTriangulation();

    void triangulate(Mesh& mesh, std::vector<glm::vec3>& normals, bool computeNormals);

private:
    std::shared_ptr<MarchingCubeGrid> _grid;
};

#endif // SURFACETRIANGULATION_H

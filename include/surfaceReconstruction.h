#ifndef SURFACERECONSTRUCTION_H
#define SURFACERECONSTRUCTION_H

#include "surfaceTriangulation.h"
#include "spatialGrid.h"
#include "modelReader.h"
#include "optionManager.h"

typedef SpatialGrid<SpatialGridPoint> SpatialGridPoints;

class SurfaceReconstruction
{
public:
    SurfaceReconstruction(int argc, char** argv);
    ~SurfaceReconstruction();

    void buildSpatialGrid(const std::vector<glm::vec3> points);
    void reconstruct();

private:
    std::shared_ptr<SurfaceTriangulation> _surfaceTriangulation;
    std::shared_ptr<OptionManager> _optionManager;
    std::shared_ptr<ModelReader> _modelReader;
    std::shared_ptr<SpatialGridPoints> _spatialGrid;

    CloudVolume getCloudVolume(std::vector<glm::vec3> points);
};

#endif // SURFACERECONSTRUCTION_H

#ifndef SURFACERECONSTRUCTION_H
#define SURFACERECONSTRUCTION_H

#include "surfaceTriangulation.h"
#include "modelReader.h"
#include "optionManager.h"

class SurfaceReconstruction
{
public:
    SurfaceReconstruction(int argc, char** argv);
    ~SurfaceReconstruction();

    void reconstruct();

private:
    std::shared_ptr<SurfaceTriangulation> _surfaceTriangulation;
    std::shared_ptr<OptionManager> _optionManager;
    std::shared_ptr<ModelReader> _modelReader;

    CloudVolume getCloudVolume(std::vector<glm::vec3> points);
};

#endif // SURFACERECONSTRUCTION_H

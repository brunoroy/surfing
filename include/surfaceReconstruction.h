#ifndef SURFACERECONSTRUCTION_H
#define SURFACERECONSTRUCTION_H

#include "surfaceTriangulation.h"
#include "spatialGrid.h"
#include "modelReader.h"
#include "optionManager.h"

#define PLY_HEADER_FIRST_PART "ply\nformat ascii 1.0\nelement vertex "
#define PLY_HEADER_SECOND_PART "property float x\nproperty float y\nproperty float z\nproperty float w\nproperty float nx\nproperty float ny\nproperty float nz"
#define PLY_HEADER_THIRD_PART "element face "
#define PLY_HEADER_LAST_PART "end_header"
#define SPLIT_CHAR " "

class SurfaceReconstruction
{
public:
    SurfaceReconstruction(int argc, char** argv);
    ~SurfaceReconstruction();

    void buildSpatialGrid(const std::vector<glm::dvec3> points);
    std::vector<unsigned int> extractSurfacePoints(const std::vector<glm::dvec3> points);
    void writeMeshOutput(Mesh mesh, const std::string filename);
    void reconstruct();

private:
    std::shared_ptr<SurfaceTriangulation> _surfaceTriangulation;
    std::shared_ptr<OptionManager> _optionManager;
    std::shared_ptr<ModelReader> _modelReader;
    std::shared_ptr<SpatialGridPoints> _spatialGrid;

    CloudVolume getCloudVolume(std::vector<glm::dvec3> points);
    void writeHeaderOutput(std::ofstream& outputFile, const unsigned int nbPoints, const unsigned int nbFaces);
};

#endif // SURFACERECONSTRUCTION_H

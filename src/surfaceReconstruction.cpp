#include "surfaceReconstruction.h"

#include "timer.h"

#include <iostream>
#include <fstream>

static bool verbose = false;
static char* modelPath = "models";
static double resolution = 0.1;
static bool computeNormals = false;

static OptionEntry entries[] =
{
    {"verbose", "v", OPTION_ARG_NONE, &verbose, "enables printing of messages"},
    {"path", "p", OPTION_ARG_STRING, &modelPath, "loading model files from this path (PLY and PCD are supported)"},
    {"resolution", "r", OPTION_ARG_DOUBLE, &resolution, "sets grid cells resolution (default 0.1)"},
    {"normals", "n", OPTION_ARG_NONE, &computeNormals, "enables face normals computing (default false)"},
};

SurfaceReconstruction::SurfaceReconstruction(int argc, char** argv)
{
    _optionManager.reset(new OptionManager(argc, argv));
    _optionManager->setOptionContext("Reconstruct a surface from unorganized points.");
    _optionManager->addUsage("-p <modelPath>");
    _optionManager->addExample("-p model/bunny");
    _optionManager->addExample("-v -n -p model/bunny -r 0.1");
    std::vector<OptionEntry> optionEntries(entries, entries + sizeof(entries)/sizeof(entries[0]));
    _optionManager->setOptionEntries(optionEntries);
}

SurfaceReconstruction::~SurfaceReconstruction()
{
}

CloudVolume SurfaceReconstruction::getCloudVolume(std::vector<glm::vec3> points)
{
    CloudVolume cloudVolume;
    cloudVolume.minimum = points.at(0);
    cloudVolume.maximum = points.at(0);
    cloudVolume.resolution = resolution;

    int nbPoints = points.size();
    for (int i = 0; i < nbPoints; ++i)
    {
        glm::vec3 point(points.at(i));

        if (point.x < cloudVolume.minimum.x)
            cloudVolume.minimum.x = point.x;
        if (point.y < cloudVolume.minimum.y)
            cloudVolume.minimum.y = point.y;
        if (point.z < cloudVolume.minimum.z)
            cloudVolume.minimum.z = point.z;

        if (point.x > cloudVolume.maximum.x)
            cloudVolume.maximum.x = point.x;
        if (point.y > cloudVolume.maximum.y)
            cloudVolume.maximum.y = point.y;
        if (point.z > cloudVolume.maximum.z)
            cloudVolume.maximum.z = point.z;
    }

    glm::vec3 offset(5.0 * resolution, 5.0 * resolution, 5.0 * resolution);
    cloudVolume.minimum -= offset;
    cloudVolume.maximum += offset;

    return cloudVolume;
}

void SurfaceReconstruction::buildSpatialGrid(const std::vector<glm::vec3> points)
{
    int nbPoints = points.size();
    for (int i = 0; i < nbPoints; ++i)
        _spatialGrid->insert(SpatialGridPoint(points.at(i), i), points.at(i));

    //_surfaceTriangulation->computeIsoValues();
}

std::vector<unsigned int> SurfaceReconstruction::extractSurfacePoints(const std::vector<glm::vec3> points)
{
    std::vector<unsigned int> surfacePoints;
    std::vector<SpatialGridPoint*> elements;

    /*std::cout << "Res: " << spatialGrid.getResX() << ", " << spatialGrid.getResY() << ", "
        << spatialGrid.getResZ() << std::endl;
    std::cout << "Cellsize: " << spatialGrid.getCellSize() << std::endl;
    std::cout << "Volume start: " << spatialGrid.getVolumeStart().x << ", " << spatialGrid.getVolumeStart().y
        << ", " << spatialGrid.getVolumeStart().z << std::endl;*/

    // Find non-empty cells (in spatial grid) with empty neighbors
    int resX = _spatialGrid->getResX();
    int resY = _spatialGrid->getResY();
    int resZ = _spatialGrid->getResZ();
    for (int xIndex = 0; xIndex < resX; ++xIndex)
    {
        for (int yIndex = 0; yIndex < resY; ++yIndex)
        {
            for (int zIndex = 0; zIndex < resZ; ++zIndex)
            {
                if (_spatialGrid->isCellEmpty(xIndex, yIndex, zIndex))
                    continue;

                // Check neighbors
                bool hasEmptyNeighbor = false;
                if (xIndex==0 || xIndex==resX-1 || yIndex==0 || yIndex==resY-1 || zIndex==0 || zIndex==resZ-1)
                {
                    hasEmptyNeighbor = true;
                }
                else
                {
                    for (int dx = -1; !hasEmptyNeighbor && dx<=1; ++dx)
                    {
                        for (int dy = -1; !hasEmptyNeighbor && dy<=1; ++dy)
                        {
                            for (int dz = -1; !hasEmptyNeighbor && dz<=1; ++dz)
                            {
                                if (dx == 0 && dy == 0 && dz == 0)
                                    continue;

                                if (_spatialGrid->isCellEmpty(xIndex + dx, yIndex + dy, zIndex + dz))
                                {
                                    hasEmptyNeighbor = true;
                                }
                            }
                        }
                    }
                }

                // Add cell's points to surface points list if it's a surface cell
                if (hasEmptyNeighbor)
                {
                    _spatialGrid->getElements(xIndex, yIndex, zIndex, elements);
                    for (int e = 0;  e < elements.size(); ++e)
                    {
                        //surfacePoints.push_back(elements[e]->id);
                        surfacePoints.push_back(elements[e]->id);
                        //_surfacePoints[elements[e]->id] = 1.0;
                    }
                }
            }
        }
    }

    return surfacePoints;
}

void SurfaceReconstruction::writeHeaderOutput(std::ofstream& outputFile, const unsigned int nbPoints, const unsigned int nbFaces)
{
    outputFile << PLY_HEADER_FIRST_PART << nbPoints << std::endl;
    outputFile << PLY_HEADER_SECOND_PART << std::endl;
    outputFile << PLY_HEADER_THIRD_PART << nbFaces << std::endl;
    outputFile << PLY_HEADER_LAST_PART << std::endl;
}

void SurfaceReconstruction::writeMeshOutput(Mesh mesh, const std::string filename)
{
    std::string meshFilename = std::string(filename).append(".mesh");
    std::ofstream meshFile(meshFilename);

    if (meshFile.is_open())
    {
        writeHeaderOutput(meshFile, mesh.points().size(), mesh.triangles().size());

        for (int i = 0; i < mesh.points().size(); ++i)
        {
            meshFile << mesh.points().at(i).x << SPLIT_CHAR << mesh.points().at(i).y << SPLIT_CHAR << mesh.points().at(i).z << "1";
            meshFile << mesh.normals().at(i).x << SPLIT_CHAR << mesh.normals().at(i).y << SPLIT_CHAR << mesh.normals().at(i).z << std::endl;
        }

        for (int i = 0; i < mesh.triangles().size(); ++i)
        {
            meshFile << "3" << SPLIT_CHAR << mesh.triangles().at(i).v[0] << SPLIT_CHAR << mesh.triangles().at(i).v[1] << SPLIT_CHAR << mesh.triangles().at(i).v[2] << std::endl;
        }

        meshFile.close();
    }
}

void SurfaceReconstruction::reconstruct()
{
    Mesh mesh;

    OptionParserError *error = NULL;
    if (_optionManager->parseOptions(&error))
    {
        _modelReader.reset(new ModelReader());
        bool modelLoaded = _modelReader->readModel(modelPath);

        if (modelLoaded)
        {
            if (verbose)
                std::clog << "model file " << modelPath << " has been loaded." << std::endl;
            std::vector<glm::vec3> points = _modelReader->getPoints();
            std::vector<glm::vec3> normals = _modelReader->getNormals();
            if (verbose)
                std::clog << points.size() << " points have been read." << std::endl;
            CloudVolume cloudVolume = getCloudVolume(points);
            if (verbose)
            {
                std::clog << "volume: minimum[" << cloudVolume.minimum.x << "," << cloudVolume.minimum.y << "," << cloudVolume.minimum.z << "]";
                std::clog << "; maximum[" << cloudVolume.maximum.x << "," << cloudVolume.maximum.y << "," << cloudVolume.maximum.z << "]";
                std::clog << "; resolution(" << cloudVolume.resolution << ")." << std::endl;
            }

            Timer spatialGridTimer(true);
            _spatialGrid.reset(new SpatialGridPoints(cloudVolume));
            buildSpatialGrid(points);
            auto elapsed = spatialGridTimer.elapsed();
            if (verbose)
                std::cout << "spatial grid: " << std::fixed << elapsed.count() << " ms." << std::endl;

            _surfaceTriangulation.reset(new SurfaceTriangulation(cloudVolume));


            std::shared_ptr<MarchingCubeGrid> grid = _surfaceTriangulation->getMarchingCubeGrid();

            Timer extractSurfacePointsTimer(true);
            std::vector<unsigned int> surfacePoints = extractSurfacePoints(points);
            elapsed = extractSurfacePointsTimer.elapsed();
            if (verbose)
                std::cout << "extract surface points: " << std::fixed << elapsed.count() << " ms." << std::endl;

            Timer computerIsoValuesTimer(true);
            double influenceRadius = 5.0 * cloudVolume.resolution;
            grid->computeIsoValues(surfacePoints, influenceRadius, *_spatialGrid.get());
            elapsed = computerIsoValuesTimer.elapsed();
            if (verbose)
                std::cout << "compute isovalues: " << std::fixed << elapsed.count() << " ms." << std::endl;

            Timer triangulateTimer(true);
            _surfaceTriangulation->getMarchingCubeGrid()->triangulate(mesh, normals, false);
            //_surfaceTriangulation->triangulate(mesh, normals, false);
            elapsed = triangulateTimer.elapsed();
            if (verbose)
                std::cout << "triangulation: " << std::fixed << elapsed.count() << " ms." << std::endl;

            Timer writeMeshTimer(true);
            writeMeshOutput(mesh, modelPath);
            elapsed = writeMeshTimer.elapsed();
            if (verbose)
                std::cout << "writing mesh file: " << std::fixed << elapsed.count() << " ms." << std::endl;
        }
        else
            std::clog << "no model has been loaded." << std::endl;
    }
}



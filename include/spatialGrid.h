#ifndef SPATIALGRID_H
#define SPATIALGRID_H

#include <vector>

#include "cloudVolume.h"

template<class T>
class SpatialGrid
{
public:
    SpatialGrid(double cellSize, const CloudVolume volume);
    ~SpatialGrid();

    void initializeGrid(double cellSize, const CloudVolume volume);
    void clear();

    void insert(const T& element, const glm::vec3& position);
    void insert(const T& element, const glm::vec3& AABBmin, const glm::vec3& AABBmax);

    void getElements(const glm::vec3& position, double radius, std::vector<T*>& elements);
    void getElements(int xIndex, int yIndex, int zIndex, std::vector<T*>& elements);

    unsigned int getNbCells() {return _grid.size();}
    unsigned int getResX() {return _resX;}
    unsigned int getResY() {return _resY;}
    unsigned int getResZ() {return _resZ;}
    double getCellSize() {return _cellSize;}
    glm::vec3 getVolumeStart() const {return _volume.minimum;}
    bool isCellEmpty(int xIndex, int yIndex, int zIndex) {return _grid[getGridIndex(xIndex, yIndex, zIndex)].empty();}

private:
    int getGridIndex(int xIndex, int yIndex, int zIndex);
    int getXIndex(const glm::vec3& position);
    int getYIndex(const glm::vec3& position);
    int getZIndex(const glm::vec3& position);
    int getXIndex(double xPos);
    int getYIndex(double yPos);
    int getZIndex(double zPos);

private:
    std::vector<std::vector<T> > _grid;
    unsigned int _resX;
    unsigned int _resY;
    unsigned int _resZ;
    double _cellSize;
    CloudVolume _volume;
};

#endif // SPATIALGRID_H

#include "spatialGrid.h"

template<class T>
SpatialGrid<T>::SpatialGrid(double cellSize, const CloudVolume volume):
    _resX(0),
    _resY(0),
    _resZ(0),
    _cellSize(0.0)
{
    initializeGrid(cellSize, volume);
}

template<class T>
SpatialGrid<T>::~SpatialGrid()
{
    clear();
}

template<class T>
void SpatialGrid<T>::initializeGrid(double cellSize, const CloudVolume volume)
{
    if (!_grid.empty())
        clear();

    _resX = static_cast<int>(ceil((volume.maximum.x-volume.minimum.x)/cellSize));
    _resY = static_cast<int>(ceil((volume.maximum.y-volume.minimum.y)/cellSize));
    _resZ = static_cast<int>(ceil((volume.maximum.z-volume.minimum.z)/cellSize));

    _cellSize = cellSize;
    _volume.minimum = volume.minimum;

    _grid.resize(_resX*_resY*_resZ);
}

template<class T>
void SpatialGrid<T>::clear()
{
    _grid.clear();
    std::vector<std::vector<T> >().swap(_grid);
}

template<class T>
void SpatialGrid<T>::insert(const T& element, const glm::vec3& position)
{
    int xIndex = getXIndex(position);
    int yIndex = getYIndex(position);
    int zIndex = getZIndex(position);
    int cellIndex = getGridIndex(xIndex, yIndex, zIndex);

    if ((xIndex>=0) && (xIndex<_resX) && (yIndex>=0) && (yIndex<_resY) &&
            (zIndex>=0) && (zIndex<_resZ))
    {
        _grid[cellIndex].push_back(element);
    }
}

template<class T>
void SpatialGrid<T>::insert(const T& element, const glm::vec3& AABBmin, const glm::vec3& AABBmax)
{
    int xMin = getXIndex(AABBmin.x);
    int yMin = getYIndex(AABBmin.y);
    int zMin = getZIndex(AABBmin.z);
    int xMax = getXIndex(AABBmax.x);
    int yMax = getYIndex(AABBmax.y);
    int zMax = getZIndex(AABBmax.z);

    if (xMin<0) xMin = 0;
    if (xMax>=_resX) xMax = _resX-1;

    if (yMin<0) yMin = 0;
    if (yMax>=_resY) yMax = _resY-1;

    if (zMin<0) zMin = 0;
    if (zMax>=_resZ) zMax = _resZ-1;

    for (int xIndex = xMin; xIndex <= xMax; ++xIndex)
    {
        for (int yIndex = yMin; yIndex <= yMax; ++yIndex)
        {
            for (int zIndex = zMin; zIndex <= zMax; ++zIndex)
            {
                int cellIndex = getGridIndex(xIndex, yIndex, zIndex);
                _grid[cellIndex].push_back(element);
            }
        }
    }
}

template<class T>
void SpatialGrid<T>::getElements(const glm::vec3& position, double radius, std::vector<T*>&	elements)
{
    int xMin = getXIndex(position.x - radius);
    int yMin = getYIndex(position.y - radius);
    int zMin = getZIndex(position.z - radius);
    int xMax = getXIndex(position.x + radius);
    int yMax = getYIndex(position.y + radius);
    int zMax = getZIndex(position.z + radius);

    if (xMin<0) xMin = 0;
    if (xMax>=_resX) xMax = _resX-1;

    if (yMin<0) yMin = 0;
    if (yMax>=_resY) yMax = _resY-1;

    if (zMin<0) zMin = 0;
    if (zMax>=_resZ) zMax = _resZ-1;

    int e = 0;
    for (int ix=xMin; ix<=xMax; ++ix)
    {
        for (int iy=yMin; iy<=yMax; ++iy)
        {
            for (int iz=zMin; iz<=zMax; ++iz)
            {
                int cellIndex = getGridIndex(ix, iy, iz);

                std::vector<T>& cell = _grid[cellIndex];
                for (int i=0; i<cell.size(); ++i, ++e)
                {
                    if (e < elements.size())
                    {
                        elements[e] = &cell[i];
                    }
                    else
                    {
                        elements.push_back(&cell[i]);
                    }
                }
            }
        }
    }

    if (e < elements.size())
    {
        elements.resize(e);
    }
}

template<class T>
void SpatialGrid<T>::getElements(int xIndex, int yIndex, int zIndex, std::vector<T*>& elements)
{
    int cellIndex = getGridIndex(xIndex, yIndex, zIndex);
    std::vector<T>& cell = _grid[cellIndex];
    int e=0;
    for (int i=0; i<cell.size(); ++i, ++e)
    {
        if (e < elements.size())
        {
            elements[e] = &cell[i];
        }
        else
        {
            elements.push_back(&cell[i]);
        }
    }

    if (e < elements.size())
    {
        elements.resize(e);
    }
}

template<class T>
int SpatialGrid<T>::getGridIndex(int xIndex, int yIndex, int zIndex)
{
    return xIndex + (yIndex*_resX) + (zIndex*_resX*_resY);
}

template<class T>
int SpatialGrid<T>::getXIndex(const glm::vec3& position)
{
    return static_cast<int>(floor( (position.x-_volume.minimum.x)/_cellSize ));
}

template<class T>
int SpatialGrid<T>::getYIndex(const glm::vec3& position)
{
    return static_cast<int>(floor( (position.y-_volume.minimum.y)/_cellSize ));
}

template<class T>
int SpatialGrid<T>::getZIndex(const glm::vec3& position)
{
    return static_cast<int>(floor( (position.z-_volume.minimum.z)/_cellSize ));
}

template<class T>
int SpatialGrid<T>::getXIndex(double xPos)
{
    return static_cast<int>(floor( (xPos-_volume.minimum.x)/_cellSize ));
}

template<class T>
int SpatialGrid<T>::getYIndex(double yPos)
{
    return static_cast<int>(floor( (yPos-_volume.minimum.y)/_cellSize ));
}

template<class T>
int SpatialGrid<T>::getZIndex(double zPos)
{
    return static_cast<int>(floor( (zPos-_volume.minimum.z)/_cellSize ));
}

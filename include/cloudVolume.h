#ifndef CLOUDVOLUME_H
#define CLOUDVOLUME_H

#include <glm/glm.hpp>

struct CloudVolume
{
    double resolution;
    glm::vec3 minimum;
    glm::vec3 maximum;
};

#endif // CLOUDVOLUME_H

#ifndef CLOUDVOLUME_H
#define CLOUDVOLUME_H

#include <glm/glm.hpp>

struct CloudVolume
{
    CloudVolume() {}
    CloudVolume(glm::vec3 minimum, glm::vec3 maximum): minimum(minimum), maximum(maximum) {}

    double resolution;
    glm::vec3 minimum;
    glm::vec3 maximum;
};

#endif // CLOUDVOLUME_H

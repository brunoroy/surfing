#ifndef CLOUDVOLUME_H
#define CLOUDVOLUME_H

#include <glm/glm.hpp>

struct CloudVolume
{
    CloudVolume() {}
    CloudVolume(glm::dvec3 minimum, glm::dvec3 maximum): minimum(minimum), maximum(maximum) {}

    double resolution;
    glm::dvec3 minimum;
    glm::dvec3 maximum;
};

#endif // CLOUDVOLUME_H

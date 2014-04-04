#ifndef COMPUTE_ISO_VALUES_CUDA_H
#define COMPUTE_ISO_VALUES_CUDA_H

#include <vector>
#include <glm/glm.hpp>

void computeIsoValues(const std::vector<glm::vec3> points, double resolution);

#endif	// COMPUTE_ISO_VALUES_CUDA_H

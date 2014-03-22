#ifndef MESH_H
#define MESH_H

#include <vector>

#include <glm/glm.hpp>

class Mesh
{
public:
    struct Triangle
    {
        Triangle() {}
        Triangle(unsigned int v0, unsigned int v1, unsigned int v2)
        {
            v[0] = v0;
            v[1] = v1;
            v[2] = v2;
        }

        unsigned int v[3];
    };

public:
    Mesh();
    ~Mesh();

    std::vector<Triangle>& triangles() {return _triangles;}
    std::vector<glm::dvec3>& points() {return _points;}
    std::vector<glm::dvec3>& normals() {return _normals;}

    const std::vector<Triangle>& triangles() const  {return _triangles;}
    const std::vector<glm::dvec3>& points() const {return _points;}
    const std::vector<glm::dvec3>& normals() const {return _normals;}

private:
    std::vector<Triangle> _triangles;
    std::vector<glm::dvec3> _points;
    std::vector<glm::dvec3> _normals;
};

#endif // MESH_H

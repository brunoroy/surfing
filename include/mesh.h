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
    std::vector<glm::vec3>& points() {return _points;}
    std::vector<glm::vec3>& normals() {return _normals;}
    std::vector<glm::vec3>& colors() {return _colors;}

    const std::vector<Triangle>& triangles() const  {return _triangles;}
    const std::vector<glm::vec3>& points() const {return _points;}
    const std::vector<glm::vec3>& normals() const {return _normals;}
    const std::vector<glm::vec3>& colors() const {return _colors;}

private:
    std::vector<Triangle> _triangles;
    std::vector<glm::vec3> _points;
    std::vector<glm::vec3> _normals;
    std::vector<glm::vec3> _colors;
};

#endif // MESH_H

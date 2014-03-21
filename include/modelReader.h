#ifndef MODELREADER_H
#define MODELREADER_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

#include <glm/glm.hpp>

#define PLY_END_HEADER "end_header"
#define PCD_FILE_EXTENSION "PCD"
#define PLY_FILE_EXTENSION "PLY"

class ModelReader
{
public:
    ModelReader();
    ~ModelReader();

    bool readPCD(const std::string filename);
    bool readPLY(const std::string filename);
    bool isFormatPCD(const std::string filename);
    bool isFormatPLY(const std::string filename);
    std::string getFileExtension(const std::string filename);
    bool readModel(std::string path);

    std::vector<glm::dvec3> getPoints();
    std::vector<glm::dvec3> getNormals();

private:
    std::vector<glm::dvec3> _points;
    std::vector<glm::dvec3> _normals;

    std::vector<std::string> split(const std::string input);
};

#endif // MODELREADER_H

#include "modelReader.h"

#include <sstream>
#include <iterator>

ModelReader::ModelReader()
{
}

ModelReader::~ModelReader()
{

}

bool ModelReader::readModel(std::string path)
{
    if (isFormatPCD(path))
        return readPCD(path);
    else if (isFormatPLY(path))
        return readPLY(path);

    return false;
}

bool ModelReader::readPCD(const std::string filename)
{
    return false;
}

bool ModelReader::readPLY(const std::string filename)
{
    bool isVertexProperties = false;
    std::string line;
    std::ifstream inputFile(filename);
    std::vector<std::string> values;

    if (inputFile.is_open())
    {
        _points.clear();
        while (std::getline(inputFile, line))
        {
            if (line.compare(PLY_END_HEADER) == 0)
            {
                isVertexProperties = true;
            }
            else if (isVertexProperties)
            {
                values = split(line);
                glm::dvec3 point(std::stof(values.at(0)), std::stof(values.at(1)), std::stof(values.at(2)));
                glm::dvec3 normal(std::stof(values.at(4)), std::stof(values.at(5)), std::stof(values.at(6)));
                _points.push_back(point);
                _normals.push_back(normal);
            }
        }
        inputFile.close();
        return true;
    }

    return false;
}

bool ModelReader::isFormatPCD(const std::string filename)
{
    std::string fileExtension = getFileExtension(filename);
    if (fileExtension.compare(PCD_FILE_EXTENSION) == 0)
        return true;
    else
        return false;
}

bool ModelReader::isFormatPLY(const std::string filename)
{
    std::string fileExtension = getFileExtension(filename);
    if (fileExtension.compare(PLY_FILE_EXTENSION) == 0)
        return true;
    else
        return false;
}

std::string ModelReader::getFileExtension(const std::string filename)
{
    std::string fileExtension = filename.substr(filename.find_last_of('.')+1, filename.length()-filename.find_last_of('.'));
    std::transform(fileExtension.begin(), fileExtension.end(), fileExtension.begin(), (int(*)(int))std::toupper);
    return fileExtension;
}

std::vector<std::string> ModelReader::split(const std::string input)
{
    std::istringstream streamInput(input);
    return {std::istream_iterator<std::string>{streamInput}, std::istream_iterator<std::string>{}};
}

std::vector<glm::dvec3> ModelReader::getPoints()
{
    return _points;
}

std::vector<glm::dvec3> ModelReader::getNormals()
{
    return _normals;
}

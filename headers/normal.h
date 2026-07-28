#ifndef NORMAL_H
#define NORMAL_H

#include <glm/glm.hpp>
#include <glm/gtc/epsilon.hpp>

#include <vector>
#include <cmath>
#include <unordered_map>

struct Vec3Hash{
    size_t operator()(const glm::vec3& v) const{
        size_t h1 = std::hash<float>{}(v.x);
        size_t h2 = std::hash<float>{}(v.y);
        size_t h3 = std::hash<float>{}(v.z);

        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

bool operator == (const glm::vec3& a, const glm::vec3& b){
    return glm::all(glm::epsilonEqual(a, b, 1e-6f));
}

std::vector<float> GetNormals(const std::vector<float>& verticies){
    std::vector<float> normals;
    std::vector<glm::vec3> positions;
    const size_t count = verticies.size() / 3;
    positions.reserve(count);

    for (size_t i = 0; i < verticies.size(); i += 3){
        positions.emplace_back(verticies[i],verticies[i+1],verticies[i+2]);
    }

    std::unordered_map<glm::vec3, glm::vec3, Vec3Hash> normalSums;
    std::unordered_map<glm::vec3, int, Vec3Hash> normalCounts;

    for (size_t i = 0; i < positions.size(); i += 3){
        glm::vec3 v0 = positions[i];
        glm::vec3 v1 = positions[i+1];
        glm::vec3 v2 = positions[i+2];

        glm::vec3 edge1 = v1 - v0;
        glm::vec3 edge2 = v2 - v0;
        glm::vec3 face_normal = glm::normalize(glm::cross(edge1,edge2));

        for (const glm::vec3& v : {v0, v1 ,v2}){
            glm::vec3 key = v;
            normalSums[key] += face_normal;
            normalCounts[key] += 1;
        }
    }

    std::vector<float> normalData;
    normalData.reserve(verticies.size());

    for (const glm::vec3& v : positions){
        glm::vec3 key = v;
        glm::vec3 normal = glm::normalize(normalSums[key] / static_cast<float>(normalCounts[key]));
        normalData.push_back(normal.x);
        normalData.push_back(normal.y);
        normalData.push_back(normal.z);
    }

    return normalData;
}

#endif 
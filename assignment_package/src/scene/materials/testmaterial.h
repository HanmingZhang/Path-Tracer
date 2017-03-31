#pragma once
#include "material.h"

class TestMaterial : public Material
{
public:
    TestMaterial(const Color3f &Kt, float roughness, float indexOfRefraction,
                   const std::shared_ptr<QImage> &roughnessMap,
                   const std::shared_ptr<QImage> &textureMap,
                   const std::shared_ptr<QImage> &normalMap)
        : Kt(Kt), roughness(roughness), indexOfRefraction(indexOfRefraction), roughnessMap(roughnessMap),
          textureMap(textureMap), normalMap(normalMap)
    {}

    void ProduceBSDF(Intersection *isect) const;


private:
    Color3f Kt;                    // The spectral reflection of this material.
                                   // This is just the base color of the material

    float roughness;               // The overall roughness of the material. Will be multiplied
                                   // with the roughness map's values.

    float indexOfRefraction;


    std::shared_ptr<QImage> roughnessMap;

    std::shared_ptr<QImage> textureMap; // The color obtained from this (assuming it is non-null) is multiplied with the base material color.
    std::shared_ptr<QImage> normalMap;
};

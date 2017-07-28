#pragma
#include "material.h"



// MetalMaterial Declarations
class MetalMaterial : public Material {

public:
    // MetalMaterial Public Methods
    MetalMaterial(const Color3f &eta,
                  const Color3f &k,
                  const float &rough,
                  const float &urough,
                  const float &vrough,
                  bool remapRoughness
                  )
        : eta(eta), k(k), roughness(rough), uRoughness(urough), vRoughness(vrough), remapRoughness(remapRoughness)
    {}

    void ProduceBSDF(Intersection *isect) const;

private:
    // MetalMaterial Private Data
    Color3f eta, k;
    float roughness, uRoughness, vRoughness;
    bool remapRoughness;
};




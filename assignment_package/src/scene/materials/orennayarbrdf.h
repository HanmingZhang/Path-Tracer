#pragma once
#include "bsdf.h"

class OrenNayarBRDF : public BxDF
{
public:
    OrenNayarBRDF(const Color3f &R, float sigma)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), R(R) {
        float sigma2 = glm::radians(sigma) * glm::radians(sigma);
        A = 1.0f - (sigma2 / (2.0f * (sigma2 + 0.33f)));
        B = 0.45f * sigma2 / (sigma2 + 0.09f);
    }

    Color3f f(const Vector3f &wo, const Vector3f &wi) const;



  private:
    const Color3f R; // The energy scattering coefficient of this BRDF (i.e. its color)
    float A, B;
};



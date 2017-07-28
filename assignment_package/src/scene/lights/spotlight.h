#pragma once
#include "light.h"


class SpotLight : public Light
{
public:
    SpotLight(const Transform &t, const Color3f& Le, float totalWidth, float falloffStart) :
        Light(t),emittedLight(Le),pLight(glm::vec3(t.T() * glm::vec4(0.f, 0.f, 0.f, 1.f))),
        cosTotalWidth(std::cos(glm::radians(totalWidth))),cosFalloffStart(std::cos(glm::radians(falloffStart)))
    {}

    virtual Color3f L(const Intersection &isect, const Vector3f &w) const;

    virtual Color3f Sample_Li(const Intersection &ref, const Point2f &xi,
                                                Vector3f *wi, float *pdf) const;

    virtual float Pdf_Li(const Intersection &ref, const Vector3f &wi) const;

    float Falloff(const Vector3f &w) const;

    virtual Ray generatePhotonRay(const Point2f &sample, Sampler *sampler, float *pdf);

    virtual Color3f GetInitialPhotonPower() const;


    const Color3f emittedLight;
    const Point3f pLight;
    const float cosTotalWidth, cosFalloffStart;
};

#pragma once
#include "integrator.h"
#include "directlightingintegrator.h"

class FullLightingIntegrator : public Integrator
{
public:
    FullLightingIntegrator(Bounds2i bounds, Scene* s, std::shared_ptr<Sampler> sampler, int recursionLimit)
        : Integrator(bounds, s, sampler, recursionLimit)
    {}

    // Evaluate the energy transmitted along the ray back to
    // its origin using multiple importance sampling
    virtual Color3f Li(const Ray &ray, const Scene &scene, std::shared_ptr<Sampler> sampler, int depth) const;


};

float BalanceHeuristic(int nf, float fPdf, int ng, float gPdf);
float PowerHeuristic(int nf, float fPdf, int ng, float gPdf);


Color3f EstimateDirect(const Intersection &it,
                       const Point2f &uScattering, const Light &light,
                       const Point2f &uLight, const Scene &scene, std::shared_ptr<Sampler> sampler,
                       const Vector3f &woW,
                       bool specular = false);

Color3f UniformSampleOneLight(const Intersection &it, const Scene &scene, std::shared_ptr<Sampler> sampler, const Vector3f &woW, bool specularBounce);


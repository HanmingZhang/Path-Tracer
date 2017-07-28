#pragma once
#include "integrator.h"
#include "directlightingintegrator.h"


class PhotonMappingIntegrator : public Integrator
{
public:
    PhotonMappingIntegrator(Bounds2i bounds, Scene* s, std::shared_ptr<Sampler> sampler, int recursionLimit)
        : Integrator(bounds, s, sampler, recursionLimit)
    {}


    virtual Color3f Li(const Ray &ray, const Scene &scene, std::shared_ptr<Sampler> sampler, int depth) const;

};



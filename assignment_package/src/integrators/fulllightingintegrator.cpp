#include "fulllightingintegrator.h"
#include <iostream>

// After how many times do we trigger Russian Roulette
static const int RussianRouletteTriggerTime = 3;


Color3f FullLightingIntegrator::Li(const Ray &ray, const Scene &scene, std::shared_ptr<Sampler> sampler, int depth) const
{
    //TODO
    Color3f AccumulatedColor(0.f);
    Color3f ThroughputColor(1.f);

    Ray traceRay = ray;
    //Intersection isect;
    bool specularBounce = false;

    Vector3f woW;


    while(depth > 0) {

        woW = -traceRay.direction;

        // Intersect ray with scene
        Intersection isect;
        bool foundIntersection = scene.Intersect(traceRay, &isect);


        // Terminate path if ray escaped or max depth is reached
        if(!foundIntersection){
            break;
        }


        // Possibly add emitted light at intersection
        if (depth == recursionLimit || specularBounce){
            // Add emitted light at path vertex or from the environment
            if(foundIntersection){
                AccumulatedColor += ThroughputColor * isect.Le(woW);
            }
            else{
                for(const auto &light : scene.lights){
                    // radiance emitted by infinite area lgiht sources, if present
                    // direction light ?
                    AccumulatedColor += ThroughputColor * light->Le(traceRay);
                }
            }
        }


        //If this ray is not acamera ray and
        //it hits a light(has not material and thus, no bsdf)
//        if(isect.objectHit->GetAreaLight()){
//             break;
//        }
        if(isect.objectHit->GetLight()){
             break;
        }

        // Compute scattering functions
        isect.ProduceBSDF();



        // check whether this ray will hit a specular surface
//        if(isect.bsdf->BxDFsHaveSpecularFlags()){
//            specularBounce = true;
//        }
//        else{
//            specularBounce = false;
//        }



        //Sample illumination from lights to find path contribution
        //if it's a specular bounce ray, we ignore the direct lighting
        AccumulatedColor += ThroughputColor * UniformSampleOneLight(isect, scene, sampler, woW, specularBounce);


        // Sample BSDF to get new path direction
        Vector3f wi(0.f);
        float pdf = 0.f;

        BxDFType flags;

        Color3f f = isect.bsdf->Sample_f(woW, &wi, sampler->Get2D(), &pdf, BSDF_ALL, &flags);

        if(IsBlack(f) || pdf == 0.f){
            break;
        }

        ThroughputColor *= (f * AbsDot(wi, glm::normalize(isect.normalGeometric)) / pdf);

        specularBounce = ((flags & BSDF_SPECULAR) != 0);


        // update the ray  from current intersection
        // and travel in the direction of the ωi just computed.
        traceRay = isect.SpawnRay(wi);


        // Possibly terminate the path with Russian roulette
        if(depth <= recursionLimit - RussianRouletteTriggerTime){
            float max_component = std::max(ThroughputColor.x, std::max(ThroughputColor.y, ThroughputColor.z));
            if(max_component < sampler->Get1D()){
                break;
            }
            ThroughputColor /= max_component;
        }

        depth--;
    }

    return AccumulatedColor;
}


float BalanceHeuristic(int nf, float fPdf, int ng, float gPdf)
{
    //TODO
    return ((float)nf * fPdf + (float)ng * gPdf == 0.f) ? 0.f : (((float)nf * fPdf) / ((float)nf * fPdf + (float)ng * gPdf));
}

float PowerHeuristic(int nf, float fPdf, int ng, float gPdf)
{
    //TODO
    float f = (float)nf * fPdf;
    float g = (float)ng * gPdf;
    return ((f * f + g * g) == 0.f) ? 0.f : ((f * f) / (f * f + g * g));
}


Color3f EstimateDirect(const Intersection &it,
                       const Point2f &uScattering, const Light &light,
                       const Point2f &uLight, const Scene &scene, std::shared_ptr<Sampler> sampler, const Vector3f &woW,
                       bool specular){

//    BxDFType bsdfFlags = specular ? BSDF_ALL : BxDFType(BSDF_ALL & ~BSDF_SPECULAR);

    Color3f Ld(0.f);


    // Sample light source with multiple importance sampling
    Vector3f wi_Direct(0.f);


    float lightPdf = 0.f;
    float scatteringPdf = 0.f;

    Color3f Li = light.Sample_Li(it, uLight, &wi_Direct, &lightPdf);


    // shadow test
    Ray shadowFeelerRay = it.SpawnRay(glm::normalize(wi_Direct));

    if(!shadowTest(shadowFeelerRay, scene, glm::length(wi_Direct), &light) && lightPdf > 0.f && !IsBlack(Li)){

        // we normalize wi here!
        // before we do normalize, we must check wiW is not (0, 0, 0)
        if(wi_Direct.x != 0.f ||
           wi_Direct.y != 0.f ||
           wi_Direct.z != 0.f){
            wi_Direct = glm::normalize(wi_Direct);
        }


        // Compute BSDF value for lgiht sample
        Color3f f(0.f);

        // Evaluate BSDF for light sampling strategy
        const Intersection &isect = (const Intersection &)it;
//        f = isect.bsdf->f(woW, wi_Direct, bsdfFlags) * AbsDot(wi_Direct, glm::normalize(isect.normalGeometric));
//        scatteringPdf = isect.bsdf->Pdf(woW, wi_Direct, bsdfFlags);

        f = isect.bsdf->f(woW, wi_Direct) * AbsDot(wi_Direct, glm::normalize(isect.normalGeometric));
        scatteringPdf = isect.bsdf->Pdf(woW, wi_Direct);




        if(!IsBlack(f)){
            // Add light's contribution to reflected radiance
            if(!IsBlack(Li)){
                float weight = PowerHeuristic(1, lightPdf, 1, scatteringPdf);
                Ld += f * Li * weight / lightPdf;
            }
        }

    }


    if(specular){
        return Ld;
    }


    // Sample BSDF with multiple importance sampling
    Color3f f(0.f);
    Vector3f wi_BSDF(0.f);

    lightPdf = 0.f;
    scatteringPdf = 0.f;

    bool sampledSpecular = false;

    // Sample scattered direction for surface interactions
    BxDFType sampledType;
    const Intersection &isect = (const Intersection &)it;

    //f = isect.bsdf->Sample_f(woW, &wi_BSDF, uScattering, &scatteringPdf, bsdfFlags, &sampledType);
    f = isect.bsdf->Sample_f(woW, &wi_BSDF, uScattering, &scatteringPdf, BSDF_ALL, &sampledType);
    f *= AbsDot(wi_BSDF, glm::normalize(isect.normalGeometric));

    sampledSpecular = sampledType & BSDF_SPECULAR;


    if(!IsBlack(f) && scatteringPdf > 0.f){
        // Account for light contributions along sampled direction wi
        float weight = 1.f;

        if(!sampledSpecular){
            lightPdf = light.Pdf_Li(it, wi_BSDF);

            if(lightPdf == 0.f){
                return Ld;
            }

            weight = PowerHeuristic(1, scatteringPdf, 1, lightPdf);
        }

        // Find intersection and compute transmittance
        Intersection lightIsect;


        Ray ray = it.SpawnRay(wi_BSDF);
//        Color3f Tr(1.f);
        bool foundIntersection = scene.Intersect(ray, &lightIsect);

        // Add light contribution form material sampling
        Color3f Li(0.f);

        if(foundIntersection){
            //if(lightIsect.objectHit->GetAreaLight() == &light){
            if(lightIsect.objectHit->GetLight() == &light){
                Li = lightIsect.Le(-wi_BSDF);
            }
        }
        else{
            Li = light.Le(ray);
        }
        if(!IsBlack(Li)){
            Ld += f * Li * weight / scatteringPdf;
        }


    }

    return Ld;
}

Color3f UniformSampleOneLight(const Intersection &it, const Scene &scene, std::shared_ptr<Sampler> sampler, const Vector3f &woW, bool specularBounce){

    // Randomly choose a single light to sample
    int nLights = scene.lights.size();
    if (nLights == 0) return Color3f(0.f);

    int lightNum = std::min((int)(sampler->Get1D() * nLights), nLights - 1);

    const std::shared_ptr<Light> &light = scene.lights[lightNum];


    Point2f uLight = sampler->Get2D();
    Point2f uScattering = sampler->Get2D();


    return (float)nLights * EstimateDirect(it, uScattering, *light, uLight, scene, sampler, woW, specularBounce);
}



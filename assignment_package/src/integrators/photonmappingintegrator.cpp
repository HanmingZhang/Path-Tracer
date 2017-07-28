#include "photonmappingintegrator.h"

Color3f PhotonMappingIntegrator::Li(const Ray &ray, const Scene &scene, std::shared_ptr<Sampler> sampler, int depth) const
{
    Color3f L(0.f);
    Color3f L_DirectLighting(0.f);
    Color3f L_CausticPhotonMap(0.f);
    Color3f L_GlobalPhotonMap(0.f);
    Color3f L_Specular(0.f);

    Vector3f woW = -ray.direction;


    // First, we get the intersection of this ray
    Intersection isect;
    scene.Intersect(ray, &isect);

    // if there is no intersection at all
    if(isect.objectHit == nullptr){
        //actually, we return background color
        return Color3f(0.0f);
    }


    // if intersection exists, add Le first
    L += isect.Le(woW);


    // if the intersection with the scene hits an object with no Material attached,
    // then L should only evaluate and return the light emitted directly from the intersection
    if(depth == 0 || isect.objectHit->material == nullptr){
        return L;
    }

    // set up BSDF
    isect.ProduceBSDF();



    // -------------------- Part I : Direct Lighting -------------------------
    Vector3f wiW(0.0f);
    float pdf = 0.0f;


    // we randomly choose a light source
    float lightChooseRandomNum = sampler->Get1D();
    int lightNum = scene.lights.size();
    int lightIdx = std::min((int)std::floor(lightChooseRandomNum * (float)lightNum),
                            lightNum - 1);



    // wiW will change here !
    Color3f Li_DirectLighting = scene.lights[lightIdx]->Sample_Li(isect, sampler->Get2D(), &wiW, &pdf);

    // after we get wiW, we do shadow test here
    // at this time, wiW is not normalized!
    // this shadow test method is based on distance

    Ray shadowFeelerRay = isect.SpawnRay(glm::normalize(wiW));
    if(!shadowTest(shadowFeelerRay, scene, glm::length(wiW), scene.lights[lightIdx].get())){

        // we normalize wiW here!
        // before we do normalize, we must check wiW is not (0, 0, 0)
        if(wiW.x != 0.f ||
           wiW.y != 0.f ||
           wiW.z != 0.f){
            wiW = glm::normalize(wiW);
        }

        Color3f f_term_DirectLighting = isect.bsdf->f(woW, wiW);

        // check pdf
        if(pdf != 0.f){
            // divide the PDF obtain from Sample_Li by the number of light sources in your scene
            pdf = pdf / (float)lightNum;

            // directly use LTE here!
            L_DirectLighting= (f_term_DirectLighting * Li_DirectLighting * AbsDot(wiW,glm::normalize(isect.normalGeometric))) / pdf;
        }
    }



    // -------------------- Part II : Global Photon -------------------------
    if(scene.globalPhotonMap.NumPhotons() > 10){
        int globalSamples = 0;
        int numOfProbes = 15;

        for(int i = 0; i < numOfProbes; i++){
            wiW = Vector3f(0.0f);
            pdf = 0.0f;
            Color3f Li_GlobalPhotonMap(0.f);


            // wiW will change here !
            Color3f f_term_GlobalPhotonMap = isect.bsdf->Sample_f(woW, &wiW, sampler->Get2D(), &pdf);

            if(pdf != 0.f){
                // we normalize wiW here!
                // before we do normalize, we must check wiW is not (0, 0, 0)
                if(wiW.x != 0.f ||
                   wiW.y != 0.f ||
                   wiW.z != 0.f){
                    wiW = glm::normalize(wiW);
                }
                Ray GlobalPhotonRay = isect.SpawnRay(wiW);
                Intersection GlobalPhotonIsect;
                scene.Intersect(GlobalPhotonRay, &GlobalPhotonIsect);

                if(GlobalPhotonIsect.objectHit != nullptr){
                    Vector3f dir(0.f);
                    scene.globalPhotonMap.EstimateIrradiance(Li_GlobalPhotonMap, dir, 0.5, GlobalPhotonIsect.point, 150, &GlobalPhotonIsect.normalGeometric);
                    L_GlobalPhotonMap += (f_term_GlobalPhotonMap * Li_GlobalPhotonMap * AbsDot(wiW,glm::normalize(isect.normalGeometric))) / pdf;
                    globalSamples++;
                }
            }
        }

        if(globalSamples != 0) {
            L_GlobalPhotonMap = L_GlobalPhotonMap / (float)globalSamples;
        }
    }



    // -------------------- Part III : Caustic Photon -------------------------
    if(scene.causticPhotonMap.NumPhotons() > 10){

        wiW = Vector3f(0.0f);
        pdf = 0.0f;
        Color3f Li_CausticPhotonMap(0.f);

        scene.causticPhotonMap.EstimateIrradiance(Li_CausticPhotonMap, wiW, 0.2, isect.point, 150, &isect.normalGeometric);

        wiW = -wiW;

        Color3f f_term_CausticPhotonMap = isect.bsdf->f(woW, wiW);

        pdf = isect.bsdf->Pdf(woW, wiW);

        L_CausticPhotonMap = (pdf == 0.f) ? Color3f(0.f) : f_term_CausticPhotonMap * Li_CausticPhotonMap * AbsDot(wiW,glm::normalize(isect.normalGeometric)) / pdf;
    }



    // -------------------- Part IV : Specular Illumination -------------------------
    if(isect.bsdf->BxDFsHaveSpecularFlags()){
        wiW = Vector3f(0.0f);
        pdf = 0.0f;

        Color3f f_term_Specular = isect.bsdf->Sample_f(woW, &wiW, sampler->Get2D(), &pdf);

        Ray SpecularRay = isect.SpawnRay(wiW);

        L_Specular = (pdf == 0.f) ? Color3f(0.f) : f_term_Specular * Li(SpecularRay, scene, sampler, depth - 1) * AbsDot(wiW,glm::normalize(isect.normalGeometric)) / pdf;
    }




    //L = L_GlobalPhotonMap;
    //L = L_CausticPhotonMap;
    L += L_DirectLighting + L_CausticPhotonMap + L_GlobalPhotonMap + L_Specular;

    if(std::isnan(L.x) ||
       std::isnan(L.y) ||
       std::isnan(L.z)){
        std::cout << "stop here, L nan!" << std::endl;
    }

    return L;
}

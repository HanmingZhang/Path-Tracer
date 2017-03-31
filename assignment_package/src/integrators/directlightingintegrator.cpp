#include "directlightingintegrator.h"


bool shadowTest(const Ray &ray, const Scene &scene, float distanceFromIsectToLight){
    Intersection test_isect;
    scene.Intersect(ray, &test_isect);

    // Since this is direct lighting, we don't consider the transmissive condition here
    // we just compare distances

    float test_distance = glm::length(test_isect.point - ray.origin);

    if(test_isect.t > 0.f && fabs(test_distance - distanceFromIsectToLight) < ShadowEpsilon){

        return false;
    }

    else return true;

//    if(test_isect.objectHit != nullptr){
//        if(test_isect.objectHit->GetAreaLight() == &light){
//            return false;
//        }
//    }
//    else return true;
}


Color3f DirectLightingIntegrator::Li(const Ray &ray, const Scene &scene, std::shared_ptr<Sampler> sampler, int depth) const
{
    //TODO
    Color3f L(0.0f);
    Vector3f woW = -ray.direction;


    // First, we get the intersection of this ray
    Intersection isect;
    scene.Intersect(ray, &isect);


    // if there is no intersection at all
    if(isect.objectHit == nullptr){
        //actually, we return background color
        return Color3f(0.0f);
    }


    // if intersection exists, compute Le first
    L += isect.Le(woW);


    // if the intersection with the scene hits an object with no Material attached,
    // then Li should only evaluate and return the light emitted directly from the intersection
    if(isect.objectHit->material == nullptr){
        return L;
    }


    // set up BSDF
    isect.ProduceBSDF();



    Vector3f wiW(0.0f);
    float pdf = 0.0f;


    // we randomly choose a light source
    float lightChooseRandomNum = sampler->Get1D();
    int lightNum = scene.lights.size();
    int lightIdx = std::min((int)std::floor(lightChooseRandomNum * (float)lightNum),
                            lightNum - 1);



    // wiW will change here !
    Color3f Li = scene.lights[lightIdx]->Sample_Li(isect, sampler->Get2D(), &wiW, &pdf);

    // after we get wiW, we do shadow test here
    // at this time, wiW is not normalized!
    // this shadow test method is based on distance

    Ray shadowFeelerRay = isect.SpawnRay(glm::normalize(wiW));
    if(shadowTest(shadowFeelerRay, scene, glm::length(wiW))){
        return Color3f(0.0f);
    }




    // we normalize wiW here!
    // before we do normalize, we must check wiW is not (0, 0, 0)
    if(wiW.x != 0.f ||
       wiW.y != 0.f ||
       wiW.z != 0.f){
        wiW = glm::normalize(wiW);
    }




    Color3f f_term = isect.bsdf->f(woW, wiW);


    // check pdf
    if(pdf == 0.f){
        return L;
    }


    // divide the PDF obtain from Sample_Li by the number of light sources in your scene
    pdf = pdf / (float)lightNum;



    // directly use LTE here!
    L += (f_term * Li * AbsDot(wiW,glm::normalize(isect.normalGeometric))) / pdf;


    return L;
}

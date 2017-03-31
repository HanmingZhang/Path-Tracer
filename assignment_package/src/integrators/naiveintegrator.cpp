#include "naiveintegrator.h"

Color3f NaiveIntegrator::Li(const Ray &ray, const Scene &scene, std::shared_ptr<Sampler> sampler, int depth) const
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


    // if intersection exists
    L += isect.Le(woW);


    // Note that if Li is invoked with a depth value of 0
    // or if the intersection with the scene hits an object with no Material attached,
    // then Li should only evaluate and return the light emitted directly from the intersection
    if(depth == 0 || isect.objectHit->material == nullptr){
        return L;
    }



    isect.ProduceBSDF();



    Vector3f wiW(0.0f);
    float pdf = 0.0f;

    // wiW will change here !
    Color3f f_term = isect.bsdf->Sample_f(woW, &wiW, sampler->Get2D(), &pdf);

    // check pdf
    if(pdf == 0.f){
        return L;
    }

    // we normalize wiW here!
    // before we do normalize, we must check wiW is not (0, 0, 0)
    if(wiW.x != 0.f ||
       wiW.y != 0.f ||
       wiW.z != 0.f){
        wiW = glm::normalize(wiW);
    }

    // we generate a new ray based on intersection and wiW
    Ray newRay = isect.SpawnRay(wiW);


    // recursion happens here !
    L += (f_term * Li(newRay, scene, sampler, depth - 1) * AbsDot(wiW,glm::normalize(isect.normalGeometric))) / pdf;



    return L;
}

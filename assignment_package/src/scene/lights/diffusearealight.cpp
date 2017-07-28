#include "diffusearealight.h"

Color3f DiffuseAreaLight::L(const Intersection &isect, const Vector3f &w) const
{
    //TODO

    if(twoSided){
        return emittedLight;
    }
    else return glm::dot(isect.normalGeometric, w) > 0.f ? emittedLight : Color3f(0.0f);
}



Color3f DiffuseAreaLight::Sample_Li(const Intersection &ref, const Point2f &xi,
                                     Vector3f *wi, float *pdf) const{

    *pdf = 0.f;

    Intersection pShape = shape->Sample(ref, xi, pdf);

    // for the sake of distance-based shadow test method,
    // we temporarily don't normalize wi here
    *wi = pShape.point - ref.point;
    //*wi = glm::normalize(pShape.point - ref.point);



    // Check if the resultant PDF is zero
    // or that the reference Intersection and the resultant Intersection are the same point in space,
    // and return black if this is the case

    if(*pdf == 0.f ||
       (fabs(pShape.point.x - ref.point.x) < 0.1f &&
        fabs(pShape.point.y - ref.point.y) < 0.1f &&
        fabs(pShape.point.z - ref.point.z) < 0.1f)){
        return Color3f(0.f);
    }



    return L(pShape, -1.0f * (*wi));

}

float DiffuseAreaLight::Pdf_Li(const Intersection &ref, const Vector3f &wi) const{

    return shape->Pdf(ref,wi);

}

Ray DiffuseAreaLight::generatePhotonRay(const Point2f &sample, Sampler *sampler, float *pdf){

    //float pdf = 0.f;
    //Point2f xi = sampler->Get2D();

    Intersection Origin = shape->Sample(sample, pdf);

    Vector3f Dir = Vector3f(0.f);

    if(twoSided){
        Dir = WarpFunctions::squareToSphereUniform(sample);
        //Dir = WarpFunctions::squareToHemisphereCosine()
    }
    else{
        Dir = WarpFunctions::squareToHemisphereCosine(sample);
    }

    Dir = glm::normalize(Dir);

    Dir = shape->transform.RotateT() * Dir;

    Dir = glm::normalize(Dir);

    return Ray(Origin.point + RayEpsilon * Origin.normalGeometric, Dir);

}


Color3f DiffuseAreaLight::GetInitialPhotonPower() const{
    return emittedLight;
}

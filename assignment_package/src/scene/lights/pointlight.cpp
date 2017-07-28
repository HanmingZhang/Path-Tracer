#include "pointlight.h"

Color3f PointLight::L(const Intersection &isect, const Vector3f &w) const{
    return emittedLight;
}


Color3f PointLight::Sample_Li(const Intersection &ref, const Point2f &xi,
                                            Vector3f *wi, float *pdf) const{

    //*wi = glm::normalize(pLight - ref.point);
    *wi = pLight - ref.point;
    *pdf = 1.f;
    return emittedLight / glm::length2(pLight - ref.point);
    //return emittedLight;
}



float PointLight::Pdf_Li(const Intersection &ref, const Vector3f &wi) const{

    return 0.f;
}


Ray PointLight::generatePhotonRay(const Point2f &sample, Sampler *sampler, float *pdf){

    float x = sampler->Get1D();
    float y = sampler->Get1D();
    float z = sampler->Get1D();

    (*pdf) = 1.0f;

    Vector3f Dir(x, y, z);
    Dir = glm::normalize(Dir);

    return Ray(pLight, Dir);
}

Color3f PointLight::GetInitialPhotonPower() const{
    return emittedLight;
}

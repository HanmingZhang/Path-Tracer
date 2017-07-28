#include "spotlight.h"

Color3f SpotLight::L(const Intersection &isect, const Vector3f &w) const{
    return emittedLight;
}


Color3f SpotLight::Sample_Li(const Intersection &ref, const Point2f &xi,
                                            Vector3f *wi, float *pdf) const{

    //*wi = glm::normalize(pLight - ref.point);
    *wi = pLight - ref.point;
    *pdf = 1.f;
    return emittedLight * Falloff(-(glm::normalize(*wi)))/ glm::length2(pLight - ref.point);
}



float SpotLight::Pdf_Li(const Intersection &ref, const Vector3f &wi) const{

    return 0.f;
}

float SpotLight::Falloff(const Vector3f &w) const{
    //Vector3f wl = glm::normalize(glm::vec3(transform.invT() * glm::vec4(w, 0.f)));

    Vector3f wl = glm::normalize(transform.RotateT() * w);

    //float cosTheta = wl.z;

    //It's a little bit tricky,
    //but we use -y value here!
    float cosTheta = -wl.y;


    if(cosTheta < cosTotalWidth) return 0.f;
    if(cosTheta > cosFalloffStart) return 1.0f;

    // Compute falloff inside spotlight cone
    float delta = (cosTheta - cosTotalWidth) / (cosFalloffStart - cosTotalWidth);

    return delta;
    //return delta * delta;
    //return delta * delta * delta * delta;
}


Ray SpotLight::generatePhotonRay(const Point2f &sample, Sampler *sampler, float *pdf){

    Vector3f Dir(0.f);

    float choice = sampler->Get1D();

    if(choice < 0.7f){
        do{
            float x = sampler->Get1D();
            float y = sampler->Get1D();
            float z = sampler->Get1D();

            Dir = Vector3f(x, y, z);
            Dir = glm::normalize(Dir);

        }while(-Dir.y >= cosFalloffStart);
    }
    else{
        do{
            float x = sampler->Get1D();
            float y = sampler->Get1D();
            float z = sampler->Get1D();

            Dir = Vector3f(x, y, z);
            Dir = glm::normalize(Dir);

        }while((-Dir.y >= cosTotalWidth) && (-Dir.y > cosFalloffStart));
    }

    (*pdf) = 1.0f;

    return Ray(pLight, Dir);
}

Color3f SpotLight::GetInitialPhotonPower() const{
    return emittedLight;
}

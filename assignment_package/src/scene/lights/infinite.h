#pragma once
#include "light.h"
#include "distribution2d.h"
//#include"mipmap.h"

// InfiniteAreaLight Declarations
class InfiniteAreaLight : public Light {
  public:
    // InfiniteAreaLight Public Methods
    InfiniteAreaLight(const Transform &t, const Color3f &power,
                      const std::shared_ptr<QImage> &texmap);



    Color3f Power() const;

    virtual Color3f Le(const Ray &ray) const;

    virtual Color3f L(const Intersection &isect, const Vector3f &woW) const;

    virtual Color3f Sample_Li(const Intersection &ref, const Point2f &u, Vector3f *wi,
                       Float *pdf) const;
    virtual float Pdf_Li(const Intersection &, const Vector3f &) const;

//    Color3f Sample_Le(const Point2f &u1, const Point2f &u2, Float time,
//                       Ray *ray, Normal3f *nLight, Float *pdfPos,
//                       Float *pdfDir) const;
//    void Pdf_Le(const Ray &, const Normal3f &, Float *pdfPos,
//                Float *pdfDir) const;

    virtual Ray generatePhotonRay(const Point2f &u, Sampler *sampler, float *pdf);

    virtual Color3f GetInitialPhotonPower() const;

  private:
    // InfiniteAreaLight Private Data
    std::shared_ptr<QImage> texmap;

    Color3f lightColor;
    Point3f worldCenter;
    Float worldRadius;
    std::unique_ptr<Distribution2D> distribution;

    Color3f Lookup(Point2f &uv_coord) const;

    Color3f tempPower;
};

//std::shared_ptr<InfiniteAreaLight> CreateInfiniteLight(
//const Transform &light2world, const ParamSet &paramSet);


inline Float SphericalTheta(const Vector3f &v){
    //return std::acos(glm::clamp(v.z, -1.f, 1.f));
    return std::acos(glm::clamp(v.y, -1.f, 1.f));
}

inline Float SphericalPhi(const Vector3f &v){
    float p = std::atan2(v.z, v.x);
    return (p < 0.f) ? (p + 2.f * Pi) : p;
}

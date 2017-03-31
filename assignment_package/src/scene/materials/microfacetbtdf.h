#pragma once
#include "bsdf.h"
#include "fresnel.h"
#include "microfacet.h"


class MicrofacetBTDF : public BxDF {
  public:
    // MicrofacetTransmission Public Methods
    MicrofacetBTDF(const Color3f &T,
                           MicrofacetDistribution *distribution,
                           float etaA,float etaB)
        : BxDF(BxDFType(BSDF_TRANSMISSION | BSDF_GLOSSY)),
          T(T),
          distribution(distribution),
          etaA(etaA),
          etaB(etaB),
          fresnel(etaA, etaB) {}

    virtual ~MicrofacetBTDF(){delete distribution;}


    Color3f f(const Vector3f &wo, const Vector3f &wi) const;
    virtual Color3f Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &xi,
                      float *pdf, BxDFType *sampledType = nullptr) const;
    virtual float Pdf(const Vector3f &wo, const Vector3f &wi) const;
    // std::string ToString() const;

  private:
    // MicrofacetTransmission Private Data
    const Color3f T;
    const MicrofacetDistribution *distribution;
    const float etaA, etaB;
    const FresnelDielectric fresnel;
    // const TransportMode mode;
};

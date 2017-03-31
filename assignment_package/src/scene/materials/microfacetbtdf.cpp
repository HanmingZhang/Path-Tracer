#include "microfacetbtdf.h"

Color3f MicrofacetBTDF::f(const Vector3f &wo,const Vector3f &wi) const {
    if (SameHemisphere(wo, wi)) return Color3f(0.f);  // transmission only

    float cosThetaO = CosTheta(wo);
    float cosThetaI = CosTheta(wi);
    if (cosThetaI == 0.f || cosThetaO == 0.f) return Color3f(0.f);

    // Compute wh from wo and wi for microfacet transmission
    float eta = CosTheta(wo) > 0 ? (etaB / etaA) : (etaA / etaB);
    Vector3f wh = glm::normalize(wo + wi * eta);
    if (wh.z < 0) wh = -wh;

    Color3f F = fresnel.Evaluate(glm::dot(wo, wh));

    float sqrtDenom = glm::dot(wo, wh) + eta * glm::dot(wi, wh);
    // float factor = (mode == TransportMode::Radiance) ? (1 / eta) : 1;
    float factor = 1.f / eta;

    return (Color3f(1.f) - F) * T *
           std::abs(distribution->D(wh) * distribution->G(wo, wi) * eta * eta *
                    AbsDot(wi, wh) * AbsDot(wo, wh) * factor * factor /
                    (cosThetaI * cosThetaO * sqrtDenom * sqrtDenom));
}


Color3f MicrofacetBTDF::Sample_f(const Vector3f &wo, Vector3f *wi,
                                 const Point2f &xi, float *pdf,
                                 BxDFType *sampledType) const {
    if (wo.z == 0) return Color3f(0.f);
    Vector3f wh = distribution->Sample_wh(wo, xi);
    float eta = CosTheta(wo) > 0 ? (etaA / etaB) : (etaB / etaA);

//    bool entering = CosTheta(wo) > 0;
//    float etaI = entering ? etaA : etaB;
//    float etaT = entering ? etaB : etaA;


  if (!Refract(wo, Faceforward((Normal3f)wh, wo), eta, wi)) return Color3f(0.f);

  // if (!Refract(wo, Faceforward(Normal3f(0.f, 0.f, 1.f), wo), eta, wi)) return Color3f(0.f);
    *pdf = Pdf(wo, *wi);
    return f(wo, *wi);
}



float MicrofacetBTDF::Pdf(const Vector3f &wo,const Vector3f &wi) const {
    if (SameHemisphere(wo, wi)) return 0;
    // Compute wh from wo and wi for microfacet transmission
    float eta = CosTheta(wo) > 0 ? (etaB / etaA) : (etaA / etaB);
    Vector3f wh = glm::normalize(wo + wi * eta);

    // Compute change of variables _dwh\_dwi_ for microfacet transmission
    float sqrtDenom = glm::dot(wo, wh) + eta * glm::dot(wi, wh);
    float dwh_dwi = std::abs((eta * eta * glm::dot(wi, wh)) / (sqrtDenom * sqrtDenom));
    return distribution->Pdf(wo, wh) * dwh_dwi;
}

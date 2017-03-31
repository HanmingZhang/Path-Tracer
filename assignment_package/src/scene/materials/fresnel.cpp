#include "fresnel.h"

Color3f FresnelDielectric::Evaluate(float cosThetaI) const
{
    //TODO
    return Color3f(FrDielectric(cosThetaI, etaI, etaT));
}

Color3f FresnelConductor::Evaluate(float cosThetaI) const
{
    return FrConductor(std::fabs(cosThetaI), etaI, etaT, k);
}

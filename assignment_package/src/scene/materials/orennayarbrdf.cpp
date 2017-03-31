#include "orennayarbrdf.h"


Color3f OrenNayarBRDF::f(const Vector3f &wo, const Vector3f &wi) const
{
    //TODO

    float sinThetaI = SinTheta(wi);
    float sinThetaO = SinTheta(wo);

    // Compute cosine term of Oren-Nayar model
    float maxCos = 0;
    if(sinThetaI > 1e-4 && sinThetaO > 1e-4){
        float sinPhiI = SinPhi(wi), cosPhiI = CosPhi(wi);
        float sinPhiO = SinPhi(wo), cosPhiO = CosPhi(wo);

        float dCos = cosPhiI * cosPhiO + sinPhiI * sinPhiO;
        maxCos = std::max(0.f, dCos);
    }


    // Compute sine and tangent terms of Oren-Nayar model
    float sinAlpha, tanBeta;
    if(AbsCosTheta(wi) > AbsCosTheta(wo)){
        sinAlpha = sinThetaO;
        tanBeta = sinThetaI / AbsCosTheta(wi);
    }
    else{
        sinAlpha = sinThetaI;
        tanBeta = sinThetaO / AbsCosTheta(wo);
    }



    return R * InvPi * (A + B * maxCos * sinAlpha * tanBeta);
}

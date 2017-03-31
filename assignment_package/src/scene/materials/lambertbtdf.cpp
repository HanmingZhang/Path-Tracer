#include "lambertbtdf.h"
#include <warpfunctions.h>


Color3f LambertBTDF::f(const Vector3f &wo, const Vector3f &wi) const
{
    //TODO
    return R * InvPi;
}

Color3f LambertBTDF::Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                        Float *pdf, BxDFType *sampledType) const
{
    //TODO


    *wi = WarpFunctions::squareToHemisphereCosine(u);

    if(wo.z > 0) wi->z *= -1;

    *pdf = Pdf(wo, *wi);

    return f(wo, *wi);



    // return BxDF::Sample_f(wo, wi, u, pdf, sampledType);
}

float LambertBTDF::Pdf(const Vector3f &wo, const Vector3f &wi) const
{
    //TODO

    return SameHemisphere(wo, wi) ? wi[2] * InvPi : 0;

    // return BxDF::Pdf(wo, wi);
}

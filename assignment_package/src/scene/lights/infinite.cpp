#include "infinite.h"

// InfiniteAreaLight Method Definitions
InfiniteAreaLight::InfiniteAreaLight(const Transform &t,
                                     const Color3f &L,
                                     const std::shared_ptr<QImage> &texmap)
    : Light(t),texmap(texmap),lightColor(L),tempPower(Color3f(0.f))
{
    // Read texel data from _texmap_ and initialize _Lmap_
    Point2i resolution(texmap->width(), texmap->height());



    // Initialize sampling PDFs for infinite area light

    // Compute scalar-valued image _img_ from environment map
    int width = 2 * resolution.x, height = 2 * resolution.y;
    std::unique_ptr<Float[]> img(new Float[width * height]);


    for (int v = 0; v < height; v++){
        float vp = (float)v / (float)height;
        float sinTheta = std::sin(Pi * float(v + .5f) / (float)height);
        for(int u = 0; u < width; u++){
            float up = (float)u / (float)width;
            Point2f temp = Point2f(up, vp);
            img[u + v * width] = Lookup(temp).y;
            img[u + v * width] *= sinTheta;
        }
    }


    // Compute sampling distributions for rows and columns of image
    distribution.reset(new Distribution2D(img.get(), width, height));


    worldCenter = t.position();
    worldRadius = t.getScale().x;
}

Color3f InfiniteAreaLight::L(const Intersection &isect, const Vector3f &woW) const{

    Vector3f w = glm::normalize(glm::vec3(transform.invT() * glm::vec4(woW,0.f)));
    Point2f st(SphericalPhi(w) * Inv2Pi, SphericalTheta(w) * InvPi);
    return glm::dot(isect.normalGeometric, w) < 0.f ? (lightColor * Color3f(Lookup(st))) : Color3f(0.f);
}


Color3f InfiniteAreaLight::Power() const {
    Point2f temp = Point2f(.5f, .5f);
    return Pi * worldRadius * worldRadius *
           lightColor * Color3f(Lookup(temp));
}

Color3f InfiniteAreaLight::Le(const Ray &ray) const {
    Vector3f w = glm::normalize(glm::vec3(transform.invT() * glm::vec4(ray.direction,0.f)));
    Point2f st(SphericalPhi(w) * Inv2Pi, SphericalTheta(w) * InvPi);
    return lightColor * Color3f(Lookup(st));
}


Color3f InfiniteAreaLight::Sample_Li(const Intersection &ref, const Point2f &u,
                                      Vector3f *wi, Float *pdf) const {
    //ProfilePhase _(Prof::LightSample);
    // Find $(u,v)$ sample coordinates in infinite light texture
    Float mapPdf;
    Point2f uv = distribution->SampleContinuous(u, &mapPdf);
    if (mapPdf == 0) return Color3f(0.f);

    // Convert infinite light sample point to direction
    Float theta = uv[1] * Pi, phi = uv[0] * 2 * Pi;
    Float cosTheta = std::cos(theta), sinTheta = std::sin(theta);
    Float sinPhi = std::sin(phi), cosPhi = std::cos(phi);
    *wi = glm::vec3(transform.T() * glm::vec4(Vector3f(sinTheta * cosPhi, cosTheta, sinTheta * sinPhi), 0.f));


    // Compute PDF for sampled infinite light direction
    *pdf = mapPdf / (2 * Pi * Pi * sinTheta);
    if (sinTheta == 0) *pdf = 0;


    return lightColor * Color3f(Lookup(uv));
}

float InfiniteAreaLight::Pdf_Li(const Intersection &, const Vector3f &w) const {

    Vector3f wi = glm::mat3(transform.invT()) * w;
    wi = glm::normalize(wi);
    Float theta = SphericalTheta(wi), phi = SphericalPhi(wi);
    Float sinTheta = std::sin(theta);
    if (sinTheta == 0) return 0;
    return distribution->Pdf(Point2f(phi * Inv2Pi, theta * InvPi)) /
           (2 * Pi * Pi * sinTheta);
}


Ray InfiniteAreaLight::generatePhotonRay(const Point2f &u, Sampler *sampler, float *pdf){

    Float mapPdf;
    Point2f uv = distribution->SampleContinuous(u, &mapPdf);
    if (mapPdf == 0) {tempPower = Color3f(0.f);}


    // Convert infinite light sample point to direction
    Float theta = uv[1] * Pi, phi = uv[0] * 2 * Pi;
    Float cosTheta = std::cos(theta), sinTheta = std::sin(theta);
    Float sinPhi = std::sin(phi), cosPhi = std::cos(phi);
    Vector3f dir = glm::vec3(transform.T() * glm::vec4(Vector3f(sinTheta * cosPhi, cosTheta, sinTheta * sinPhi), 0.f));

    Point3f pointOnSphere = dir;
    dir = -dir;

    pointOnSphere = glm::vec3(transform.T() * glm::vec4(pointOnSphere, 1.0f));


    // Compute PDF for sampled infinite light direction
    *pdf = mapPdf / (2 * Pi * Pi * sinTheta);
    if (sinTheta == 0) *pdf = 0;

    tempPower = lightColor * Color3f(Lookup(uv));

    return Ray(pointOnSphere, glm::normalize(dir));
}

Color3f InfiniteAreaLight::GetInitialPhotonPower() const{
   return tempPower;
}


//Color3f InfiniteAreaLight::Sample_Le(const Point2f &u1, const Point2f &u2,
//                                      Float time, Ray *ray, Normal3f *nLight,
//                                      Float *pdfPos, Float *pdfDir) const {
//    ProfilePhase _(Prof::LightSample);
//    // Compute direction for infinite light sample ray
//    Point2f u = u1;

//    // Find $(u,v)$ sample coordinates in infinite light texture
//    Float mapPdf;
//    Point2f uv = distribution->SampleContinuous(u, &mapPdf);
//    if (mapPdf == 0) return Color3f(0.f);
//    Float theta = uv[1] * Pi, phi = uv[0] * 2.f * Pi;
//    Float cosTheta = std::cos(theta), sinTheta = std::sin(theta);
//    Float sinPhi = std::sin(phi), cosPhi = std::cos(phi);
//    Vector3f d =
//        -LightToWorld(Vector3f(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta));
//    *nLight = (Normal3f)d;

//    // Compute origin for infinite light sample ray
//    Vector3f v1, v2;
//    CoordinateSystem(-d, &v1, &v2);
//    Point2f cd = ConcentricSampleDisk(u2);
//    Point3f pDisk = worldCenter + worldRadius * (cd.x * v1 + cd.y * v2);
//    *ray = Ray(pDisk + worldRadius * -d, d, Infinity, time);

//    // Compute _InfiniteAreaLight_ ray PDFs
//    *pdfDir = sinTheta == 0 ? 0 : mapPdf / (2 * Pi * Pi * sinTheta);
//    *pdfPos = 1 / (Pi * worldRadius * worldRadius);
//    return Color3f(Lmap->Lookup(uv), Color3fType::Illuminant);
//}

//void InfiniteAreaLight::Pdf_Le(const Ray &ray, const Normal3f &, Float *pdfPos,
//                               Float *pdfDir) const {
//    ProfilePhase _(Prof::LightPdf);
//    Vector3f d = -WorldToLight(ray.d);
//    Float theta = SphericalTheta(d), phi = SphericalPhi(d);
//    Point2f uv(phi * Inv2Pi, theta * InvPi);
//    Float mapPdf = distribution->Pdf(uv);
//    *pdfDir = mapPdf / (2 * Pi * Pi * std::sin(theta));
//    *pdfPos = 1 / (Pi * worldRadius * worldRadius);
//}



//std::shared_ptr<InfiniteAreaLight> CreateInfiniteLight(
//    const Transform &light2world, const ParamSet &paramSet) {
//    Color3f L = paramSet.FindOneColor3f("L", Color3f(1.0));
//    Color3f sc = paramSet.FindOneColor3f("scale", Color3f(1.0));
//    std::string texmap = paramSet.FindOneFilename("mapname", "");
//    int nSamples = paramSet.FindOneInt("samples",
//                                       paramSet.FindOneInt("nsamples", 1));
//    if (PbrtOptions.quickRender) nSamples = std::max(1, nSamples / 4);
//    return std::make_shared<InfiniteAreaLight>(light2world, L * sc, nSamples,
//                                               texmap);
//}


Color3f InfiniteAreaLight::Lookup(Point2f &uv_coord) const
{
    if(texmap)
    {
        int X = glm::min(texmap->width() * uv_coord.x, texmap->width() - 1.0f);
        int Y = glm::min(texmap->height() * (1.0f - uv_coord.y), texmap->height() - 1.0f);
        QColor color = texmap->pixel(X, Y);
        return Color3f(color.red(), color.green(), color.blue())/255.0f;
    }
    return Color3f(1.f, 1.f, 1.f);
}

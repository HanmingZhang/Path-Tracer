#include "shape.h"
#include <QDateTime>

pcg32 Shape::colorRNG = pcg32(QDateTime::currentMSecsSinceEpoch());


void Shape::InitializeIntersection(Intersection *isect, float t, Point3f pLocal) const
{
    isect->point = Point3f(transform.T() * glm::vec4(pLocal, 1));
    ComputeTBN(pLocal, &(isect->normalGeometric), &(isect->tangent), &(isect->bitangent));
    isect->uv = GetUVCoordinates(pLocal);
    isect->t = t;
}

Intersection Shape::Sample(const Intersection &ref, const Point2f &xi, float *pdf) const
{
    //TODO
    // First, we invoke two-input Sample of subclasses
    Intersection it = Sample(xi,pdf);

    // if pdf is not zero, we need to convert
    // from a PDF with respect to surface area
    // to a PDF with respect to solid angle at the reference point of intersection

    // PDFa = 1 / Area;
    // PDFshape = (r * r / cos(theta)) * PDFa

    if(*pdf != 0.f){
        float distanceSquared = glm::distance2(it.point, ref.point);
        *pdf *= distanceSquared / AbsDot(it.normalGeometric, glm::normalize(ref.point - it.point));
    }


    return it;
}

float Shape::Pdf(const Intersection &ref, const Vector3f &wi) const{
    // Intersect sample ray with area light geometry
    Ray ray = ref.SpawnRay(wi);
    Intersection isectLight;
    if(!Intersect(ray, &isectLight)) return 0;

    // Convert light sample weight to solid angle measure
    float pdf = glm::distance2(ref.point, isectLight.point) /
                (AbsDot(isectLight.normalGeometric, -wi) * (Area() * 2.0));

    return pdf;
}

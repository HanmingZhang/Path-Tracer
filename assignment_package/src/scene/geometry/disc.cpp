#include "disc.h"

float Disc::Area() const
{
    //TODO

    // Actually, it should be something like this
    // return phiMax * 0.5 * (radius * radius - innerRadius * innerRadius);

    glm::vec3 scale = transform.getScale();

    float scaleX = scale.x;

    float scaleY = scale.y;

    return Pi * scaleX * scaleY;
}

bool Disc::Intersect(const Ray &ray, Intersection *isect) const
{
    //Transform the ray
    Ray r_loc = ray.GetTransformedCopy(transform.invT());

    //Ray-plane intersection
    float t = glm::dot(glm::vec3(0,0,1), (glm::vec3(0.5f, 0.5f, 0) - r_loc.origin)) / glm::dot(glm::vec3(0,0,1), r_loc.direction);
    Point3f P = Point3f(t * r_loc.direction + r_loc.origin);
    //Check that P is within the bounds of the disc (not bothering to take the sqrt of the dist b/c we know the radius)
    float dist2 = (P.x * P.x + P.y * P.y);
    if(t > 0 && dist2 <= 1.f)
    {
        InitializeIntersection(isect, t, P);
        return true;
    }
    return false;
}

void Disc::ComputeTBN(const Point3f &P, Normal3f *nor, Vector3f *tan, Vector3f *bit) const
{
    *nor = glm::normalize(transform.invTransT() * Normal3f(0,0,1));
    //TODO: Compute tangent and bitangent

    Vector3f tangent = glm::normalize(transform.T3() * glm::vec3(1,0,0));
    *tan = tangent;

    Vector3f bitangent = glm::normalize(transform.T3() * glm::vec3(0,1,0));
    *bit = bitangent;


}


Point2f Disc::GetUVCoordinates(const Point3f &point) const
{
    return glm::vec2((point.x + 1)/2.f, (point.y + 1)/2.f);
}


Intersection Disc::Sample(const Point2f &xi, Float *pdf) const{


    Point2f pd = Point2f(WarpFunctions::squareToDiskConcentric(xi));

    //now, our disc's height is 0 and ridius is 1 based on its creation func in glshapecreation.cpp
    float height = 0.f;
    float radius = 1.0f;

    Point3f pObj = glm::vec3(pd.x * radius, pd.y * radius, height);


    Intersection it;
    it.normalGeometric = glm::vec3(glm::normalize(transform.T() * glm::vec4(0.f, 0.f, 1.f,0.f)));
    it.point = glm::vec3(transform.T() * glm::vec4(pObj, 1.0f));

    *pdf = 1.0f / Area();

    return it;
}

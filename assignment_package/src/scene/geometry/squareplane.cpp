#include "squareplane.h"

float SquarePlane::Area() const
{
    //TODO

    glm::vec3 scale = transform.getScale();

    float scaleX = scale.x;

    float scaleY = scale.y;

    return scaleX * scaleY;
}

bool SquarePlane::Intersect(const Ray &ray, Intersection *isect) const
{
    //Transform the ray
    Ray r_loc = ray.GetTransformedCopy(transform.invT());

    //Ray-plane intersection
    float t = glm::dot(glm::vec3(0,0,1), (glm::vec3(0.5f, 0.5f, 0) - r_loc.origin)) / glm::dot(glm::vec3(0,0,1), r_loc.direction);
    Point3f P = Point3f(t * r_loc.direction + r_loc.origin);
    //Check that P is within the bounds of the square
    if(t > 0 && P.x >= -0.5f && P.x <= 0.5f && P.y >= -0.5f && P.y <= 0.5f)
    {
        InitializeIntersection(isect, t, P);
        return true;
    }
    return false;
}

void SquarePlane::ComputeTBN(const Point3f &P, Normal3f *nor, Vector3f *tan, Vector3f *bit) const
{
    *nor = glm::normalize(transform.invTransT() * Normal3f(0,0,1));
    //TODO: Compute tangent and bitangent

    Vector3f tangent = glm::normalize(transform.T3() * glm::vec3(1,0,0));
    *tan = tangent;

    Vector3f bitangent = glm::normalize(transform.T3() * glm::vec3(0,1,0));
    *bit = bitangent;
}


Point2f SquarePlane::GetUVCoordinates(const Point3f &point) const
{
    return Point2f(point.x + 0.5f, point.y + 0.5f);
}

Intersection SquarePlane::Sample(const Point2f &xi, Float *pdf) const{

    Point3f pObj = glm::vec3(xi.x - 0.5f, xi.y - 0.5f, 0.f);

    Intersection it;
    it.normalGeometric = glm::vec3(glm::normalize(transform.T() * glm::vec4(0.f, 0.f, 1.f,0.f)));
    it.point = glm::vec3(transform.T() * glm::vec4(pObj, 1.0f));

    *pdf = 1.0f / Area();

    return it;
}

Bounds3f SquarePlane::WorldBound() const{

    Bounds3f local_bounding_box(Point3f(-0.5f, -0.5f, -0.01f),
                                Point3f(0.5f, 0.5f, 0.01f));

    return local_bounding_box.Apply(transform);

}

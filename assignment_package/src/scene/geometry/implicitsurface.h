#pragma once
#include <scene/geometry/shape.h>

static const float Epsilon = 0.01f;
static const int MAX_MARCHING_STEPS = 1000;



class ImplicitSurface : public Shape
{
public:
    virtual bool Intersect(const Ray &ray, Intersection *isect) const;
    virtual Point2f GetUVCoordinates(const Point3f &point) const;
    virtual void ComputeTBN(const Point3f& P, Normal3f* nor, Vector3f* tan, Vector3f* bit) const;

    virtual float Area() const;

    virtual Intersection Sample(const Point2f &xi, Float *pdf) const;

    virtual Bounds3f WorldBound() const;

    void create();

private:

    float sceneSDF(Point3f p) const;
};



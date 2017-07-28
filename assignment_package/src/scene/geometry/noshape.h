#pragma once

#include <scene/geometry/shape.h>


class NoShape : public Shape
{
public:
    virtual bool Intersect(const Ray &ray, Intersection *isect) const {return false;}

    virtual Point2f GetUVCoordinates(const Point3f &point) const {return Point2f(0.f);}

    virtual void ComputeTBN(const Point3f& P, Normal3f* nor, Vector3f* tan, Vector3f* bit) const {}

    virtual float Area() const {return 0.f;}

    virtual Intersection Sample(const Point2f &xi, Float *pdf) const {return Intersection();}

    virtual Bounds3f WorldBound() const {return Bounds3f();}

    void create();
};

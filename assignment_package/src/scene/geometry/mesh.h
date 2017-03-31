#pragma once

#include "shape.h"
#include <QList>

class Triangle : public Shape
{
public:
    Triangle(const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3);
    Triangle(const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3, const glm::vec3 &n1, const glm::vec3 &n2, const glm::vec3 &n3);
    Triangle(const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3, const glm::vec3 &n1, const glm::vec3 &n2, const glm::vec3 &n3, const glm::vec2 &t1, const glm::vec2 &t2, const glm::vec2 &t3);


    virtual bool Intersect(const Ray &ray, Intersection *isect) const;
    virtual Point2f GetUVCoordinates(const Point3f &point) const;
    virtual void ComputeTBN(const Point3f& P, Normal3f* nor, Vector3f* tan, Vector3f* bit) const;
    virtual void ComputeTriangleTBN(const Point3f& P, Normal3f* nor, Vector3f* tan, Vector3f* bit, const Point2f& uv) const;

    virtual float Area() const;
    virtual void InitializeIntersection(Intersection* isect, float t, Point3f pLocal) const;

    // Sample a point on the surface of the shape and return the PDF with
    // respect to area on the surface.
    virtual Intersection Sample(const Point2f &xi, Float *pdf) const;

    Point3f points[3];
    Normal3f normals[3];
    Point2f uvs[3];
    Normal3f planeNormal;

    void create();//This does nothing because individual triangles are not rendered with OpenGL;
                            //they are rendered all together in their Mesh.

    Normal3f GetNormal(const Point3f &position) const;
};

//A mesh just holds a collection of triangles against which one can test intersections.
//Its primary purpose is to store VBOs for rendering the triangles in OpenGL.
class Mesh : public Shape
{
public:
    virtual bool Intersect(const Ray &ray, Intersection *isect) const;
    virtual Point2f GetUVCoordinates(const Point3f &point) const;
    virtual void ComputeTBN(const Point3f& P, Normal3f* nor, Vector3f* tan, Vector3f* bit) const;
    // Should return the sum of all triangles' areas
    virtual float Area() const;
    virtual void InitializeIntersection(Intersection* isect, float t, Point3f pLocal) const;

    // Sample a point on the surface of the shape and return the PDF with
    // respect to area on the surface.
    virtual Intersection Sample(const Point2f &xi, Float *pdf) const;

    void create();

    void LoadOBJ(const QStringRef &filename, const QStringRef &local_path);

private:
    QList<Triangle*> faces;
};







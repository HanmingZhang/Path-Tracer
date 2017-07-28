#include "implicitsurface.h"

float ImplicitSurface::Area() const
{
    return 1.0f;
}


void ImplicitSurface::ComputeTBN(const Point3f &P, Normal3f *nor, Vector3f *tan, Vector3f *bit) const
{

    Normal3f n = glm::normalize(glm::vec3(sceneSDF(Point3f(P.x + Epsilon, P.y, P.z)) - sceneSDF(Point3f(P.x - Epsilon, P.y, P.z)),
                                          sceneSDF(Point3f(P.x, P.y + Epsilon, P.z)) - sceneSDF(Point3f(P.x, P.y - Epsilon, P.z)),
                                          sceneSDF(Point3f(P.x, P.y, P.z + Epsilon)) - sceneSDF(Point3f(P.x, P.y, P.z - Epsilon))));


    *nor = glm::normalize(transform.invTransT() * n);
    CoordinateSystem((Vector3f)(*nor), tan, bit);

    return;
}


bool ImplicitSurface::Intersect(const Ray &ray, Intersection *isect) const
{
    //Transform the ray
    Ray r_loc = ray.GetTransformedCopy(transform.invT());


    float start = 0.f;
    float end = 100000.f;

    float depth = start;

    for(int i = 0; i < MAX_MARCHING_STEPS; i++){

        Point3f marchingPoint = r_loc.origin + depth * r_loc.direction;

        float dist = sceneSDF(marchingPoint);

        if(dist < Epsilon){
            // We're inside the scene surface!

            // If it's inside, we need to move it a little bit outside
            // considering shadow test!
            while(dist < 0.f){
                marchingPoint = marchingPoint - 0.1f * r_loc.direction;
                dist = sceneSDF(marchingPoint);
            }

            InitializeIntersection(isect, depth, marchingPoint);
            return true;
        }

        //depth += dist;

        // How about using a naive way.
        depth += 0.05f;

        if(depth >= end){
            // Gone too far; give up!
            return false;
        }

    }
    return false;
}



Point2f ImplicitSurface::GetUVCoordinates(const Point3f &point) const
{

    return Point2f(point);
}


Intersection ImplicitSurface::Sample(const Point2f &xi, Float *pdf) const
{

    return Intersection();
}


Bounds3f ImplicitSurface::WorldBound() const
{
   Bounds3f local_bounding_box(Point3f(-4.0f, -4.0f, -4.0f),
                               Point3f( 4.0f,  4.0f,  4.0f));

   return local_bounding_box.Apply(transform);

}


float ImplicitSurface::sceneSDF(Point3f p) const
{
     float x = p.x;
     float y = p.y;
     float z = p.z;

     //Sphere
//     return sqrt(x * x + y * y + z * z) - 1.0f;


     //Heart
//     float temp = x * x + 2.25f * y * y + z * z - 1.0f;

//     float temp2 = temp * temp * temp - x * x * z * z * z - 9.0f / 80.f * y * y * z * z * z;

//     return (temp2 < 0.f) ? -pow(-temp2, 1.0 / 6.0) : pow(temp2, 1.0 / 6.0);


     //Kiss Surface
//     float temp = x * x + y * y - (1.0f - z) * z * z * z * z;
//     return pow(temp, 1.0f / 5.0f);

     //Mobius Strip
//     float R = 1.0f;

//     float temp =  -R * R * y + x * x * y + y * y * y - 2.0f * R * x * z
//             -2.0f * x * x * z - 2.0f * y * y * z + y * z * z;

//     return (temp < 0.f) ? -pow(-temp, 1.0 / 3.0) : pow(temp, 1.0 / 3.0);


     //Chair Surface
//     float k = 5.0f;
//     float a = 0.95f;
//     float b = 0.8f;

//     float temp = x * x + y * y + z * z - a * k * k;
//     float temp2 = temp * temp - b * ((z - k) * (z - k) - 2.0f * x * x) * ((z + k) * (z + k) - 2.0f * y * y);

//     //return temp2;
//     return (temp2 < 0.f) ? -pow(-temp2, 1.0 / 16.0) : pow(temp2, 1.0 / 16.0);

     //Tanglecube
     float temp = x * x * x * x- 5.0f * x * x + y * y * y * y - 5.0f * y * y
                + z * z * z * z - 5.0f * z * z + 11.8f;

     return (temp < 0.f) ? -pow(-temp, 1.0 / 8.0) : pow(temp, 1.0 / 8.0);
}

#include "bounds.h"

bool Bounds3f::Intersect(const Ray &r, float* t) const
{
    //TODO
    Ray r_loc = r;

    float t_n = -1000000.0f;
    float t_f = 1000000.f;
    for(int i = 0; i < 3; i++){

        float minValue, maxValue;

        if(i == 0){
            minValue = min.x;
            maxValue = max.x;
        }
        if(i == 1){
            minValue = min.y;
            maxValue = max.y;
        }
        if(i == 2){
            minValue = min.z;
            maxValue = max.z;
        }


        //Ray parallel to slab check
        if(r_loc.direction[i] == 0){
            if(r_loc.origin[i] < minValue || r_loc.origin[i] > maxValue){
                return false;
            }
        }
        //If not parallel, do slab intersect check
        float t0 = (minValue - r_loc.origin[i])/r_loc.direction[i];
        float t1 = (maxValue - r_loc.origin[i])/r_loc.direction[i];
        if(t0 > t1){
            float temp = t1;
            t1 = t0;
            t0 = temp;
        }
        if(t0 > t_n){
            t_n = t0;
        }
        if(t1 < t_f){
            t_f = t1;
        }
    }
    if(t_n < t_f)
    {
        float t_result = t_n > 0 ? t_n : t_f;
        if(t_result < 0){

            // negative t values are valid if and only if the ray's origin lies within the bounding box
            if(Inside(r_loc.origin)){
                (*t) = t_result;
                return true;
            }

            return false;
        }

        (*t) = t_result;
        return true;
    }
    else{
        //If t_near was greater than t_far, we did not hit the cube
        return false;
    }
}

Bounds3f Bounds3f::Apply(const Transform &tr)
{
    //TODO
    std::vector<Point3f> vertexPos;


    vertexPos.push_back(Point3f(min.x, min.y, min.z));
    vertexPos.push_back(Point3f(min.x, min.y, max.z));
    vertexPos.push_back(Point3f(min.x, max.y, min.z));
    vertexPos.push_back(Point3f(min.x, max.y, max.z));
    vertexPos.push_back(Point3f(max.x, min.y, min.z));
    vertexPos.push_back(Point3f(max.x, min.y, max.z));
    vertexPos.push_back(Point3f(max.x, max.y, min.z));
    vertexPos.push_back(Point3f(max.x, max.y, max.z));


    float new_min_x = 9999999.0f;
    float new_min_y = 9999999.0f;
    float new_min_z = 9999999.0f;

    float new_max_x = -9999999.0f;
    float new_max_y = -9999999.0f;
    float new_max_z = -9999999.0f;

    for(Point3f each : vertexPos){

        Point3f newVertexPos = glm::vec3(tr.T() * glm::vec4(each, 1.0f));

        // process new X pos
        if(newVertexPos.x < new_min_x){
            new_min_x = newVertexPos.x;
        }
        if(newVertexPos.x > new_max_x){
            new_max_x = newVertexPos.x;
        }

        // process new Y pos
        if(newVertexPos.y < new_min_y){
            new_min_y = newVertexPos.y;
        }
        if(newVertexPos.y > new_max_y){
            new_max_y = newVertexPos.y;
        }

        // process new Z pos
        if(newVertexPos.z < new_min_z){
            new_min_z = newVertexPos.z;
        }
        if(newVertexPos.z > new_max_z){
            new_max_z = newVertexPos.z;
        }
    }

    return Bounds3f(Point3f(new_min_x, new_min_y, new_min_z),
                    Point3f(new_max_x, new_max_y, new_max_z));
}

float Bounds3f::SurfaceArea() const
{
    //TODO
    Vector3f d = Diagonal();

    return 2.0f * (d.x * d.y + d.x * d.z + d.y * d.z);
}

Bounds3f Union(const Bounds3f& b1, const Bounds3f& b2)
{
    if(b1.min.x == 0.f && b1.max.x == 0.f &&
       b1.min.y == 0.f && b1.max.y == 0.f &&
       b1.min.z == 0.f && b1.max.z == 0.f){
        return b2;
    }


    else{
        return Bounds3f(Point3f(std::min(b1.min.x, b2.min.x),
                                std::min(b1.min.y, b2.min.y),
                                std::min(b1.min.z, b2.min.z)),
                        Point3f(std::max(b1.max.x, b2.max.x),
                                std::max(b1.max.y, b2.max.y),
                                std::max(b1.max.z, b2.max.z)));
    }
}

Bounds3f Union(const Bounds3f& b1, const Point3f& p)
{
    return Bounds3f(Point3f(std::min(b1.min.x, p.x),
                            std::min(b1.min.y, p.y),
                            std::min(b1.min.z, p.z)),
                    Point3f(std::max(b1.max.x, p.x),
                            std::max(b1.max.y, p.y),
                            std::max(b1.max.z, p.z)));
}

Bounds3f Union(const Bounds3f& b1, const glm::vec4& p)
{
    return Union(b1, Point3f(p));
}



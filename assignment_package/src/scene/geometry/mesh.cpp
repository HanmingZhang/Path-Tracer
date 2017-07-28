#include <scene/geometry/mesh.h>
#include <la.h>
#include <tinyobj/tiny_obj_loader.h>
#include <iostream>

Bounds3f Triangle::WorldBound() const
{

//    // Transform to local space
//    Point3f point0 = glm::vec3(transform.invT() * glm::vec4(points[0], 1.0f));
//    Point3f point1 = glm::vec3(transform.invT() * glm::vec4(points[1], 1.0f));
//    Point3f point2 = glm::vec3(transform.invT() * glm::vec4(points[2], 1.0f));

    float minX = std::min(points[0].x, std::min(points[1].x, points[2].x));
    float minY = std::min(points[0].y, std::min(points[1].y, points[2].y));
    float minZ = std::min(points[0].z, std::min(points[1].z, points[2].z));

    float maxX = std::max(points[0].x, std::max(points[1].x, points[2].x));
    float maxY = std::max(points[0].y, std::max(points[1].y, points[2].y));
    float maxZ = std::max(points[0].z, std::max(points[1].z, points[2].z));


    if(minX == maxX){
        std::cout << "triange bounding box x is the same!" << std::endl;
        minX -= 0.01f;
        maxX += 0.01f;
    }

    if(minY == maxY){
        std::cout << "triange bounding box y is the same!" << std::endl;
        minY -= 0.01f;
        maxY += 0.01f;
    }

    if(minZ == maxZ){
        std::cout << "triange bounding box z is the same!" << std::endl;
        minZ -= 0.01f;
        maxZ += 0.01f;
    }


    return Bounds3f(Point3f(minX, minY, minZ),
                    Point3f(maxX, maxY, maxZ));

//    Bounds3f local_bounding_box(Point3f(minX, minY, minZ),
//                                Point3f(maxX, maxY, maxZ));

//    return local_bounding_box.Apply(transform);

}

float Triangle::Area() const
{
    return glm::length(glm::cross(points[0] - points[1], points[2] - points[1])) * 0.5f;
}

Triangle::Triangle(const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3):
    Triangle(p1, p2, p3, glm::vec3(1), glm::vec3(1), glm::vec3(1), glm::vec2(0), glm::vec2(0), glm::vec2(0))
{
    for(int i = 0; i < 3; i++)
    {
        normals[i] = planeNormal;
    }
}


Triangle::Triangle(const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3, const glm::vec3 &n1, const glm::vec3 &n2, const glm::vec3 &n3):
    Triangle(p1, p2, p3, n1, n2, n3, glm::vec2(0), glm::vec2(0), glm::vec2(0))
{}


Triangle::Triangle(const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3, const glm::vec3 &n1, const glm::vec3 &n2, const glm::vec3 &n3, const glm::vec2 &t1, const glm::vec2 &t2, const glm::vec2 &t3)
{
    planeNormal = glm::normalize(glm::cross(p2 - p1, p3 - p2));
    points[0] = p1;
    points[1] = p2;
    points[2] = p3;
    normals[0] = n1;
    normals[1] = n2;
    normals[2] = n3;
    uvs[0] = t1;
    uvs[1] = t2;
    uvs[2] = t3;
}

float TriArea(const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3)
{
    return glm::length(glm::cross(p1 - p2, p3 - p2)) * 0.5f;
}

//Returns the interpolation of the triangle's three normals based on the point inside the triangle that is given.
Normal3f Triangle::GetNormal(const Point3f &P) const
{
    float A = TriArea(points[0], points[1], points[2]);
    float A0 = TriArea(points[1], points[2], P);
    float A1 = TriArea(points[0], points[2], P);
    float A2 = TriArea(points[0], points[1], P);
    return glm::normalize(normals[0] * A0/A + normals[1] * A1/A + normals[2] * A2/A);
}


bool Triangle::Intersect(const Ray& r, Intersection* isect) const
{

    //1. Ray-plane intersection
    float t =  glm::dot(planeNormal, (points[0] - r.origin)) / glm::dot(planeNormal, r.direction);
    if(t < 0) return false;

    glm::vec3 P = r.origin + t * r.direction;
    //2. Barycentric test
    float S = 0.5f * glm::length(glm::cross(points[0] - points[1], points[0] - points[2]));
    float s1 = 0.5f * glm::length(glm::cross(P - points[1], P - points[2]))/S;
    float s2 = 0.5f * glm::length(glm::cross(P - points[2], P - points[0]))/S;
    float s3 = 0.5f * glm::length(glm::cross(P - points[0], P - points[1]))/S;
    float sum = s1 + s2 + s3;

    if(s1 >= 0 && s1 <= 1 && s2 >= 0 && s2 <= 1 && s3 >= 0 && s3 <= 1 && fequal(sum, 1.0f)){
        isect->t = t;
        InitializeIntersection(isect, t, Point3f(P));
        return true;
    }
    return false;
}

void Triangle::InitializeIntersection(Intersection *isect, float t, Point3f pLocal) const
{
    isect->point = pLocal;
    isect->uv = GetUVCoordinates(pLocal);
    ComputeTriangleTBN(pLocal, &(isect->normalGeometric), &(isect->tangent), &(isect->bitangent), isect->uv);
    isect->t = t;
}

void Triangle::ComputeTBN(const Point3f &P, Normal3f *nor, Vector3f *tan, Vector3f *bit) const
{
    //Triangle uses ComputeTriangleTBN instead of this function.

    ComputeTriangleTBN(P, nor, tan, bit, GetUVCoordinates(P));
}

void Triangle::ComputeTriangleTBN(const Point3f &P, Normal3f *nor, Vector3f *tan, Vector3f *bit, const Point2f &uv) const
{
    *nor = GetNormal(P);
    //TODO: Compute tangent and bitangent based on UV coordinates.


    Point3f pos0 = P;
    Point3f pos1 = points[1];
    Point3f pos2 = points[2];

    Point2f uv0 = uv;
    Point2f uv1 = uvs[1];
    Point2f uv2 = uvs[2];


    // Edges of the triangle : postion delta
    Vector3f deltaPos1 = pos1 - pos0;
    Vector3f deltaPos2 = pos2 - pos0;;

    // UV delta
    Vector2f deltaUV1 = uv1-uv0;
    Vector2f deltaUV2 = uv2-uv0;


    float r = 1.0f / (deltaUV1.x * deltaUV2.y - deltaUV1.y * deltaUV2.x);
    Vector3f tangent = (deltaPos1 * deltaUV2.y   - deltaPos2 * deltaUV1.y)*r;
    Vector3f bitangent = (deltaPos2 * deltaUV1.x   - deltaPos1 * deltaUV2.x)*r;

    tangent = glm::normalize(tangent);
    bitangent = glm::normalize(bitangent);

    if(std::isnan(tangent.x) ||
       std::isnan(tangent.y) ||
       std::isnan(tangent.z)){
//        std::cout << "TBN triange stop here" << std::endl;
        CoordinateSystem((*nor), tan, bit);
    }
    if(std::isnan((*nor).x) ||
       std::isnan((*nor).y) ||
       std::isnan((*nor).z)){
//        std::cout << "TBN triange stop here" << std::endl;
        *nor = planeNormal;
        CoordinateSystem((*nor), tan, bit);
    }
    if(std::isnan(bitangent.x) ||
       std::isnan(bitangent.y) ||
       std::isnan(bitangent.z)){
//        std::cout << "TBN triange stop here" << std::endl;
        CoordinateSystem((*nor), tan, bit);
    }


    *tan = tangent;

    *bit = bitangent;


}


Intersection Triangle::Sample(const Point2f &xi, Float *pdf) const
{
    //TODO for extra credit
    return Intersection();
}


Point2f Triangle::GetUVCoordinates(const Point3f &point) const
{
    float A = TriArea(points[0], points[1], points[2]);
    float A0 = TriArea(points[1], points[2], point);
    float A1 = TriArea(points[0], points[2], point);
    float A2 = TriArea(points[0], points[1], point);
    return uvs[0] * A0/A + uvs[1] * A1/A + uvs[2] * A2/A;
}

void Mesh::LoadOBJ(const QStringRef &filename, const QStringRef &local_path, const Transform &transform)
{
    QString filepath = local_path.toString(); filepath.append(filename);
    std::vector<tinyobj::shape_t> shapes; std::vector<tinyobj::material_t> materials;
    std::string errors = tinyobj::LoadObj(shapes, materials, filepath.toStdString().c_str());
    std::cout << errors << std::endl;
    if(errors.size() == 0)
    {
        //Read the information from the vector of shape_ts
        for(unsigned int i = 0; i < shapes.size(); i++)
        {
            std::vector<float> &positions = shapes[i].mesh.positions;
            std::vector<float> &normals = shapes[i].mesh.normals;
            std::vector<float> &uvs = shapes[i].mesh.texcoords;
            std::vector<unsigned int> &indices = shapes[i].mesh.indices;
            for(unsigned int j = 0; j < indices.size(); j += 3)
            {
                glm::vec3 p1 = glm::vec3(transform.T() * glm::vec4(positions[indices[j]*3], positions[indices[j]*3+1], positions[indices[j]*3+2], 1));
                glm::vec3 p2 = glm::vec3(transform.T() * glm::vec4(positions[indices[j+1]*3], positions[indices[j+1]*3+1], positions[indices[j+1]*3+2], 1));
                glm::vec3 p3 = glm::vec3(transform.T() * glm::vec4(positions[indices[j+2]*3], positions[indices[j+2]*3+1], positions[indices[j+2]*3+2], 1));

                auto t = std::make_shared<Triangle>(p1, p2, p3);
                if(normals.size() > 0)
                {
                    glm::vec3 n1 = transform.invTransT() * glm::vec3(normals[indices[j]*3], normals[indices[j]*3+1], normals[indices[j]*3+2]);
                    glm::vec3 n2 = transform.invTransT() * glm::vec3(normals[indices[j+1]*3], normals[indices[j+1]*3+1], normals[indices[j+1]*3+2]);
                    glm::vec3 n3 = transform.invTransT() * glm::vec3(normals[indices[j+2]*3], normals[indices[j+2]*3+1], normals[indices[j+2]*3+2]);
                    t->normals[0] = n1;
                    t->normals[1] = n2;
                    t->normals[2] = n3;
                }
                if(uvs.size() > 0)
                {
                    glm::vec2 t1(uvs[indices[j]*2], uvs[indices[j]*2+1]);
                    glm::vec2 t2(uvs[indices[j+1]*2], uvs[indices[j+1]*2+1]);
                    glm::vec2 t3(uvs[indices[j+2]*2], uvs[indices[j+2]*2+1]);
                    t->uvs[0] = t1;
                    t->uvs[1] = t2;
                    t->uvs[2] = t3;
                }
                this->faces.append(t);
            }
        }
        std::cout << "" << std::endl;
        //TODO: .mtl file loading
    }
    else
    {
        //An error loading the OBJ occurred!
        std::cout << errors << std::endl;
    }
}

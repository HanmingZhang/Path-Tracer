#include "thinlenperspectivecamera.h"
#include <warpfunctions.h>

ThinLenPerspectiveCamera::ThinLenPerspectiveCamera():
    Camera(400, 400)
{
    look = Vector3f(0,0,-1);
    up = Vector3f(0,1,0);
    right = Vector3f(1,0,0);

    lensRadius = 5.f;
    focalDistance = glm::length(ref - eye);
}


ThinLenPerspectiveCamera::ThinLenPerspectiveCamera(const Camera &c):
    Camera(c)
{
    lensRadius = 5.f;
    focalDistance = glm::length(ref - eye);
}

Ray ThinLenPerspectiveCamera::Raycast(const Point2f &pt) const{

    float x = pt.x;
    float y = pt.y;

    // From screen space to Camera Space
    float ndc_x = (2.f*x/width - 1);
    float ndc_y = (1 - 2.f*y/height);
    Point3f ndc_point(ndc_x, ndc_y, 0.f);


    // Set this ray in Camera Space
    glm::mat4 persp_mat = glm::perspective(glm::radians(fovy), width / (float)height, near_clip, far_clip);

    glm::mat4 inv_persp_mat = glm::inverse(persp_mat);

    glm::vec4 camera_point_vec4 = inv_persp_mat * glm::vec4(ndc_point, 1.f);

    glm::vec3 camera_point = glm::vec3(camera_point_vec4);

    Ray r = Ray(Point3f(0.f, 0.f, 0.f), glm::normalize(camera_point));


    // Adjust ray in Camera Space
//    float lensRadius = 5.0f;
//    float focalDistance = 30.0f;

    Point2f pLens = Point2f(lensRadius * WarpFunctions::squareToDiskConcentric(sampler->Get2D()));
    float ft = focalDistance / fabs(r.direction.z);

    Point3f pFocus = r.origin + ft * r.direction;

    r.origin = Point3f(pLens, 0.f);
    r.direction = glm::normalize(pFocus - r.origin);



    // Finally, transform this ray from Camera Space to World space
    glm::mat4 cameraToWorld = glm::lookAt(eye, ref, up);

    r.origin = glm::vec3(glm::inverse(cameraToWorld) * glm::vec4(r.origin, 1.0f));
    r.direction = glm::normalize(glm::vec3(glm::inverse(cameraToWorld) * glm::vec4(r.direction, 0.f)));


    return r;
}



void  ThinLenPerspectiveCamera::setLenPara(const float &lenR){

    lensRadius = lenR;
    focalDistance = glm::length(ref - eye);

}

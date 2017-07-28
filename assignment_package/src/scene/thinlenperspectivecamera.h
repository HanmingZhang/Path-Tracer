#pragma once

#include "camera.h"
#include <samplers/sampler.h>

class ThinLenPerspectiveCamera : public Camera {

public:
    ThinLenPerspectiveCamera();
    ThinLenPerspectiveCamera(const Camera &c);
    virtual Ray Raycast(const Point2f &pt) const;
    void setLenPara(const float &lenR);


    std::shared_ptr<Sampler> sampler;
    float lensRadius, focalDistance;
};

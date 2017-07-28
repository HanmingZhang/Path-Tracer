#define _USE_MATH_DEFINES
#include "warpfunctions.h"
#include <math.h>

Point3f WarpFunctions::squareToDiskUniform(const Point2f &sample)
{
    //TODO
    //throw std::runtime_error("You haven't yet implemented uniform disk warping!");


    //PBRT
    //Map uniform random numbers to [-1,1]^2
    glm::vec2 uOffset = 2.0f * sample - glm::vec2(1.f,1.f);

    //handle degeneracy at the origin
    if (fabs(uOffset.x) < 0.0001 && fabs(uOffset.y) < 0.0001){
        return glm::vec3(0, 0, 0);
    }


    //Apply concentric mapping to point
    float theta, r;

    if(fabs((double)uOffset.x) - fabs((double)uOffset.y) > 0.0001) {
        r = uOffset.x;
        theta = M_PI_4 * (uOffset.y / uOffset.x);
    }
    else{
        r = uOffset.y;
        theta = M_PI_2 - M_PI_4 * (uOffset.x / uOffset.y);
    }

    return Point3f(r * cos(theta), r * sin(theta), 0);
}

Point3f WarpFunctions::squareToDiskConcentric(const Point2f &sample)
{
    //TODO
    //throw std::runtime_error("You haven't yet implemented concentric disk warping!");

    float r = sample.x;
    float theta = 2.0f * M_PI * sample.y;

    return Point3f(r * cos(theta), r * sin(theta), 0);
}

float WarpFunctions::squareToDiskPDF(const Point3f &sample)
{
    //TODO
    return InvPi;
}

Point3f WarpFunctions::squareToSphereUniform(const Point2f &sample)
{
    //TODO
    //throw std::runtime_error("You haven't yet implemented uniform sphere warping!");
    float z = 1.0f - 2 * sample.x;

    float x = cos(2* M_PI * sample.y) * sqrt(1.0 - z * z);
    float y = sin(2* M_PI * sample.y) * sqrt(1.0 - z * z);

    return Point3f(x, y, z);
}

float WarpFunctions::squareToSphereUniformPDF(const Point3f &sample)
{
    //TODO
    return Inv4Pi;
}

Point3f WarpFunctions::squareToSphereCapUniform(const Point2f &sample, float thetaMin)
{
    //TODO
    //throw std::runtime_error("You haven't yet implemented sphere cap warping!");
    double c = glm::cos(glm::radians(180 - thetaMin));

    float z = 1.0f - (1.0f - c) * sample.x;


    float r = sqrt(std::max(0.0f, 1.0f - z * z));


    float x = cos(2* M_PI * sample.y) * r;
    float y = sin(2* M_PI * sample.y) * r;

    return Point3f(x, y, z);
}

float WarpFunctions::squareToSphereCapUniformPDF(const Point3f &sample, float thetaMin)
{
    //TODO
    double c = glm::cos(glm::radians(180 - thetaMin));

    return ( 1.0 / (1.0 - c)) * Inv2Pi;
}

Point3f WarpFunctions::squareToHemisphereUniform(const Point2f &sample)
{
    //TODO
    //throw std::runtime_error("You haven't yet implemented uniform hemisphere warping!");
    float z = sample.x;
    float r = sqrt(std::max(0.0f, 1.0f - z * z));


    float x = cos(2* M_PI * sample.y) * r;
    float y = sin(2* M_PI * sample.y) * r;

    return Point3f(x, y, z);

}

float WarpFunctions::squareToHemisphereUniformPDF(const Point3f &sample)
{
    //TODO
    return Inv2Pi;
}

Point3f WarpFunctions::squareToHemisphereCosine(const Point2f &sample)
{
    //TODO
    //throw std::runtime_error("You haven't yet implemented cosine-weighted hemisphere warping!");
    glm::vec2 d = glm::vec2(squareToDiskUniform(sample));

    float z = sqrt(std::max(0.0, 1.0 - d.x * d.x - d.y * d.y));

    return Point3f(d.x, d.y, z);
}

float WarpFunctions::squareToHemisphereCosinePDF(const Point3f &sample)
{
    //TODO
    return sample[2] * InvPi;
}

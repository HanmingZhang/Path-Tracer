#include "volumeintegrator.h"
//#include "grid3d.h"


Color3f VolumeIntegrator::Li(const Ray &ray, const Scene &scene, std::shared_ptr<Sampler> sampler, int depth) const
{

    // Ray Marching
    Color3f C(0.f);
    Point3f xi = ray.origin;
    Vector3f dir = ray.direction;
    float T = 1.f;


    float ds = 0.01f;
    float k = 0.1f;
    bool hasEnteredVolume = false;

    float deltaT;
    Color3f c, Q;

    while(T > 0.001f){
        xi += ds * dir;

        bool isPtOutVol = scene.DensityGrid.isOutOfGrid(xi);

        if(hasEnteredVolume && isPtOutVol){
            break;
        }

        if(isPtOutVol){
            continue;
        }
        else{
            if(!hasEnteredVolume){
                hasEnteredVolume = true;
            }
        }

        deltaT = glm::exp((-k) * ds * scene.DensityGrid.GetTriLinearInteporateValue(xi));


        T *= deltaT;

        c = Color3f(scene.ColorGridR.GetTriLinearInteporateValue(xi),
                    scene.ColorGridG.GetTriLinearInteporateValue(xi),
                    scene.ColorGridB.GetTriLinearInteporateValue(xi));
        Q = Color3f(scene.LightTrGridR.GetTriLinearInteporateValue(xi),
                    scene.LightTrGridG.GetTriLinearInteporateValue(xi),
                    scene.LightTrGridB.GetTriLinearInteporateValue(xi));

        C += ((1.f - deltaT) / k)  * T * c * Q;
    }


    if(std::isnan(C[0]) ||
       std::isnan(C[1]) ||
       std::isnan(C[2])){
        std::cout << "C is nan" << std::endl;
    }

    return (1.f / 1000.f) * C ;
}

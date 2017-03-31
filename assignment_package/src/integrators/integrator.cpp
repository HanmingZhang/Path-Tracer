#include "integrator.h"

void Integrator::run()
{
    Render();
}

void Integrator::Render()
{
    // Compute the bounds of our sample, clamping to screen's max bounds if necessary
    // Instantiate a FilmTile to store this thread's pixel colors
    std::vector<Point2i> tilePixels = bounds.GetPoints();
    // For every pixel in the FilmTile:
    for(Point2i pixel : tilePixels)
    {
        ///Uncomment this to debug a particular pixel within this tile
//        if(pixel.x != 200 && pixel.y != 200)
//        {
//            continue;
//        }
        Color3f pixelColor(0.f);
        // Ask our sampler for a collection of stratified samples, then raycast through each sample
        std::vector<Point2f> pixelSamples = sampler->GenerateStratifiedSamples();
        for(Point2f sample : pixelSamples)
        {
            sample = sample + Point2f(pixel); // _sample_ is [0, 1), but it needs to be translated to the pixel's origin.
            // Generate a ray from this pixel sample
            Ray ray = camera->Raycast(sample);
            // Get the L (energy) for the ray by calling Li(ray, scene, tileSampler, arena)
            // Li is implemented by Integrator subclasses, like DirectLightingIntegrator
            Color3f L = Li(ray, *scene, sampler, recursionLimit);
            // Accumulate color in the pixel
            pixelColor += L;
        }
        // Average all samples' energies
        pixelColor /= pixelSamples.size();
        film->SetPixelColor(pixel, glm::clamp(pixelColor, 0.f, 1.f));
    }
    //We're done here! All pixels have been given an averaged color.
}


void Integrator::ClampBounds()
{
    Point2i max = bounds.Max();
    max = Point2i(std::min(max.x, film->bounds.Max().x), std::min(max.y, film->bounds.Max().y));
    bounds = Bounds2i(bounds.Min(), max);
}



Color3f GetEnvironmentMapColor(const Ray &ray, const Scene &scene, std::shared_ptr<Sampler> sampler){


    //TODO
    Color3f L(0.0f);
    Vector3f woW = -ray.direction;


    // First, we get the intersection of this ray
    // It must have intersection
    Intersection isect;
    scene.environmentBox->Intersect(ray, &isect);



    // if intersection exists, compute Le first
    L += isect.Le(woW);


    // set up BSDF
    isect.ProduceBSDF();


    // reverse the normal
    isect.normalGeometric *= -1.0f;
    isect.bsdf->normal *= -1.0f;



    Vector3f wiW(0.0f);
    float pdf = 0.0f;


    // we randomly choose a light source
    float lightChooseRandomNum = sampler->Get1D();
    int lightNum = scene.lights.size();
    int lightIdx = std::min((int)std::floor(lightChooseRandomNum * (float)lightNum),
                            lightNum - 1);


    // wiW will change here !
    Color3f Li = scene.lights[lightIdx]->Sample_Li(isect, sampler->Get2D(), &wiW, &pdf);


    // we normalize wiW here!
    // before we do normalize, we must check wiW is not (0, 0, 0)
    if(wiW.x != 0.f ||
       wiW.y != 0.f ||
       wiW.z != 0.f){
        wiW = glm::normalize(wiW);
    }


    Color3f f_term = isect.bsdf->f(woW, wiW);


    // divide the PDF obtain from Sample_Li by the number of light sources in your scene
    pdf = pdf / (float)lightNum;



    // directly use LTE here!
    L += (f_term * Li * AbsDot(wiW,glm::normalize(isect.normalGeometric))) / pdf;


    return L;


}

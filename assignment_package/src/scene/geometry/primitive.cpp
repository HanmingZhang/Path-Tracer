#include "primitive.h"

bool Primitive::Intersect(const Ray &r, Intersection *isect) const
{
    if(!shape->Intersect(r, isect)) return false;
    isect->objectHit = this;
    // We create a BSDF for this intersection in our Integrator classes
    return true;
}

bool Primitive::ProduceBSDF(Intersection *isect) const
{
    if(material)
    {
        material->ProduceBSDF(isect);
        return true;
    }
    return false;
}

const Light* Primitive::GetLight() const
{
    return Light.get();
}


const Material* Primitive::GetMaterial() const
{
    return material.get();
}

Bounds3f Primitive::WorldBound() const{

    if(shape != nullptr){
        return shape->WorldBound();
    }

    else qDebug() << "the shape ptr of this Primitive is null!";
}

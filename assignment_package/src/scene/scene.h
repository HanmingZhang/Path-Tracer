#pragma once
#include <QList>
#include <raytracing/film.h>
#include <scene/camera.h>
#include <scene/thinlenperspectivecamera.h>
#include <scene/lights/light.h>
#include <scene/geometry/shape.h>
#include <photonmapper.h>
#include "bvh.h"

#include "kdtree.h"
#include "grid3d.h"


class Primitive;
class BVHAccel;
class Material;
class Light;

class Scene
{
public:
    Scene();
    ~Scene();
    QList<std::shared_ptr<Primitive>> primitives;
    QList<std::shared_ptr<Material>> materials;
    QList<std::shared_ptr<Light>> lights;

    std::shared_ptr<Primitive> environmentBox;



    Camera camera;
    ThinLenPerspectiveCamera lenCamera;

    bool UseThinLenCam;


    Film film;

    BVHAccel* bvh;

    KdTreeAccel* kdtree;

    QList<std::shared_ptr<Drawable>> drawables;

    bool isPhotonMapperExist;
    PhotonMapper globalPhotonMap, causticPhotonMap;

    int AARate;

    void SetCamera(Camera c);

    void CreateTestScene();
    void Clear();

    bool Intersect(const Ray& ray, Intersection* isect) const;

    void clearBVH();
    void clearKdtree();

    void EmitPhotonsFromLights();


    Grid3D DensityGrid;
    Grid3D ColorGridR;
    Grid3D ColorGridG;
    Grid3D ColorGridB;
    Grid3D LightTrGridR;
    Grid3D LightTrGridG;
    Grid3D LightTrGridB;

    void InitializeVolGrids();
    void PreCompLightTrGrid(float k, float ds);


};

#include <scene/scene.h>

#include <scene/geometry/cube.h>
#include <scene/geometry/sphere.h>
#include <scene/geometry/mesh.h>
#include <scene/geometry/squareplane.h>
#include <scene/materials/mattematerial.h>
#include <scene/lights/diffusearealight.h>

#include <grid3d.h>

#include <QTime>

Scene::Scene() : UseThinLenCam(false),bvh(nullptr), kdtree(nullptr), isPhotonMapperExist(false), AARate(1)
{}

Scene::~Scene()
{
    delete bvh;
    delete kdtree;
}

void Scene::SetCamera(Camera c)
{
    c.width *= AARate;
    c.height *= AARate;

    camera = Camera(c);

    lenCamera = ThinLenPerspectiveCamera(c);

    camera.create();

    film.SetDimensions(c.width, c.height);
}

bool Scene::Intersect(const Ray &ray, Intersection *isect) const
{
    if(bvh)
    {
        return bvh->Intersect(ray, isect);
    }


    else if(kdtree)
    {
        return kdtree->Intersect(ray, isect);
    }


    else
    {
        bool result = false;
        for(std::shared_ptr<Primitive> p : primitives)
        {
            Intersection testIsect;
            if(p->Intersect(ray, &testIsect))
            {
                if(testIsect.t < isect->t || isect->t < 0)
                {
                    *isect = testIsect;
                    result = true;
                }
            }
        }
        return result;
    }
    return false;
}

void Scene::CreateTestScene()
{
    //Floor
    //Area light
    //Figure in front of light

    auto matteWhite = std::make_shared<MatteMaterial>(Color3f(1,1,1), 0, nullptr, nullptr);
    auto matteRed = std::make_shared<MatteMaterial>(Color3f(1,0,0), 0, nullptr, nullptr);
    auto matteGreen = std::make_shared<MatteMaterial>(Color3f(0,1,0), 0, nullptr, nullptr);

    // Floor, which is a large white plane
    auto floor = std::make_shared<SquarePlane>();
    floor->transform = Transform(Vector3f(0,0,0), Vector3f(-90,0,0), Vector3f(10,10,1));
    auto floorPrim = std::make_shared<Primitive>(floor);
    floorPrim->material = matteWhite;
    floorPrim->name = QString("Floor");

    // Light source, which is a diffuse area light with a large plane as its shape
    auto lightSquare = std::make_shared<SquarePlane>();
    lightSquare->transform = Transform(Vector3f(0,2.5f,5), Vector3f(0,180,0), Vector3f(8, 5, 1));
    auto lightSource = std::make_shared<DiffuseAreaLight>(lightSquare->transform, Color3f(1,1,1) * 2.f, lightSquare);
    auto lightPrim = std::make_shared<Primitive>(lightSquare, nullptr, lightSource);
    lightPrim->name = QString("Light Source");
    lightSource->name = QString("Light Source 1");

    // Light source 2, which is a diffuse area light with a large plane as its shape
    auto lightSquare2 = std::make_shared<SquarePlane>();
    lightSquare2->transform = Transform(Vector3f(5,2.5f,0), Vector3f(0,90,0), Vector3f(8, 5, 1));
    auto lightSource2 = std::make_shared<DiffuseAreaLight>(lightSquare2->transform, Color3f(0.9,1,0.7) * 2.f, lightSquare2, true);
    auto lightPrim2 = std::make_shared<Primitive>(lightSquare2, nullptr, lightSource2);
    lightPrim2->name = QString("Light Source 2");
    lightSource2->name = QString("Light Source 2");

    // Shadow casting shape, which is a red sphere
    auto sphere = std::make_shared<Sphere>();
    sphere->transform = Transform(Vector3f(0,1,0), Vector3f(0,0,0), Vector3f(1,1,1));
    auto spherePrim = std::make_shared<Primitive>(sphere);
    spherePrim->material = matteRed;
    spherePrim->name = QString("Red Sphere");

    // Back wall, which is a green rectangle
    auto greenWall = std::make_shared<SquarePlane>();
    greenWall->transform = Transform(Vector3f(-5,2.5f,0), Vector3f(0,90,0), Vector3f(10, 5, 1));
    auto greenWallPrim = std::make_shared<Primitive>(greenWall);
    greenWallPrim->material = matteGreen;
    greenWallPrim->name = QString("Wall");


    primitives.append(floorPrim);
    primitives.append(lightPrim);
    primitives.append(lightPrim2);
    primitives.append(spherePrim);
    primitives.append(greenWallPrim);

    lights.append(lightSource);
    lights.append(lightSource2);

    for(std::shared_ptr<Primitive> p : primitives)
    {
        p->shape->create();
    }

    camera = Camera(400, 400, Point3f(5, 8, -5), Point3f(0,0,0), Vector3f(0,1,0));
    camera.near_clip = 0.1f;
    camera.far_clip = 100.0f;
    camera.create();
    film = Film(400, 400);
}

void Scene::Clear()
{
    // These lists contain shared_ptrs
    // so the pointers will be freed
    // if appropriate when we clear the lists.
    primitives.clear();
    lights.clear();
    materials.clear();
    drawables.clear();
    camera = Camera();

    lenCamera = ThinLenPerspectiveCamera();

    film = Film();
}

void Scene::clearBVH()
{
    if(bvh)
    {
        delete bvh;
        bvh = nullptr;
    }
}

void Scene::clearKdtree()
{
    if(kdtree)
    {
        delete kdtree;
        kdtree = nullptr;
    }
}

void Scene::EmitPhotonsFromLights(){

    if(lights.size() == 0){
        return;
    }

    QTime myTimer;
    std::cout << "Emit photons begins!" << std::endl;
    myTimer.start();


    int numberOfPhotons = 800000;
    int photonBounceTime = 5;
    float lightPowerCofficient = 13.0f;


    Sampler sampler = Sampler(numberOfPhotons, 0);


    // LightSample is [0, 1)
    std::vector<Point2f> lightSamples = sampler.GenerateUniformSamples();


    // We select a light where photon shoot from base on
    // this light's emittedColor / Power
    std::vector<float> lightPowerRatio;
    lightPowerRatio.reserve(lights.size());
    float sumOfAllPower = 0.f;
    float sumOfAllRatio = 0.f;

    for(int i = 0; i < lights.size(); i++){
        Color3f each_power = lights[i]->GetInitialPhotonPower();
        lightPowerRatio[i] = each_power.x + each_power.y + each_power.z;
        sumOfAllPower += each_power.x + each_power.y + each_power.z;
    }

    for(int i = 0; i < lights.size(); i++){
        lightPowerRatio[i] = (lightPowerRatio[i] / sumOfAllPower) + sumOfAllRatio;
        sumOfAllRatio = lightPowerRatio[i];
    }




    // Shoot photon begains
    for(int i = 0; i < lightSamples.size(); i++){

            //select a light
            float lightChoice = sampler.Get1D();
            int lightIdx = 0;

            for(int i = 0; i < lights.size(); i++){
                 float largerThan = (i == 0) ? 0.f : lightPowerRatio[i - 1];
                 float smallerThan = lightPowerRatio[i];
                 if(lightChoice >= largerThan && lightChoice < smallerThan){
                     lightIdx = i;
                     break;
                 }
            }

            Light *light_ptr = lights[lightIdx].get();



            // Initialize a photon
            // we use uniform light area sampling to compute this light source's PDF
            float light_pdf = 0.f;
            Ray photonBounceRay = light_ptr->generatePhotonRay(lightSamples[i], &sampler, &light_pdf);
            Color3f power = lightPowerCofficient * light_ptr->GetInitialPhotonPower() / ((float) numberOfPhotons);
            Color3f photon_weight(1.0f);

            bool isPreviousBounceSpecular = false;
            bool isHitDiffuse = false;

            // start bounce a photon
            for(int j = 0; j < photonBounceTime; j++){

                Intersection isect;
                Intersect(photonBounceRay, &isect);

                // if this photon hit nothing,
                // no need to do rest bounce
                if(isect.objectHit == nullptr){
                    break;
                }


                // if this photon hit something has no material
                // then it won't bounce
                if(isect.objectHit->material == nullptr){
                    break;
                }


                isect.ProduceBSDF();

                Vector3f woW = -photonBounceRay.direction;
                Vector3f wiW(0.f);
                float pdf = 0.f;
                BxDFType samplerType;


                // wiW and pdf Will change here!
                Color3f f_term = isect.bsdf->Sample_f(woW, &wiW, sampler.Get2D(), &pdf, BSDF_ALL, &samplerType);
                photon_weight *= f_term * AbsDot(wiW,glm::normalize(isect.normalGeometric)) / pdf;

                // check pdf
                if(pdf == 0.f){
                    break;
                }

                // we normalize wiW here!
                // before we do normalization, we must check wiW is not (0, 0, 0)
                if(wiW.x != 0.f ||
                   wiW.y != 0.f ||
                   wiW.z != 0.f){
                    wiW = glm::normalize(wiW);
                }


                // if this is the first bounce,
                // we don't record this photon
                if(j == 0){
                    photonBounceRay = isect.SpawnRay(wiW);
                    //if((samplerType & BSDF_SPECULAR) == BSDF_SPECULAR){isPreviousBounceSpecular = true;}
                    if(isect.bsdf->BxDFsHaveSpecularFlags()){isPreviousBounceSpecular = true;}
                    else{
                        isPreviousBounceSpecular = false;
                        isHitDiffuse = true;
                    }
                    continue;
                }


                // Conduct Russian Roulette
                float max_component = std::max(photon_weight.x, std::max(photon_weight.y, photon_weight.z));
                if(max_component < sampler.Get1D()){
                    break;
                }
                photon_weight /= max_component;


                // Let's say whether should record photon
                // at this bounce
                // Check whether this surface is specular or not
                // photons will only stay on non-specular surfaces

                //if(!((samplerType & BSDF_SPECULAR) == BSDF_SPECULAR)){
                if(!isect.bsdf->BxDFsHaveSpecularFlags()){
                    if(isPreviousBounceSpecular == true &&
                                   isHitDiffuse == false){
                        causticPhotonMap.AddPhoton(isect.point,
                                                   glm::normalize(photonBounceRay.direction),
                                                   power * photon_weight);
                    }
                    else{
                        globalPhotonMap.AddPhoton(isect.point,
                                                  glm::normalize(photonBounceRay.direction),
                                                  power * photon_weight);
                    }
                }


                // initialize a new bounce
                photonBounceRay = isect.SpawnRay(wiW);
                //if((samplerType & BSDF_SPECULAR) == BSDF_SPECULAR){isPreviousBounceSpecular = true;}
                if(isect.bsdf->BxDFsHaveSpecularFlags()){isPreviousBounceSpecular = true;}
                else{
                    isPreviousBounceSpecular = false;
                    isHitDiffuse = true;
                }
            }
    }


    causticPhotonMap.PrepareForIrradianceEstimation();

    globalPhotonMap.PrepareForIrradianceEstimation();

    std::cout << "Emit photons ends, and total time used: "
              << myTimer.elapsed() << " milliseconds." << std::endl;

}



void Scene::InitializeVolGrids(){
//    DensityGrid = Grid3D(100, 100, 350, Point3f(5.f,-2.5f, -30.f), Point3f(-5.f, 7.5f, 5.f));

//    ColorGridR = Grid3D(100, 100, 350, Point3f(5.f,-2.5f, -30.f), Point3f(-5.f, 7.5f, 5.f));
//    ColorGridG = Grid3D(100, 100, 350, Point3f(5.f,-2.5f, -30.f), Point3f(-5.f, 7.5f, 5.f));
//    ColorGridB = Grid3D(100, 100, 350, Point3f(5.f,-2.5f, -30.f), Point3f(-5.f, 7.5f, 5.f));

//    LightTrGridR = Grid3D(100, 100, 350, Point3f(5.f,-2.5f, -30.f), Point3f(-5.f, 7.5f, 5.f), 0.f);
//    LightTrGridG = Grid3D(100, 100, 350, Point3f(5.f,-2.5f, -30.f), Point3f(-5.f, 7.5f, 5.f), 0.f);
//    LightTrGridB = Grid3D(100, 100, 350, Point3f(5.f,-2.5f, -30.f), Point3f(-5.f, 7.5f, 5.f), 0.f);

    int n = 100;
    Point3f start = Point3f(5.f, 7.5f, 5.f);
    Point3f end = Point3f(-5.f,-2.5f, -30.f);

    DensityGrid = Grid3D(n, n, 350, end, start, 2.0f);

    ColorGridR = Grid3D(n, n, 350, end, start, 0.6f);
    ColorGridG = Grid3D(n, n, 350, end, start, 0.6f);
    ColorGridB = Grid3D(n, n, 350, end, start, 0.6f);

    LightTrGridR = Grid3D(n, n, 350, end, start, 0.f);
    LightTrGridG = Grid3D(n, n, 350, end, start, 0.f);
    LightTrGridB = Grid3D(n, n, 350, end, start, 0.f);
}

void Scene::PreCompLightTrGrid(float k, float ds){

    if(lights.size() == 0){
        return;
    }

    QTime myTimer;
    std::cout << "Shooting light beams(Vol Integrator) begins!" << std::endl;
    myTimer.start();

    int numberOfBeams = 250000;
    int maxBounceTime = 5;
    float endThreshold = 0.001f;

    Sampler sampler = Sampler(numberOfBeams, 0);

    // LightSample is [0, 1)
    std::vector<Point2f> lightSamples = sampler.GenerateStratifiedSamples();


    // We select a light where photon shoot from base on
    // this light's emittedColor / Power
    std::vector<float> lightPowerRatio;
    lightPowerRatio.reserve(lights.size());
    float sumOfAllPower = 0.f;
    float sumOfAllRatio = 0.f;

    for(int i = 0; i < lights.size(); i++){
        Color3f each_power = lights[i]->GetInitialPhotonPower();
        lightPowerRatio[i] = each_power.x + each_power.y + each_power.z;
        sumOfAllPower += each_power.x + each_power.y + each_power.z;
    }

    for(int i = 0; i < lights.size(); i++){
        lightPowerRatio[i] = (lightPowerRatio[i] / sumOfAllPower) + sumOfAllRatio;
        sumOfAllRatio = lightPowerRatio[i];
    }



    // Shoot light beam begains
    for(int i = 0; i < lightSamples.size(); i++){

            //select a light
            float lightChoice = sampler.Get1D();
            int lightIdx = 0;

            for(int i = 0; i < lights.size(); i++){
                 float largerThan = (i == 0) ? 0.f : lightPowerRatio[i - 1];
                 float smallerThan = lightPowerRatio[i];
                 if(lightChoice >= largerThan && lightChoice < smallerThan){
                     lightIdx = i;
                     break;
                 }
            }

            Light *light_ptr = lights[lightIdx].get();


            // Initialize a light beam
            // we use uniform light area sampling to compute this light source's PDF
            float light_pdf = 0.f;
            Ray BounceRay = light_ptr->generatePhotonRay(lightSamples[i], &sampler, &light_pdf);

            Point3f Pos = BounceRay.origin;
            Point3f NewPos;
            Vector3f Dir= BounceRay.direction;
            Vector3f NewDir;
            int BounceTime = 0;
            float Q = 1.f;

            Intersection isect;
            Intersect(BounceRay, &isect);
            float t = 0.f;
            float t_max = Infinity;

            if(!(isect.objectHit == nullptr || isect.objectHit->material == nullptr)){
                    isect.ProduceBSDF();
                    Vector3f woW = -BounceRay.direction;
                    Vector3f wiW(0.f);
                    float pdf = 0.f;
                    BxDFType samplerType;

                    // wiW and pdf Will change here!
                    Color3f f_term = isect.bsdf->Sample_f(woW, &wiW, sampler.Get2D(), &pdf, BSDF_ALL, &samplerType);

                    BounceRay = isect.SpawnRay(wiW);

                    NewPos = BounceRay.origin;
                    NewDir = BounceRay.direction;
                    t_max = glm::distance(NewPos, Pos);
            }


            Color3f initialPower = light_ptr->GetInitialPhotonPower();
            LightTrGridR.SetValue(Pos, LightTrGridR.GetValue(Pos) + initialPower[0]);
            LightTrGridG.SetValue(Pos, LightTrGridG.GetValue(Pos) + initialPower[1]);
            LightTrGridB.SetValue(Pos, LightTrGridB.GetValue(Pos) + initialPower[2]);


            // Ray-marching start!
            while(BounceTime < maxBounceTime && Q > endThreshold){

                t += ds;
                Pos += ds * Dir;

                if(DensityGrid.isOutOfGrid(Pos)){
                    break;
                }

                if(t >= t_max){
                    BounceTime++;
                    Pos = NewPos;
                    Dir = NewDir;

                    if(DensityGrid.isOutOfGrid(Pos)){
                        break;
                    }


                    Ray isectTestRay = Ray(Pos, Dir);
                    Intersect(isectTestRay, &isect);
                    if(!(isect.objectHit == nullptr || isect.objectHit->material == nullptr)){
                            isect.ProduceBSDF();
                            Vector3f woW = -isectTestRay.direction;
                            Vector3f wiW(0.f);
                            float pdf = 0.f;
                            BxDFType samplerType;

                            // wiW and pdf Will change here!
                            Color3f f_term = isect.bsdf->Sample_f(woW, &wiW, sampler.Get2D(), &pdf, BSDF_ALL, &samplerType);

                            BounceRay = isect.SpawnRay(wiW);

                            NewPos = BounceRay.origin;
                            NewDir = BounceRay.direction;
                            t_max = glm::distance(NewPos, Pos);
                            t = 0.f;
                    }
                    else if((isect.objectHit != nullptr) && (isect.objectHit->name == QString("Light Source"))){
                        break;
                    }
                    else{
                        t_max = Infinity;
                        t = 0.f;
                    }
                }


                Q *= glm::exp(-k * DensityGrid.GetValue(Pos) * ds);
                LightTrGridR.SetValue(Pos, LightTrGridR.GetValue(Pos) + Q * initialPower[0]);
                LightTrGridG.SetValue(Pos, LightTrGridG.GetValue(Pos) + Q * initialPower[1]);
                LightTrGridB.SetValue(Pos, LightTrGridB.GetValue(Pos) + Q * initialPower[2]);
            }

    }




    float divNum = (float)numberOfBeams * 2.0f;

    LightTrGridR.GridDivideBy(divNum);
    LightTrGridG.GridDivideBy(divNum);
    LightTrGridB.GridDivideBy(divNum);


    std::cout << "Shooting light beams(Vol Integrator) ends, and total time used: "
              << myTimer.elapsed() << " milliseconds." << std::endl;
}

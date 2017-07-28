#ifndef GRID3D_H
#define GRID3D_H

#include <vector>
#include "globals.h"
#include <cassert>

using namespace std;


class Grid3D
{
public:
    Grid3D() : i(10), j(10), k(10), LLPt(Point3f(0.f)), RHPt(Point3f(0.f)),
               lengthUnit(1.f), heightUnit(1.f), depthUnit(1.f){}
    Grid3D(int i, int j, int k, Point3f LeftLowestPt, Point3f RightHighestPt, float defaultValue = 1.f);
    void SetValue(Point3f p, float value);
    float GetValue(Point3f p);
    float GetValue(Point3f p) const;
    float GetTriLinearInteporateValue(Point3f p) const;
    bool isOutOfGrid(Point3f p) const;
    void GridDivideBy(float x);

    int i;
    int j;
    int k;
    Point3f LLPt;
    Point3f RHPt;
    float lengthUnit;
    float heightUnit;
    float depthUnit;
    vector<vector<vector <float>>> grid;
};

#endif // GRID3D_H

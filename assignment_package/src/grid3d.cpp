#include "grid3d.h"


Grid3D::Grid3D(int i, int j, int k, Point3f LeftLowestPt, Point3f RightHighestPt, float defaultValue)
    : i(i), j(j), k(k), LLPt(LeftLowestPt), RHPt(RightHighestPt)

{
    for(int m = 0; m < i; m++){
        vector<vector<float>> temp1;
        for(int n = 0; n < j; n++){
            vector<float> temp2;
            for(int p = 0; p < k; p++){
                   temp2.push_back(defaultValue);
            }
            temp1.push_back(temp2);
        }
        grid.push_back(temp1);
    }

    lengthUnit = fabs(RHPt[0] - LLPt[0]) / (float)i;
    heightUnit = fabs(RHPt[1] - LLPt[1]) / (float)j;
    depthUnit  = fabs(RHPt[2] - LLPt[2]) / (float)k;

}

void Grid3D::SetValue(Point3f p, float value){

//    assert(p[0] >= LLPt[0] && p[0] <= RHPt[0]
//        && p[1] >= LLPt[1] && p[1] <= RHPt[1]
//        && p[2] >= LLPt[2] && p[0] <= RHPt[2]);

    if(!(p[0] >= LLPt[0] && p[0] <= RHPt[0]
        && p[1] >= LLPt[1] && p[1] <= RHPt[1]
        && p[2] >= LLPt[2] && p[2] <= RHPt[2])){
        std::cout << "Set Value error!" << std::endl;
    }

//    int idxI = (p[0] - LLPt[0] + 0.5f * lengthUnit) / lengthUnit;
//    int idxJ = (p[1] - LLPt[1] + 0.5f * heightUnit) / heightUnit;
//    int idxK = (p[2] - LLPt[2] + 0.5f * depthUnit) / depthUnit;

    int idxI = (p[0] - LLPt[0]) / lengthUnit;
    int idxJ = (p[1] - LLPt[1]) / heightUnit;
    int idxK = (p[2] - LLPt[2]) / depthUnit;

    if(idxI == i){
        idxI--;
    }
    if(idxJ == j){
        idxJ--;
    }
    if(idxK == k){
        idxK--;
    }

    grid[idxI][idxJ][idxK] = value;
}

float Grid3D::GetValue(Point3f p){
//    assert(p[0] >= LLPt[0] && p[0] <= RHPt[0]
//        && p[1] >= LLPt[1] && p[1] <= RHPt[1]
//        && p[2] >= LLPt[2] && p[0] <= RHPt[2]);


    if(!(p[0] >= LLPt[0] && p[0] <= RHPt[0]
        && p[1] >= LLPt[1] && p[1] <= RHPt[1]
        && p[2] >= LLPt[2] && p[2] <= RHPt[2])){
        std::cout << "Get Value error!" << std::endl;
    }

//    int idxI = (p[0] - LLPt[0] + 0.5f * lengthUnit) / lengthUnit;
//    int idxJ = (p[1] - LLPt[1] + 0.5f * heightUnit) / heightUnit;
//    int idxK = (p[2] - LLPt[2] + 0.5f * depthUnit) / depthUnit;

    int idxI = (p[0] - LLPt[0]) / lengthUnit;
    int idxJ = (p[1] - LLPt[1]) / heightUnit;
    int idxK = (p[2] - LLPt[2]) / depthUnit;

    if(idxI == i){
        idxI--;
    }
    if(idxJ == j){
        idxJ--;
    }
    if(idxK == k){
        idxK--;
    }


    return grid[idxI][idxJ][idxK];
}


float Grid3D::GetValue(Point3f p) const{
    if(!(p[0] >= LLPt[0] && p[0] <= RHPt[0]
        && p[1] >= LLPt[1] && p[1] <= RHPt[1]
        && p[2] >= LLPt[2] && p[2] <= RHPt[2])){
        std::cout << "Get Value error!" << std::endl;
    }


    int idxI = (p[0] - LLPt[0]) / lengthUnit;
    int idxJ = (p[1] - LLPt[1]) / heightUnit;
    int idxK = (p[2] - LLPt[2]) / depthUnit;

    if(idxI == i){
        idxI--;
    }
    if(idxJ == j){
        idxJ--;
    }
    if(idxK == k){
        idxK--;
    }


    return grid[idxI][idxJ][idxK];

}

float Grid3D::GetTriLinearInteporateValue(Point3f p) const{
//    assert(p[0] >= LLPt[0] && p[0] <= RHPt[0]
//        && p[1] >= LLPt[1] && p[1] <= RHPt[1]
//        && p[2] >= LLPt[2] && p[0] <= RHPt[2]);


    if(!(p[0] >= LLPt[0] && p[0] <= RHPt[0]
        && p[1] >= LLPt[1] && p[1] <= RHPt[1]
        && p[2] >= LLPt[2] && p[2] <= RHPt[2])){
        std::cout << "Get Tri inteporate Value error!" << std::endl;
    }


//    int idxI = (p[0] - LLPt[0] + 0.5f * lengthUnit) / lengthUnit;
//    int idxJ = (p[1] - LLPt[1] + 0.5f * heightUnit) / heightUnit;
//    int idxK = (p[2] - LLPt[2] + 0.5f * depthUnit) / depthUnit;




    int idxI = (p[0] - LLPt[0]) / lengthUnit;
    int idxJ = (p[1] - LLPt[1]) / heightUnit;
    int idxK = (p[2] - LLPt[2]) / depthUnit;


    if(idxI == i){
        idxI--;
    }
    if(idxJ == j){
        idxJ--;
    }
    if(idxK == k){
        idxK--;
    }

    int idxIAdd1 = (idxI == i - 1) ? idxI : (idxI + 1);
    int idxJAdd1 = (idxJ == j - 1) ? idxJ : (idxJ + 1);
    int idxKAdd1 = (idxK == k - 1) ? idxK : (idxK + 1);

    float wx = (lengthUnit - fabs(p[0] - (double)idxI * lengthUnit - LLPt[0])) / lengthUnit;
    float wy = (heightUnit - fabs(p[1] - (double)idxJ * heightUnit - LLPt[1])) / heightUnit;
    float wz = (depthUnit - fabs(p[2] - (double)idxK * depthUnit - LLPt[2])) / depthUnit;


    return grid[idxI][idxJ][idxK] * wx * wy * wz
         + grid[idxI][idxJ][idxKAdd1] * wx * wy * (1.f - wz)
         + grid[idxI][idxJAdd1][idxK] * wx * (1.f - wy) * wz
         + grid[idxI][idxJAdd1][idxKAdd1] * wx * (1.f - wy) * (1.f - wz)
         + grid[idxIAdd1][idxJ][idxK] * (1.f - wx) * wy * wz
         + grid[idxIAdd1][idxJ][idxKAdd1] * (1.f - wx) * wy * (1.f - wz)
         + grid[idxIAdd1][idxJAdd1][idxK] * (1.f - wx) * (1.f - wy) * wz
         + grid[idxIAdd1][idxJAdd1][idxKAdd1] * (1.f - wx) * (1.f - wy) * (1.f - wz);

}


bool Grid3D::isOutOfGrid(Point3f p) const{
    return !(p[0] >= LLPt[0] && p[0] <= RHPt[0]
          && p[1] >= LLPt[1] && p[1] <= RHPt[1]
          && p[2] >= LLPt[2] && p[2] <= RHPt[2]);
}

void Grid3D::GridDivideBy(float x){
    for(vector<vector <float>> t1 : grid){
        for(vector<float> t2 : t1){
            for(float t3 : t2){
                t3 = t3 / x;
            }
        }
    }
}

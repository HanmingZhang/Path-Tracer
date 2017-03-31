#pragma once
#include <globals.h>

class MicrofacetDistribution {
public:
    virtual ~MicrofacetDistribution();

    // Computes the differential area of microfacets on a
    // surface that are aligned with the given surface normal
    // vector wh (the half-vector between wo and its specular
    // reflection wi)
    virtual float D(const Vector3f &wh) const = 0;

    //
    virtual float Lambda(const Vector3f &w) const = 0;

    // Computes the geometric self-shadowing and interreflection term
    float G(const Vector3f &wo, const Vector3f &wi) const {
        return 1 / (1 + Lambda(wo) + Lambda(wi));
    }

    // Samples the distribution of microfacet normals to generate one
    // about which to reflect wo to create a wi.
    virtual Vector3f Sample_wh(const Vector3f &wo, const Point2f &xi) const = 0;

    // Computes the PDF of the given half-vector normal based on the
    // given incident ray direction
    virtual float Pdf(const Vector3f &wo, const Vector3f &wh) const;
};

class TrowbridgeReitzDistribution : public MicrofacetDistribution {
public:
    TrowbridgeReitzDistribution(float alphax, float alphay)
        : alphax(alphax), alphay(alphay) {}

    float D(const Vector3f &wh) const;
    Vector3f Sample_wh(const Vector3f &wo, const Point2f &xi) const;

private:
    float Lambda(const Vector3f &w) const;

    const float alphax, alphay;
};

// Converts a [0, 1] scalar measurement of surface roughness
// to the alpha term used in the microfacet distribution equation.
// It is not presently used in the materials we have provided you,
// but feel free to add it in if you want to see how it alters your renders.
inline float RoughnessToAlpha(float roughness) {
    //    return roughness;
    roughness = std::max(roughness, (float)1e-3);
    float x = std::log(roughness);
    return 1.62142f + 0.819955f * x + 0.1734f * x * x + 0.0171201f * x * x * x +
            0.000640711f * x * x * x * x;
}


// hw06 EC
// BeckmannDistribution


class BeckmannDistribution : public MicrofacetDistribution {
  public:
    // BeckmannDistribution Public Methods
    static float RoughnessToAlpha(float roughness) {
        roughness = std::max(roughness, (float)1e-3);
        float x = std::log(roughness);
        return 1.62142f + 0.819955f * x + 0.1734f * x * x +
               0.0171201f * x * x * x + 0.000640711f * x * x * x * x;
    }

    BeckmannDistribution(float alphax, float alphay, bool samplevis = true)
        :alphax(alphax), alphay(alphay) {}
    float D(const Vector3f &wh) const;
    Vector3f Sample_wh(const Vector3f &wo, const Point2f &u) const;

    //std::string ToString() const;

  private:
    // BeckmannDistribution Private Methods
    float Lambda(const Vector3f &w) const;

    // BeckmannDistribution Private Data
    const float alphax, alphay;
};

inline Vector3f SphericalDirection(float sinTheta, float cosTheta, float phi) {
    return Vector3f(sinTheta * std::cos(phi), sinTheta * std::sin(phi),
                    cosTheta);
}

//inline Vector3f SphericalDirection(Float sinTheta, Float cosTheta, Float phi,
//                                   const Vector3f &x, const Vector3f &y,
//                                   const Vector3f &z) {
//    return sinTheta * std::cos(phi) * x + sinTheta * std::sin(phi) * y +
//           cosTheta * z;
//}


// Blinn part
class Blinn : public MicrofacetDistribution {
    public:
     Blinn(float e, float alphax, float alphay): alphax(alphax), alphay(alphay)
     { if (e > 10000.f || std::isnan(e)) e = 10000.f;
                      exponent = e; }
     // Blinn Public Methods
     float D(const Vector3f &wh) const {
         float costhetah = AbsCosTheta(wh);
         return (exponent + 2.0f) * Inv2Pi * powf(costhetah, exponent);
     }

     Vector3f Sample_wh(const Vector3f &wo, const Point2f &u) const;

     // virtual void Sample_f(const Vector &wi, Vector *sampled_f, float u1, float u2, float *pdf) const;
     virtual float Pdf(const Vector3f &wi, const Vector3f &wo) const;
 private:
    float exponent;

    float Lambda(const Vector3f &w) const;
    const float alphax, alphay;
 };

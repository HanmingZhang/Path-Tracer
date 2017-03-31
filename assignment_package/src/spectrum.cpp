#include "spectrum.h"

inline float Lerp(float t, float v1, float v2) { return (1 - t) * v1 + t * v2; }

inline void XYZToRGB(const float xyz[3], float rgb[3]) {
    rgb[0] = 3.240479f * xyz[0] - 1.537150f * xyz[1] - 0.498535f * xyz[2];
    rgb[1] = -0.969256f * xyz[0] + 1.875991f * xyz[1] + 0.041556f * xyz[2];
    rgb[2] = 0.055648f * xyz[0] - 0.204043f * xyz[1] + 1.057311f * xyz[2];
}

bool SpectrumSamplesSorted(const float *lambda, const float *vals, int n) {
    for (int i = 0; i < n - 1; ++i)
        if (lambda[i] > lambda[i + 1]) return false;
    return true;
}

void SortSpectrumSamples(float *lambda, float *vals, int n) {

    std::vector<std::pair<float, float>> sortVec;
    sortVec.reserve(n);

    for (int i = 0; i < n; ++i)
        sortVec.push_back(std::make_pair(lambda[i], vals[i]));



    std::sort(sortVec.begin(), sortVec.end());
    for (int i = 0; i < n; ++i) {
        lambda[i] = sortVec[i].first;
        vals[i] = sortVec[i].second;
    }
}


float AverageSpectrumSamples(const float *lambda, const float *vals, int n, float lambdaStart, float lambdaEnd) {



    // Handle cases with out-of-bounds range or single sample only
    if (lambdaEnd <= lambda[0]) return vals[0];
    if (lambdaStart >= lambda[n - 1]) return vals[n - 1];
    if (n == 1) return vals[0];
    float sum = 0;

    // Add contributions of constant segments before/after samples
    if (lambdaStart < lambda[0]) sum += vals[0] * (lambda[0] - lambdaStart);
    if (lambdaEnd > lambda[n - 1])
        sum += vals[n - 1] * (lambdaEnd - lambda[n - 1]);

    // Advance to first relevant wavelength segment
    int i = 0;
    while (lambdaStart > lambda[i + 1]) ++i;

    // Loop over wavelength sample segments and add contributions
    auto interp = [lambda, vals](float w, int i) {
        return Lerp((w - lambda[i]) / (lambda[i + 1] - lambda[i]), vals[i],
                    vals[i + 1]);
    };
    for (; i + 1 < n && lambdaEnd >= lambda[i]; ++i) {
        float segLambdaStart = std::max(lambdaStart, lambda[i]);
        float segLambdaEnd = std::min(lambdaEnd, lambda[i + 1]);
        sum += 0.5 * (interp(segLambdaStart, i) + interp(segLambdaEnd, i)) *
               (segLambdaEnd - segLambdaStart);
    }
    return sum / (lambdaEnd - lambdaStart);
}



// Spectrum Public Methods
Spectrum::Spectrum(float v) {
    for (int i = 0; i < nSpectrumSamples; ++i) c[i] = v;

}




// Define basic operations for class Spectrum

// Add operations
Spectrum& Spectrum::operator+=(const Spectrum &s2) {
    for (int i = 0; i < nSpectrumSamples; ++i) c[i] += s2.c[i];
    return *this;
}
Spectrum Spectrum::operator+(const Spectrum &s2) const {
    Spectrum ret = *this;
    for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] += s2.c[i];
    return ret;
}

// Subtract operations
Spectrum& Spectrum::operator-=(const Spectrum &s2){
    for (int i = 0; i < nSpectrumSamples; ++i) c[i] -= s2.c[i];
    return *this;
}
Spectrum Spectrum::operator-(const Spectrum &s2) const {
    Spectrum ret = *this;
    for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] -= s2.c[i];
    return ret;
}

// Multiply operations
Spectrum Spectrum::operator*(const Spectrum &sp) const {
    Spectrum ret = *this;
    for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] *= sp.c[i];
    return ret;
}
Spectrum& Spectrum::operator*=(const Spectrum &sp) {
    for (int i = 0; i < nSpectrumSamples; ++i) c[i] *= sp.c[i];
    return *this;
}
Spectrum Spectrum::operator*(float a) const {
    Spectrum ret = *this;
    for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] *= a;
    return ret;
}
Spectrum& Spectrum::operator*=(float a) {
    for (int i = 0; i < nSpectrumSamples; ++i) c[i] *= a;
    return *this;
}


// Divide operations
Spectrum Spectrum::operator/(const Spectrum &s2) const {
    Spectrum ret = *this;
    for (int i = 0; i < nSpectrumSamples; ++i) {
      ret.c[i] /= s2.c[i];
    }
    return ret;
}
Spectrum Spectrum::operator/(float a) const {
    Spectrum ret = *this;
    for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] /= a;
    return ret;
}
Spectrum& Spectrum::operator/=(float a) {
    for (int i = 0; i < nSpectrumSamples; ++i) c[i] /= a;
    return *this;
}


// bool opeartions to check whether is equal or not
bool Spectrum::operator==(const Spectrum &sp) const {
    for (int i = 0; i < nSpectrumSamples; ++i)
        if (c[i] != sp.c[i]) return false;
    return true;
}
bool Spectrum::operator!=(const Spectrum &sp) const {
    return !(*this == sp);
}

// sometimes we need to check whether it's black or not
bool Spectrum::IsBlack() const {
    for (int i = 0; i < nSpectrumSamples; ++i)
        if (c[i] != 0.) return false;
    return true;
}


// we need [] operation to get certain wavelength value
float& Spectrum::operator[](int i) {
    // first check idx
    if(i >= 0 && i < nSpectrumSamples){
        return c[i];
    }
}

float Spectrum::operator[](int i) const {
    if(i >= 0 && i < nSpectrumSamples){
        return c[i];
    }
    else return -1.0f;
}


Spectrum Spectrum::FromSampled(const float *lambda, const float *v,int n) {
    // Sort samples if unordered, use sorted for returned spectrum
    if (!SpectrumSamplesSorted(lambda, v, n)) {
        std::vector<float> slambda(&lambda[0], &lambda[n]);
        std::vector<float> sv(&v[0], &v[n]);
        SortSpectrumSamples(&slambda[0], &sv[0], n);
        return FromSampled(&slambda[0], &sv[0], n);
    }

    Spectrum r;
    for (int i = 0; i < nSpectrumSamples; ++i) {
        // Compute average value of given SPD over $i$th sample's range
        float lambda0 = Lerp(float(i) / float(nSpectrumSamples),
                             sampledLambdaStart, sampledLambdaEnd);
        float lambda1 = Lerp(float(i + 1) / float(nSpectrumSamples),
                             sampledLambdaStart, sampledLambdaEnd);
        r.c[i] = AverageSpectrumSamples(lambda, v, n, lambda0, lambda1);
    }
    return r;
}


void Spectrum::ToXYZ(float xyz[3]) const {
    xyz[0] = xyz[1] = xyz[2] = 0.f;

    for (int i = 0; i < nSpectrumSamples; ++i) {
        xyz[0] += X.c[i] * c[i];
        xyz[1] += Y.c[i] * c[i];
        xyz[2] += Z.c[i] * c[i];
    }


    float scale = float(sampledLambdaEnd - sampledLambdaStart) / float(CIE_Y_integral * nSpectrumSamples);
    xyz[0] *= scale;
    xyz[1] *= scale;
    xyz[2] *= scale;
}


float Spectrum::y() const {
    float yy = 0.f;
    for (int i = 0; i < nSpectrumSamples; ++i) yy += Y.c[i] * c[i];

    return yy * float(sampledLambdaEnd - sampledLambdaStart) / float(CIE_Y_integral * nSpectrumSamples);
}

void Spectrum::ToRGB(float rgb[3]) const {
    float xyz[3];
    ToXYZ(xyz);
    XYZToRGB(xyz, rgb);
}

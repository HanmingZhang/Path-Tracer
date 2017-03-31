#ifndef SPECTRUM_H
#define SPECTRUM_H


#endif // SPECTRUM_H
#include <vector>
#include <algorithm>

// Spectrum Utility Declarations
static const int sampledLambdaStart = 400;
static const int sampledLambdaEnd = 700;
static const int nSpectrumSamples = 60;

bool SpectrumSamplesSorted(const float *lambda, const float *vals, int n);
void SortSpectrumSamples(float *lambda, float *vals, int n);
float AverageSpectrumSamples(const float *lambda, const float *vals, int n, float lambdaStart, float lambdaEnd);
inline float Lerp(float t, float v1, float v2);
inline void XYZToRGB(const float xyz[3], float rgb[3]);


// Spectral Data Declarations
static const int nCIESamples = 471;
static const float CIE_Y_integral = 106.856895;
static const int nRGB2SpectSamples = 32;




// Spectrum Declarations

class Spectrum {
  public:
    // Spectrum Public Methods
    Spectrum(float v = 0.f);


    // Define basic operations for class Spectrum

    // Add operations
    Spectrum &operator+=(const Spectrum &s2);
    Spectrum operator+(const Spectrum &s2) const;

    // Subtract operations
    Spectrum &operator-=(const Spectrum &s2);
    Spectrum operator-(const Spectrum &s2) const;

    // Multiply operations
    Spectrum operator*(const Spectrum &sp) const;
    Spectrum &operator*=(const Spectrum &sp);
    Spectrum operator*(float a) const;
    Spectrum &operator*=(float a);

    // Divide operations
    Spectrum operator/(const Spectrum &s2) const;
    Spectrum operator/(float a) const;
    Spectrum &operator/=(float a);


    // bool opeartions to check whether is equal or not
    bool operator==(const Spectrum &sp) const;
    bool operator!=(const Spectrum &sp) const;

    // sometimes we need to check whether it's black or not
    bool IsBlack() const;


    // we need [] operation to get certain wavelength value
    float &operator[](int i);

    float operator[](int i) const;


    // form Spectrum from sample
    Spectrum FromSampled(const float *lambda, const float *v,int n);

    void ToXYZ(float xyz[3]) const;
    float y() const;
    void ToRGB(float rgb[3]) const;


    // samples we take
    static const int nSamples = nSpectrumSamples;

  private:
    // Spectrum Protected Data
    float c[nSpectrumSamples];

};


static Spectrum X, Y, Z;

#ifndef DISTRIBUTION2D_H
#define DISTRIBUTION2D_H

#include "globals.h"



struct Distribution1D {
    // Distribution1D Public Methods
    Distribution1D(const Float *f, int n) : func(f, f + n), cdf(n + 1) {
        // Compute integral of step function at $x_i$
        cdf[0] = 0;
        for (int i = 1; i < n + 1; ++i) cdf[i] = cdf[i - 1] + func[i - 1] / n;

        // Transform step function integral into CDF
        funcInt = cdf[n];
        if (funcInt == 0) {
            for (int i = 1; i < n + 1; ++i) cdf[i] = Float(i) / Float(n);
        } else {
            for (int i = 1; i < n + 1; ++i) cdf[i] /= funcInt;
        }
    }
    int Count() const { return func.size(); }
    Float SampleContinuous(Float u, Float *pdf, int *off = nullptr) const {
        // Find surrounding CDF segments and _offset_
        int offset = FindInterval(cdf.size(),
                                  [&](int index) { return cdf[index] <= u; });
        if (off) *off = offset;
        // Compute offset along CDF segment
        Float du = u - cdf[offset];
        if ((cdf[offset + 1] - cdf[offset]) > 0) {
            //CHECK_GT(cdf[offset + 1], cdf[offset]);
            du /= (cdf[offset + 1] - cdf[offset]);
        }
        //DCHECK(!std::isnan(du));

        // Compute PDF for sampled offset
        if (pdf) *pdf = (funcInt > 0) ? func[offset] / funcInt : 0;

        // Return $x\in{}[0,1)$ corresponding to sample
        return (offset + du) / Count();
    }
    int SampleDiscrete(Float u, Float *pdf = nullptr,
                       Float *uRemapped = nullptr) const {
        // Find surrounding CDF segments and _offset_
        int offset = FindInterval(cdf.size(),
                                  [&](int index) { return cdf[index] <= u; });
        if (pdf) *pdf = (funcInt > 0) ? func[offset] / (funcInt * Count()) : 0;
        if (uRemapped)
            *uRemapped = (u - cdf[offset]) / (cdf[offset + 1] - cdf[offset]);
        //if (uRemapped) CHECK(*uRemapped >= 0.f && *uRemapped <= 1.f);
        return offset;
    }
    Float DiscretePDF(int index) const {
        //CHECK(index >= 0 && index < Count());
        return func[index] / (funcInt * Count());
    }

    // Distribution1D Public Data
    std::vector<Float> func, cdf;
    Float funcInt;
};



class Distribution2D {
  public:
    // Distribution2D Public Methods
    Distribution2D(const Float *data, int nu, int nv);
    Point2f SampleContinuous(const Point2f &u, Float *pdf) const {
        Float pdfs[2];
        int v;
        Float d1 = pMarginal->SampleContinuous(u[1], &pdfs[1], &v);
        Float d0 = pConditionalV[v]->SampleContinuous(u[0], &pdfs[0]);
        *pdf = pdfs[0] * pdfs[1];
        return Point2f(d0, d1);
    }
    Float Pdf(const Point2f &p) const {
        int iu = glm::clamp(int(p[0] * pConditionalV[0]->Count()), 0,
                       pConditionalV[0]->Count() - 1);
        int iv =
            glm::clamp(int(p[1] * pMarginal->Count()), 0, pMarginal->Count() - 1);
        return pConditionalV[iv]->func[iu] / pMarginal->funcInt;
    }

  private:
    // Distribution2D Private Data
    std::vector<std::unique_ptr<Distribution1D>> pConditionalV;
    std::unique_ptr<Distribution1D> pMarginal;
};

#endif // DISTRIBUTION2D_H

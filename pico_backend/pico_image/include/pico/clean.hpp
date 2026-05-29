#pragma once
// Classic Hogbom CLEAN. The dirty image is searched iteratively for its peak
// pixel; a scaled, shifted copy of the PSF is subtracted; the position and
// flux are recorded as a CLEAN component. After convergence the components
// are convolved with a Gaussian fit to the PSF main lobe and added back to
// the residual.

#include <vector>
#include <cstdint>

namespace pico {

struct CleanComponent { int ix, iy; double flux; };

struct CleanParams {
    int    niter           = 1000;
    float  gain            = 0.1f;
    float  threshold_sigma = 3.0f;  // stop when |peak| < threshold_sigma * noise,
                                    // noise = robust (MAD) sigma of the dirty map.
};

struct CleanOut {
    std::vector<double>          residual;   // n1*n2
    std::vector<double>          restored;   // n1*n2 (residual + clean*beam)
    std::vector<CleanComponent>  components;
    double bmaj_pix = 0, bmin_pix = 0, bpa_rad = 0; // fitted Gaussian beam
};

// dirty, psf:  size n1*n2, row-major (y outer, x inner).
// PSF must be centered at (n1/2, n2/2) with peak = 1.
int hogbom_clean(const std::vector<double>& dirty,
                 const std::vector<double>& psf,
                 int n1, int n2,
                 const CleanParams& p,
                 CleanOut& out);

} // namespace pico

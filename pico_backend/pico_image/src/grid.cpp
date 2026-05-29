// pico_image/src/grid.cpp
// MFS gridder built on FINUFFT 2D type-1. Two NUFFT calls per image:
//   - dirty:  c_j = vis_j * wt_j      → ifft → dirty(x,y)
//   - PSF :   c_j = wt_j              → ifft → psf  (x,y)
//
// The "MFS" part lives entirely in how the caller fills GridSamples:
// we push (u_lambda_chan, v_lambda_chan, vis_chan) for every channel of
// every record onto the SAME sample list. All channels then land on the
// SAME grid via one NUFFT call, with each channel's uvw already in its
// own wavelengths — i.e. multifrequency synthesis at the image cellsize.

#include "pico/grid.hpp"
#include <complex>
#include <cmath>
#include <cstdio>
#include <vector>

#include <finufft.h>

namespace pico {

void to_finufft_coords(std::vector<double>& u_lam, std::vector<double>& v_lam,
                       double cellsize_rad) {
    // FINUFFT non-uniform points must lie in [-π, π).  Map by
    //     x = 2π · u_λ · Δθ        (Δθ = cellsize in radians)
    // The resulting image has pixel spacing Δθ and N1×N2 modes covers
    // the field [-N1·Δθ/2, +N1·Δθ/2).
    //
    // The u axis is NEGATED: FINUFFT type-1 (+isign) puts l increasing with
    // pixel index, but our FITS header uses CDELT1 < 0 (CASA/WSClean RA-east-
    // left convention) where RA — hence l — DECREASES with pixel index. Flip
    // u so the gridded sky matches the WCS (else the image is mirrored in x).
    // v is left as-is: CDELT2 > 0 already agrees with the gridder.
    const double s = 2.0 * M_PI * cellsize_rad;
    const double pi = M_PI;
    for (auto& u : u_lam) {
        double x = -u * s;
        x = std::fmod(x + pi, 2.0 * pi);
        if (x < 0) x += 2.0 * pi;
        u = x - pi;
    }
    for (auto& v : v_lam) {
        double y = v * s;
        y = std::fmod(y + pi, 2.0 * pi);
        if (y < 0) y += 2.0 * pi;
        v = y - pi;
    }
}

int run_finufft_dirty(const GridSamples& g, int n1, int n2, double tol,
                      DirtyOut& out) {
    if (g.size() == 0) {
        std::fprintf(stderr, "grid: no samples to image\n");
        return -1;
    }
    out.n1 = n1; out.n2 = n2;
    out.grid.assign(static_cast<std::size_t>(n1) * n2, std::complex<double>{0,0});
    out.img .assign(static_cast<std::size_t>(n1) * n2, 0.0);

    finufft_opts opts;
    finufft_default_opts(&opts);
    opts.upsampfac = 1.25;
    opts.modeord   = 0;          // CMCL: zero-frequency mode at (n1/2, n2/2)
                                 // — matches Hogbom/PSF-fit centering assumption.

    // We pass non-const pointers (finufft API). Local copies because g is const.
    std::vector<double> u = g.u_lambda;
    std::vector<double> v = g.v_lambda;
    std::vector<std::complex<double>> c = g.c;

    {
        std::size_t bad_uv = 0, bad_c = 0;
        for (std::size_t i = 0; i < c.size(); ++i) {
            if (!std::isfinite(u[i]) || !std::isfinite(v[i])) ++bad_uv;
            if (!std::isfinite(c[i].real()) || !std::isfinite(c[i].imag())) ++bad_c;
        }
        if (bad_uv || bad_c)
            std::fprintf(stderr,
                "finufft input: %zu/%zu uv non-finite, %zu/%zu c non-finite\n",
                bad_uv, c.size(), bad_c, c.size());
    }

    const int ier = finufft2d1(static_cast<int64_t>(c.size()), u.data(), v.data(),
                               c.data(), +1, tol,
                               n1, n2, out.grid.data(), &opts);
    if (ier != 0) {
        std::fprintf(stderr, "finufft2d1 returned ier=%d\n", ier);
        return ier;
    }

    // sum of weights for normalisation
    out.sumw = 0.0;
    for (double w : g.wt) out.sumw += w;
    const double norm = (out.sumw > 0.0) ? (1.0 / out.sumw) : 1.0;

    for (std::size_t i = 0; i < out.grid.size(); ++i)
        out.img[i] = out.grid[i].real() * norm;

    return 0;
}

} // namespace pico

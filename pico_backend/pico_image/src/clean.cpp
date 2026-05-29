// pico_image/src/clean.cpp
// Hogbom CLEAN. The dirty image and PSF are real, n1×n2, row-major
// (idx = iy*n1 + ix). The PSF is assumed centred at (n1/2, n2/2) with
// peak=1; subtraction wraps using offsets relative to the PSF centre.
//
// Restoring beam: a Gaussian fit to the PSF main lobe — width estimated
// from the half-maximum contour ellipse. We currently fit only an
// axis-aligned width (good enough for a basic skeleton).

#include "pico/clean.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace pico {

namespace {

inline std::size_t idx(int ix, int iy, int n1) {
    return static_cast<std::size_t>(iy) * n1 + ix;
}

// Robust (MAD) noise estimate of an image: 1.4826 * median(|x - median(x)|).
// Robust to bright sources/sidelobes, unlike a plain RMS.
double mad_sigma(const std::vector<double>& a) {
    if (a.empty()) return 0.0;
    std::vector<double> tmp(a);
    auto med = [](std::vector<double>& v) -> double {
        const std::size_t m = v.size() / 2;
        std::nth_element(v.begin(), v.begin() + m, v.end());
        return v[m];
    };
    const double median = med(tmp);
    for (double& x : tmp) x = std::fabs(x - median);
    return 1.4826 * med(tmp);
}

// Robust restoring-beam fit: least-squares 2D Gaussian to the PSF main lobe
// via log-quadratic linearisation. ln(psf) = a + d·x² + e·y² + f·xy over pixels
// above half-max gives the inverse covariance directly, from which we recover
// FWHM major/minor and position angle (handles rotated beams; gives a real BPA).
// Falls back to an axis-aligned half-max estimate if the fit is ill-conditioned.
void fit_psf_gauss(const std::vector<double>& psf, int n1, int n2,
                   double& bmaj_pix, double& bmin_pix, double& bpa_rad) {
    const int cx = n1 / 2, cy = n2 / 2;
    const double peak = psf[idx(cx, cy, n1)];

    // Axis-aligned half-max widths: used both to size the fit box and as a
    // fallback if the LSQ fit fails.
    auto hwhm_along = [&](int dx, int dy) -> double {
        double hw = 0.5;
        for (int k = 1; k < std::min(n1, n2) / 2; ++k) {
            const int x = cx + k * dx, y = cy + k * dy;
            if (x < 0 || x >= n1 || y < 0 || y >= n2) break;
            if (psf[idx(x, y, n1)] < 0.5 * peak) {
                const double prev = psf[idx(x - dx, y - dy, n1)];
                const double cur  = psf[idx(x, y, n1)];
                const double frac = (0.5 * peak - cur) / (prev - cur + 1e-30);
                hw = double(k) - frac;
                break;
            }
        }
        return hw;
    };
    const double hx = hwhm_along(1, 0);
    const double hy = hwhm_along(0, 1);
    auto fallback = [&]() {
        bmaj_pix = 2.0 * std::max(hx, hy);
        bmin_pix = 2.0 * std::min(hx, hy);
        bpa_rad  = (hx >= hy) ? 0.0 : M_PI / 2.0;
    };

    if (peak <= 0.0) { fallback(); return; }

    const int box = std::max(3, static_cast<int>(std::ceil(3.0 * std::max(hx, hy))));
    // Normal equations for basis {1, x², y², xy} (peak fixed at the centre).
    double S[4][4] = {{0}};
    double rhs[4]  = {0};
    int npix = 0;
    for (int dy = -box; dy <= box; ++dy) {
        const int y = cy + dy;
        if (y < 0 || y >= n2) continue;
        for (int dx = -box; dx <= box; ++dx) {
            const int x = cx + dx;
            if (x < 0 || x >= n1) continue;
            const double v = psf[idx(x, y, n1)];
            if (v <= 0.5 * peak) continue;           // main lobe only
            const double z = std::log(v / peak);
            const double g[4] = { 1.0, double(dx)*dx, double(dy)*dy, double(dx)*dy };
            for (int i = 0; i < 4; ++i) {
                rhs[i] += g[i] * z;
                for (int j = 0; j < 4; ++j) S[i][j] += g[i] * g[j];
            }
            ++npix;
        }
    }
    if (npix < 6) { fallback(); return; }

    // Solve 4x4 S c = rhs by Gaussian elimination with partial pivoting.
    double c[4];
    {
        double M[4][5];
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) M[i][j] = S[i][j];
            M[i][4] = rhs[i];
        }
        bool ok = true;
        for (int col = 0; col < 4 && ok; ++col) {
            int piv = col;
            for (int r = col + 1; r < 4; ++r)
                if (std::fabs(M[r][col]) > std::fabs(M[piv][col])) piv = r;
            if (std::fabs(M[piv][col]) < 1e-300) { ok = false; break; }
            for (int j = 0; j <= 4; ++j) std::swap(M[col][j], M[piv][j]);
            for (int r = 0; r < 4; ++r) {
                if (r == col) continue;
                const double f = M[r][col] / M[col][col];
                for (int j = col; j <= 4; ++j) M[r][j] -= f * M[col][j];
            }
        }
        if (!ok) { fallback(); return; }
        for (int i = 0; i < 4; ++i) c[i] = M[i][4] / M[i][i];
    }

    // c = {a, d, e, f}. Inverse covariance M_inv = [[-2d, -f], [-f, -2e]].
    const double d = c[1], e = c[2], f = c[3];
    const double detMinv = (-2.0 * d) * (-2.0 * e) - f * f;   // = 4de - f²
    if (!(detMinv > 0.0) || d >= 0.0 || e >= 0.0) { fallback(); return; }
    // Covariance Sigma = M_inv^{-1}.
    const double Sxx = (-2.0 * e) / detMinv;
    const double Syy = (-2.0 * d) / detMinv;
    const double Sxy =  f         / detMinv;

    const double half_sum = 0.5 * (Sxx + Syy);
    const double rad = std::sqrt(0.25 * (Sxx - Syy) * (Sxx - Syy) + Sxy * Sxy);
    const double lam_maj = half_sum + rad;   // larger eigenvalue = major variance
    const double lam_min = half_sum - rad;
    if (!(lam_min > 0.0)) { fallback(); return; }

    const double FWHM = 2.0 * std::sqrt(2.0 * std::log(2.0));   // 2.3548
    bmaj_pix = FWHM * std::sqrt(lam_maj);
    bmin_pix = FWHM * std::sqrt(lam_min);
    bpa_rad  = 0.5 * std::atan2(2.0 * Sxy, Sxx - Syy);          // major axis from +x
}

void add_restoring_beam(std::vector<double>& restored, int n1, int n2,
                        const std::vector<CleanComponent>& comps,
                        double bmaj_pix, double bmin_pix, double bpa_rad) {
    // Build a small Gaussian kernel sized to ~3*FWHM and stamp it.
    const double sx = bmaj_pix / (2.0 * std::sqrt(2.0 * std::log(2.0)));
    const double sy = bmin_pix / (2.0 * std::sqrt(2.0 * std::log(2.0)));
    const int radius = static_cast<int>(std::ceil(3.0 * std::max(sx, sy)));
    const double cos_pa = std::cos(bpa_rad), sin_pa = std::sin(bpa_rad);
    for (const auto& cc : comps) {
        for (int dy = -radius; dy <= radius; ++dy) {
            for (int dx = -radius; dx <= radius; ++dx) {
                const int ix = cc.ix + dx, iy = cc.iy + dy;
                if (ix < 0 || ix >= n1 || iy < 0 || iy >= n2) continue;
                const double xr =  dx * cos_pa + dy * sin_pa;
                const double yr = -dx * sin_pa + dy * cos_pa;
                const double g = std::exp(-0.5 * ((xr*xr)/(sx*sx) + (yr*yr)/(sy*sy)));
                restored[idx(ix, iy, n1)] += cc.flux * g;
            }
        }
    }
}

} // namespace

int hogbom_clean(const std::vector<double>& dirty,
                 const std::vector<double>& psf,
                 int n1, int n2,
                 const CleanParams& p,
                 CleanOut& out) {
    if (static_cast<int>(dirty.size()) != n1 * n2 ||
        static_cast<int>(psf.size())   != n1 * n2) {
        std::fprintf(stderr, "clean: dimension mismatch\n");
        return -1;
    }
    out.residual = dirty;
    out.components.clear();
    const int cx = n1 / 2, cy = n2 / 2;

    // Robust noise from the dirty map, then an absolute flux floor in Jy.
    // CLEAN stops once the brightest residual peak drops below it, so we never
    // dig into the noise (which is what produced scattered spurious components).
    const double sigma     = mad_sigma(dirty);
    const double flux_floor = static_cast<double>(p.threshold_sigma) * sigma;
    std::fprintf(stderr,
        "clean: noise(MAD)=%.6g Jy, threshold=%.2f sigma => stop at %.6g Jy\n",
        sigma, p.threshold_sigma, flux_floor);

    int it = 0;
    for (; it < p.niter; ++it) {
        // find peak |residual|
        std::size_t pk_i = 0;
        double pk_v = 0.0;
        for (std::size_t i = 0; i < out.residual.size(); ++i) {
            const double a = std::fabs(out.residual[i]);
            if (a > pk_v) { pk_v = a; pk_i = i; }
        }
        const int pk_x = static_cast<int>(pk_i % n1);
        const int pk_y = static_cast<int>(pk_i / n1);
        const double pk_val = out.residual[pk_i];
        if (pk_v < flux_floor) break;

        const double comp = p.gain * pk_val;
        out.components.push_back({pk_x, pk_y, comp});

        // subtract gain*peak*PSF shifted so PSF centre lands at pk
        for (int iy = 0; iy < n2; ++iy) {
            const int sy = iy - pk_y + cy;
            if (sy < 0 || sy >= n2) continue;
            for (int ix = 0; ix < n1; ++ix) {
                const int sx = ix - pk_x + cx;
                if (sx < 0 || sx >= n1) continue;
                out.residual[idx(ix, iy, n1)] -= comp * psf[idx(sx, sy, n1)];
            }
        }
    }

    std::fprintf(stderr, "clean: ran %d/%d iterations (%zu components)\n",
                 it, p.niter, out.components.size());

    fit_psf_gauss(psf, n1, n2, out.bmaj_pix, out.bmin_pix, out.bpa_rad);
    out.restored = out.residual;
    add_restoring_beam(out.restored, n1, n2, out.components,
                       out.bmaj_pix, out.bmin_pix, out.bpa_rad);
    return 0;
}

} // namespace pico

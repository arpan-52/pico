#pragma once
// MFS gridder built on FINUFFT 2D type-1 (nonuniform → uniform). For every
// (baseline, channel, integration) sample we push:
//     u_lambda = u_seconds * freq_channel
//     v_lambda = v_seconds * freq_channel
//     c       = vis * weight
// onto a single common Fourier-mode grid sized to the image (N1 × N2).
//
// FINUFFT type-1 evaluates  F[k] = sum_j c_j exp(+i k·x_j) ; we take that as
// the (unnormalised) dirty image's Fourier transform — i.e. the gridded
// visibility plane. A subsequent inverse FFT (built into FINUFFT) hands us
// the dirty image.
//
// We do TWO finufft calls per image: one with c = vis*wt → dirty image,
// one with c = wt → PSF (point spread function).

#include "pico/types.hpp"
#include <vector>
#include <complex>
#include <cmath>

namespace pico {

struct GridSamples {
    // FINUFFT wants doubles for the points. Visibilities are float; we
    // convert at push time. All arrays are the same length.
    std::vector<double>           u_lambda;    // u in wavelengths at ref freq grid
    std::vector<double>           v_lambda;
    std::vector<std::complex<double>> c;       // vis * wt for dirty, wt for PSF
    std::vector<double>           wt;          // raw weight (for PSF + sumw)

    void reserve(std::size_t n) {
        u_lambda.reserve(n); v_lambda.reserve(n);
        c.reserve(n); wt.reserve(n);
    }
    std::size_t size() const { return u_lambda.size(); }
};

// Push one channel's worth of one baseline's worth of one record into the
// grid samples, scaling uvw_sec by freq_hz to get uvw_lambda. Skips flagged
// (wt <= 0) samples.
inline void push_sample(GridSamples& g, const UvwSec& uvw_sec,
                        double freq_hz, const Vis& v) {
    if (v.wt <= 0.0f) return;
    // Reject any non-finite sample (NaN-decoded halfs, calib divide-by-zero,
    // etc.). One NaN entering FINUFFT spreads across the whole image grid.
    if (!std::isfinite(v.r) || !std::isfinite(v.i) || !std::isfinite(v.wt)) return;
    const double u = uvw_sec.u * freq_hz;
    const double vv = uvw_sec.v * freq_hz;
    if (!std::isfinite(u) || !std::isfinite(vv)) return;
    g.u_lambda.push_back(u);
    g.v_lambda.push_back(vv);
    g.c.emplace_back(double(v.r) * double(v.wt),
                     double(v.i) * double(v.wt));
    g.wt.push_back(double(v.wt));
}

// Convert u,v (wavelengths) → FINUFFT non-uniform-point coordinates in
// [-pi, pi). Image cellsize sets the conversion:
//     x_finufft = 2*pi * u_lambda * cellsize_rad
// (Standard radio-astronomy MFS convention.)
void to_finufft_coords(std::vector<double>& u_lam, std::vector<double>& v_lam,
                       double cellsize_rad);

struct DirtyOut {
    int n1 = 0, n2 = 0;
    std::vector<std::complex<double>> grid;  // size n1*n2, FINUFFT output
    std::vector<double>               img;   // real dirty image (n1*n2)
    double sumw = 0.0;
};

// Robust MAD clip over the sample range [begin, size). Ports svfits clip()
// (svsubs.c:1322), which runs per (file,slice) on that slice's dumped
// visibilities — call this after each slice with `begin` = bag size before
// the slice, so the med/MAD scope matches svfits. One median + MAD over the
// range's amplitudes (|c|/wt), drop samples with |amp - med| > thresh*mad
// from BOTH bags in lockstep (they must be the same length / index-aligned).
// Two-sided, raw MAD — matches svfits robust_stats. Returns flagged count.
std::size_t global_clip_samples(GridSamples& dirty, GridSamples& psf,
                                double thresh, std::size_t begin = 0);

// Run FINUFFT 2D type-1 on samples g (already in [-pi,pi)) into a (n1 x n2)
// grid, then take the magnitude/real part as the dirty image. tol is the
// FINUFFT tolerance (e.g. 1e-6).
int run_finufft_dirty(const GridSamples& g, int n1, int n2, double tol,
                      DirtyOut& out);

} // namespace pico

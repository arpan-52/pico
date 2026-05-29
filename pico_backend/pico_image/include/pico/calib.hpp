#pragma once
// Off-source mean spectrum / amplitude bandpass / baseline (mean off-source
// visibility) subtraction. Ports svfits make_bpass() + the per-record
// application in make_onechan_group().

#include "pico/types.hpp"
#include "pico/raw_io.hpp"
#include "pico/antsamp.hpp"
#include "pico/config.hpp"
#include <vector>
#include <cmath>

namespace pico {

struct Bandpass {
    int n_base   = 0;
    int n_chan   = 0;
    // amplitude bandpass per baseline (normalized so geomean = 1)
    std::vector<std::vector<float>>   abp;       // [n_base][n_chan]
    // mean off-source complex visibility per (baseline, channel)
    std::vector<std::vector<Complex>> off_src;   // [n_base][n_chan]
    int file_idx = -1;
    int slice    = -1;
};

// Build bandpass + off-source-mean from one slice of raw data, by averaging
// records that do NOT contain the burst's time-freq region.
// rbuf points to rec_per_slice * recl bytes as returned by read_slice().
int make_bandpass(const Config& cfg, const RawSet& rs, const AntSamp& as,
                  const void* rbuf, int idx, int slice, Bandpass& bp);

// Subtract bp.off_src and divide by bp.abp from the visibility (in place).
// Operates on one (baseline, channel) sample.
inline void apply_calib(Vis& v, const Bandpass& bp, int b, int c) {
    if (bp.abp.empty() || b >= bp.n_base) return;
    const Complex off = bp.off_src[b][c];
    if (std::isfinite(off.real()) && std::isfinite(off.imag())) {
        v.r -= off.real();
        v.i -= off.imag();
    }
    const float a = bp.abp[b][c];
    if (std::isfinite(a) && a > 0.0f) { v.r /= a; v.i /= a; }
    if (!std::isfinite(v.r) || !std::isfinite(v.i)) { v.wt = -1.0f; }
}

// MAD-based robust clip (ports svfits clip()). Marks v.wt = -1 if flagged.
int clip_record(std::vector<Vis>& chans, float mad_thresh);

} // namespace pico

// pico_image/src/calib.cpp
// Off-source bandpass + baseline construction (svfits make_bpass).
// Strategy mirrors svfits:
//   - For records that do NOT overlap the DM track, accumulate per-(baseline,
//     channel) sum/mean of complex visibilities → off_src[b][c]
//   - Per-baseline mean amplitude over channels → normalize to 1 → abp[b][c]
//   - Median absolute deviation flagger for residual RFI (clip)

#include "pico/calib.hpp"
#include "pico/half.hpp"
#include "pico/dm.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

namespace pico {

namespace {

inline void decode_vis(const float* slot, Vis& v) {
    // svfits packs (re, im) as two uint16 halfs inside the float bit pattern.
    std::uint32_t bits;
    std::memcpy(&bits, slot, sizeof(bits));
    const std::uint16_t re_h = static_cast<std::uint16_t>(bits & 0xFFFFu);
    const std::uint16_t im_h = static_cast<std::uint16_t>((bits >> 16) & 0xFFFFu);
    v.r  = half_to_float(re_h);
    v.i  = half_to_float(im_h);
    v.wt = 1.0f;
}

} // namespace

int make_bandpass(const Config& cfg, const RawSet& rs, const AntSamp& as,
                  const void* rbuf, int idx, int slice, Bandpass& bp) {
    const int nch  = rs.channels;
    const int nbase = as.nbase;
    bp.n_base = nbase;
    bp.n_chan = nch;
    bp.abp.assign(nbase, std::vector<float>(nch, 1.0f));
    bp.off_src.assign(nbase, std::vector<Complex>(nch, Complex{0,0}));
    bp.file_idx = idx;
    bp.slice    = slice;

    // Iterate records in this slice. For each record, compute its time and
    // check whether the burst lies in this record's [cs,ce) range. If not,
    // accumulate into off_src.
    const auto* base = static_cast<const float*>(rbuf);
    const double dt_rec = rs.t_slice / rs.rec_per_slice;
    const double t_slice_start = idx * rs.t_slice + slice * rs.slice_interval;

    std::vector<int> n_off(nbase * nch, 0);

    for (int r = 0; r < rs.rec_per_slice; ++r) {
        const double trec = t_slice_start + r * dt_rec;
        int cs = 0, ce = 0;
        burst_chans_for_record(cfg, trec, dt_rec, nch, &cs, &ce);
        // off-source = channels NOT in [cs, ce)
        const float* rec_base = base + static_cast<std::ptrdiff_t>(r) * (nch * nbase);
        for (int b = 0; b < nbase; ++b) {
            const float* row = rec_base + static_cast<std::ptrdiff_t>(b) * nch;
            for (int c = 0; c < nch; ++c) {
                if (c >= cs && c < ce) continue;
                Vis v; decode_vis(row + c, v);
                bp.off_src[b][c] += Complex(v.r, v.i);
                ++n_off[b * nch + c];
            }
        }
    }
    // Normalise off_src to means
    for (int b = 0; b < nbase; ++b) {
        for (int c = 0; c < nch; ++c) {
            const int n = n_off[b * nch + c];
            if (n > 0) bp.off_src[b][c] /= static_cast<float>(n);
        }
    }
    // Amplitude bandpass: |off_src[b][c]| normalised per-baseline so that
    // its geomean across channels is 1
    for (int b = 0; b < nbase; ++b) {
        double logsum = 0.0; int n = 0;
        for (int c = 0; c < nch; ++c) {
            const float a = std::abs(bp.off_src[b][c]);
            if (a > 0) { logsum += std::log(double(a)); ++n; bp.abp[b][c] = a; }
            else       { bp.abp[b][c] = 1.0f; }
        }
        if (n > 0) {
            const float gm = static_cast<float>(std::exp(logsum / n));
            if (gm > 0) for (int c = 0; c < nch; ++c) bp.abp[b][c] /= gm;
        }
    }
    return 0;
}

int clip_record(std::vector<Vis>& ch, float mad_thresh) {
    if (ch.empty()) return 0;
    // Compute amplitude vector
    std::vector<float> amp; amp.reserve(ch.size());
    for (const auto& v : ch) if (v.wt > 0) amp.push_back(std::hypot(v.r, v.i));
    if (amp.size() < 4) return 0;
    std::nth_element(amp.begin(), amp.begin() + amp.size()/2, amp.end());
    const float med = amp[amp.size() / 2];
    // MAD
    std::vector<float> dev; dev.reserve(amp.size());
    for (float a : amp) dev.push_back(std::fabs(a - med));
    std::nth_element(dev.begin(), dev.begin() + dev.size()/2, dev.end());
    const float mad = dev[dev.size() / 2];
    if (mad <= 0) return 0;
    const float thr = med + mad_thresh * 1.4826f * mad;
    int flagged = 0;
    for (auto& v : ch) {
        if (std::hypot(v.r, v.i) > thr) { v.wt = -1.0f; ++flagged; }
    }
    return flagged;
}

} // namespace pico

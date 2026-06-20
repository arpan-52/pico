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
#include <cstdio>
#include <cstdlib>
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
    const double t_slice_start = rs.files[idx].t_start + slice * rs.slice_interval;

    // Per-record burst channel range, so we know which records are off-source
    // for a given channel.
    std::vector<int> rec_cs(rs.rec_per_slice), rec_ce(rs.rec_per_slice);
    for (int r = 0; r < rs.rec_per_slice; ++r) {
        const double trec = t_slice_start + r * dt_rec;
        burst_chans_for_record(cfg, trec, dt_rec, nch, &rec_cs[r], &rec_ce[r]);
    }

    if (!cfg.do_band && !cfg.do_base) return 0;  // nominal abp=1/off=0 stand

    // Per-channel burst record span, exactly as svfits (svsubs.c:1478-1488):
    // channel c is in-burst for record r iff start_chan < c < end_chan
    // (STRICT both ends), and the off-source region excludes the full
    // contiguous span [r0,r1] of such records.
    std::vector<int> ch_r0(nch, rs.rec_per_slice), ch_r1(nch, -1);
    for (int r = 0; r < rs.rec_per_slice; ++r)
        for (int c = rec_cs[r] + 1; c < rec_ce[r]; ++c) {
            if (r < ch_r0[c]) ch_r0[c] = r;
            if (r > ch_r1[c]) ch_r1[c] = r;
        }

    // off_src per (baseline, channel) = arithmetic MEAN of the off-burst
    // records (svfits make_bpass, svsubs.c:1494-1499). The MEAN — not the
    // median — is what removes the steady coherent component (continuum /
    // persistent correlated RFI) common to on- and off-burst samples; svfits
    // subtracts exactly this and the reference image comes out blank when no
    // burst is present. abp[b][c] = MEAN amplitude <sqrt(re^2+im^2)> over the
    // same records (svsubs.c:1495,1500). NOTE: abp must be <|V|> (mean of the
    // amplitudes), NOT |<V>| (= |off_src|): for noise |<V>| -> 0 and dividing
    // by it explodes the visibilities. n==0 channels are flagged (-1) and
    // filled by neighbour-carry below (svsubs.c:1502-1503,1529-1530).
    #pragma omp parallel for num_threads(cfg.num_threads)
    for (int b = 0; b < nbase; ++b) {
        // Skip autocorrelations and user-ANTMASK-excluded baselines — never
        // imaged, only pollute off_src. (b indexes the on-disk record block.)
        if (as.baseline[b].s0.ant_id == as.baseline[b].s1.ant_id) continue;
        if (as.baseline[b].drop) continue;
        for (int c = 0; c < nch; ++c) {
            double sre = 0, sim = 0, samp = 0; long n = 0;
            for (int r = 0; r < rs.rec_per_slice; ++r) {
                if (r >= ch_r0[c] && r <= ch_r1[c]) continue;  // burst span
                const float* slot = base +
                    (static_cast<std::ptrdiff_t>(r) * (nch * nbase) +
                     static_cast<std::ptrdiff_t>(b) * nch + c);
                Vis v; decode_vis(slot, v);
                // Guard non-finite decoded halfs (svfits svsubs.c:1493).
                if (!std::isfinite(v.r) || !std::isfinite(v.i)) continue;
                sre += v.r; sim += v.i; samp += std::hypot(v.r, v.i); ++n;
            }
            if (n > 0) {
                bp.off_src[b][c] = Complex(float(sre / n), float(sim / n));
                bp.abp[b][c]     = float(samp / n);
            } else {
                bp.off_src[b][c] = Complex{0, 0};
                bp.abp[b][c]     = -1.0f;  // flag for neighbour-carry below
            }
        }
    }
    // off_src bias diagnostic: this complex constant is subtracted from EVERY
    // burst sample, so a large coherent +ve real here drives a central source.
    {
        double sre = 0, sim = 0, samp = 0, amax = 0; std::size_t n = 0;
        for (int b = 0; b < nbase; ++b)
            for (int c = 0; c < nch; ++c) {
                const Complex o = bp.off_src[b][c];
                const double a = std::abs(o);
                sre += o.real(); sim += o.imag(); samp += a;
                if (a > amax) amax = a; ++n;
            }
        if (n > 0)
            std::fprintf(stderr,
                "make_bandpass[file=%d slice=%d]: off_src mean re=%.6g im=%.6g "
                "|.|=%.6g max|.|=%.6g over %zu (base,chan)\n",
                idx, slice, sre / n, sim / n, samp / n, amax, n);
    }
    // do_band off → nominal flat bandpass, skip interpolation/normalisation
    // (svfits svsubs.c:1512-1519).
    if (!cfg.do_band) {
        for (int b = 0; b < nbase; ++b)
            for (int c = 0; c < nch; ++c) bp.abp[b][c] = 1.0f;
        return 0;
    }

    // Amplitude bandpass normalisation (svfits svsubs.c:1521-1536): carry the
    // first valid value backward, forward-fill flagged channels, then divide
    // by the per-baseline median across channels so the bandpass sits near 1.
    for (int b = 0; b < nbase; ++b) {
        if (as.baseline[b].s0.ant_id == as.baseline[b].s1.ant_id) continue;
        if (as.baseline[b].drop) continue;
        auto& abp = bp.abp[b];
        if (abp[0] < 0.0f)
            for (int c = 1; c < nch; ++c) if (abp[c] > 0.0f) { abp[0] = abp[c]; break; }
        if (abp[0] < 0.0f) {                       // all channels flagged
            for (int c = 0; c < nch; ++c) abp[c] = 1.0f;
            continue;
        }
        for (int c = 1; c < nch; ++c) if (abp[c] < 0.0f) abp[c] = abp[c - 1];
        // svfits quick_median averages the two middle elements for even n
        // (stats.c:140-151) — match it, nch is 4096.
        std::vector<float> tmp(abp.begin(), abp.end());
        std::nth_element(tmp.begin(), tmp.begin() + nch / 2, tmp.end());
        float medc = tmp[nch / 2];
        if (nch % 2 == 0) {
            const float lo = *std::max_element(tmp.begin(), tmp.begin() + nch / 2);
            medc = 0.5f * (medc + lo);
        }
        if (medc > 0.0f) for (int c = 0; c < nch; ++c) abp[c] /= medc;
    }

    // do_base off → no off-source subtraction (svfits svsubs.c:1538-1545).
    if (!cfg.do_base)
        for (int b = 0; b < nbase; ++b)
            for (int c = 0; c < nch; ++c) bp.off_src[b][c] = Complex{0, 0};

    // --- per-channel bandpass debug (env PICO_BP_DEBUG) -------------------
    // For the first cross baseline of the first slice: recompute the RAW (pre-
    // normalization) abp = <|raw|> over off-source records, the per-baseline
    // median used to normalize, the abp distribution, and the channels with the
    // smallest normalized abp (the ÷abp blow-up). This shows whether the tiny
    // abp comes from a genuinely quiet channel, too-few off-source records, or
    // an inflated median.
    if (std::getenv("PICO_BP_DEBUG")) {
        static bool done = false;
        if (!done) {
            done = true;
            int bd = -1;
            for (int b = 0; b < nbase; ++b)
                if (as.baseline[b].s0.ant_id != as.baseline[b].s1.ant_id) { bd = b; break; }
            if (bd >= 0) {
                const auto* base = static_cast<const float*>(rbuf);
                std::vector<float> pre(nch, -1.f); std::vector<int> noff(nch, 0);
                for (int c = 0; c < nch; ++c) {
                    double s = 0; int n = 0;
                    for (int r = 0; r < rs.rec_per_slice; ++r) {
                        if (r >= ch_r0[c] && r <= ch_r1[c]) continue;
                        Vis v; decode_vis(base + (std::ptrdiff_t)r*(nch*nbase) + (std::ptrdiff_t)bd*nch + c, v);
                        if (std::isfinite(v.r) && std::isfinite(v.i)) { s += std::hypot(v.r, v.i); ++n; }
                    }
                    noff[c] = n; if (n > 0) pre[c] = float(s / n);
                }
                // forward-fill flagged then median (same as the real path)
                std::vector<float> ff(pre);
                if (ff[0] < 0) for (int c = 1; c < nch; ++c) if (ff[c] > 0) { ff[0] = ff[c]; break; }
                for (int c = 1; c < nch; ++c) if (ff[c] < 0) ff[c] = ff[c-1];
                std::vector<float> tmp(ff); std::nth_element(tmp.begin(), tmp.begin()+nch/2, tmp.end());
                float med = tmp[nch/2];
                if (nch%2==0) med = 0.5f*(med + *std::max_element(tmp.begin(), tmp.begin()+nch/2));
                // normalized abp = ff/med ; distribution + smallest
                std::vector<float> nrm(nch); for (int c=0;c<nch;++c) nrm[c]=ff[c]/med;
                int lt01=0,lt05=0,gt2=0,flagged=0;
                for (int c=0;c<nch;++c){ if(pre[c]<0)++flagged; float a=nrm[c]; if(a<0.1f)++lt01; else if(a<0.5f)++lt05; else if(a>2.f)++gt2; }
                std::fprintf(stderr,
                    "BP_DEBUG b=%d a%d-a%d slice=%d: median(abp_pre)=%.3f  flagged(n_off==0)=%d  "
                    "abp_norm dist <0.1:%d 0.1-0.5:%d >2:%d /%d\n",
                    bd, as.baseline[bd].s0.ant_id, as.baseline[bd].s1.ant_id, slice,
                    med, flagged, lt01, lt05, gt2, nch);
                std::vector<int> ord(nch); for(int c=0;c<nch;++c)ord[c]=c;
                std::sort(ord.begin(),ord.end(),[&](int a,int b){return nrm[a]<nrm[b];});
                std::fprintf(stderr,"  smallest-abp channels (c: n_off abp_pre abp_norm off_re):\n");
                for (int k=0;k<8;++k){ int c=ord[k];
                    std::fprintf(stderr,"   c=%4d  n_off=%2d  abp_pre=%9.4f  abp_norm=%.5f  off_re=%9.3f\n",
                        c, noff[c], pre[c], nrm[c], bp.off_src[bd][c].real()); }
            }
        }
    }
    return 0;
}

} // namespace pico

// pico_image/src/dm.cpp
// DM-aware channel-range and record-range selection. Ports
//   svfits get_chan_num()  — per-record (cs,ce) channel range
//   svfits get_rec_num()   — per-file (start_rec,n_rec) selection
//
// Time convention: all times are seconds since rs.mjd_ref. Each file's first
// record time comes from its own slice-0 timestamp + embedded file index
// (rs.files[i].t_start, set in raw_io — svfits get_slice_time convention; the
// supplied file order is NOT trusted). Files are time-multiplexed at the
// slice level; within a file, consecutive slices are slice_interval =
// nfile*t_slice apart.
//
// Cold-plasma dispersion delay (s at freq_hz):
//     Δt = K0 · DM / f_GHz²    with K0 = 4.148808e-3 s·GHz²·(pc/cm³)⁻¹

#include "pico/dm.hpp"
#include <cmath>
#include <algorithm>
#include <cstdio>

namespace pico {

// svfits uses 4.15e-3 (svsubs.c:35); match it exactly so record/channel
// selection is identical record-for-record.
static constexpr double K0 = 4.15e-3;

int burst_chans_for_record(const Config& cfg, double trec, double integ_sec,
                           int nchan_full, int* cs, int* ce) {
    const auto& b = cfg.burst;
    const double fb   = b.freq * 1e-9;                       // GHz
    const double freq = cfg.freq0_hz * 1e-9;                 // GHz (chan 0)
    const double cw   = std::fabs(cfg.freq1_hz - cfg.freq0_hz) / cfg.nchan * 1e-9;
    const double net_sign = (cfg.freq1_hz >= cfg.freq0_hz) ? +1.0 : -1.0;
    const double df   = cw * net_sign;

    const double tr1 = trec;
    const double tr2 = trec + integ_sec;
    const double tb  = b.t;
    const double wd  = b.int_wd;
    const double DM  = b.DM;

    const double fn1 = 1.0 / (fb * fb) + (tr1 - tb - wd / 2.0) / (K0 * DM);
    const double fn2 = 1.0 / (fb * fb) + (tr2 - tb + wd / 2.0) / (K0 * DM);
    const double fs  = (fn1 > 0.0) ? 1.0 / std::sqrt(fn1) : -1.0;   // highest freq
    const double fe  = (fn2 > 0.0) ? 1.0 / std::sqrt(fn2) : -1.0;   // lowest freq

    int c0, c1;
    if (fs > 0.0 && fe > 0.0) {
        c0 = static_cast<int>(std::floor((fs - freq) / df));
        c1 = static_cast<int>(std::floor((fe - freq) / df));
    } else {
        c0 = c1 = nchan_full;
    }
    if (c0 > c1) std::swap(c0, c1);
    *cs = (c0 > 0) ? c0 : 0;
    *ce = (c1 < nchan_full) ? c1 : nchan_full;
    if (c1 < 0) *cs = *ce = nchan_full;
    return 0;
}

// Mirrors svfits's get_rec_num: t_burst = arrival of trailing edge of the
// intrinsic burst at the highest band frequency, then walk records forward
// until t_rec > t_burst + burst_width.
int compute_record_ranges(const Config& cfg, RawSet& rs) {
    const auto& b = cfg.burst;
    const double f_hi = std::max(cfg.freq0_hz, cfg.freq1_hz) * 1e-9;   // GHz
    const double fb   = b.freq * 1e-9;

    const double t_burst = b.t - b.int_wd / 2.0
                         + K0 * b.DM * (1.0 / (f_hi * f_hi) - 1.0 / (fb * fb));
    const double t_end   = t_burst + b.width;

    const double integ          = rs.t_slice / rs.rec_per_slice;
    const double slice_interval = rs.slice_interval;
    const double t_slice        = rs.t_slice;
    const int    rec_per_slice  = rs.rec_per_slice;

    std::fprintf(stderr,
        "dm: t_burst=%.6f s (since mjd_ref), t_end=%.6f s, integ=%.6f ms,"
        " t_slice=%.6f s, slice_interval=%.6f s\n",
        t_burst, t_end, integ * 1e3, t_slice, slice_interval);

    for (int i = 0; i < cfg.nfile; ++i) {
        auto& f = rs.files[i];
        const double t_start = f.t_start;  // from per-file timestamp + embedded idx

        int start_slice = static_cast<int>(std::floor((t_burst - t_start) / slice_interval));
        if (start_slice < 0) start_slice = 0;
        int rec_num;
        if (t_burst - (t_start + start_slice * slice_interval) > t_slice) {
            ++start_slice;
            rec_num = 0;
        } else {
            rec_num = static_cast<int>(std::floor(
                (t_burst - (t_start + start_slice * slice_interval)) / integ));
            if (rec_num < 0) rec_num = 0;
        }

        const int start_rec = start_slice * rec_per_slice + rec_num;
        double t_rec        = t_start + start_slice * slice_interval + rec_num * integ;
        f.b_start = t_rec;

        int r = 0;
        while (t_rec <= t_end) {
            if (rec_num == rec_per_slice - 1) {
                t_rec  += slice_interval - t_slice;  // jump to next slice on same file
                rec_num = 0;
            } else {
                t_rec  += integ;
                ++rec_num;
            }
            ++r;
        }
        f.start_rec = start_rec;
        f.n_rec     = r;
    }
    return 0;
}

} // namespace pico

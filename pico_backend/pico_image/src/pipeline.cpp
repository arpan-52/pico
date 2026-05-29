// pico_image/src/pipeline.cpp
// Top-level driver. Orchestrates: load configs → open raw files → DM-aware
// record selection → per-slice loop (parallel across files) → for each
// record: decode visibilities, bandpass + baseline subtract, MAD clip,
// compute uvw (apparent), rotate to J2000, push (u_λ, v_λ, vis*wt) into the
// MFS sample list, scaled per-channel by freq → after the entire burst is
// gathered run two FINUFFT 2D-type-1 calls (dirty + PSF), then Hogbom CLEAN,
// then write the FITS image.

#include "pico/pipeline.hpp"
#include "pico/antsamp.hpp"
#include "pico/raw_io.hpp"
#include "pico/calib.hpp"
#include "pico/dm.hpp"
#include "pico/uvw.hpp"
#include "pico/grid.hpp"
#include "pico/clean.hpp"
#include "pico/fits_out.hpp"
#include "pico/half.hpp"
#include <cstdio>
#include <cstring>
#include <complex>
#include <vector>
#include <cmath>
#include <cstdint>
#include <omp.h>

namespace pico {

namespace {

inline void decode_packed(const float* slot, Vis& v) {
    std::uint32_t bits;
    std::memcpy(&bits, slot, sizeof(bits));
    v.r  = half_to_float(static_cast<std::uint16_t>(bits & 0xFFFFu));
    v.i  = half_to_float(static_cast<std::uint16_t>((bits >> 16) & 0xFFFFu));
    v.wt = 1.0f;
}

} // namespace

// Ports svfits update_burst(): recompute burst.t (in seconds since mjd_ref)
// from burst.mjd, and burst.width as the full DM sweep across the band plus
// intrinsic width. Run after raw_io opens (we need mjd_ref).
static void update_burst(Config& cfg, double mjd_ref) {
    if (!cfg.update_burst) return;
    auto& b = cfg.burst;
    constexpr double K0 = 4.148808e-3;
    const double fb = b.freq * 1e-9;                                  // GHz
    const double f_hi = std::max(cfg.freq0_hz, cfg.freq1_hz) * 1e-9;
    const double f_lo = std::min(cfg.freq0_hz, cfg.freq1_hz) * 1e-9;
    const double th = -b.int_wd / 2.0
                    + K0 * b.DM * (1.0/(f_hi*f_hi) - 1.0/(fb*fb));
    const double tl =  b.int_wd / 2.0
                    + K0 * b.DM * (1.0/(f_lo*f_lo) - 1.0/(fb*fb));
    b.width = tl - th;
    // burst.t in seconds since mjd_ref (the time origin used everywhere else
    // in pico_image — t_start[idx] = idx*t_slice, trec measured from the same
    // origin). burst.mjd is the burst arrival time at burst.freq, in MJD.
    b.t = (b.mjd - mjd_ref) * 86400.0;
    std::fprintf(stderr,
        "update_burst: burst.mjd=%.10f mjd_ref=%.10f → burst.t=%.6f s ;"
        " burst.width=%.6f s (full DM sweep + intrinsic)\n",
        b.mjd, mjd_ref, b.t, b.width);
}

int run_pipeline(const Config& cfg_in) {
    Config cfg = cfg_in;  // local mutable copy (so update_burst can adjust)

    AntSamp as;
    if (load_antsamp(cfg.ant_hdr_path, as) < 0) return -1;
    rebuild_baselines(as, cfg.antmask);
    std::fprintf(stderr, "antsamp: %d antennas, %d samplers, %d baselines (RR+LL)\n",
                 __builtin_popcount(as.antmask), as.nsamp, as.nbase);

    RawSet rs;
    if (open_raw_set(cfg, as, rs) < 0) return -1;
    std::fprintf(stderr, "raw_io: mjd_ref=%.6f, recl=%zu bytes, %d rec/slice\n",
                 rs.mjd_ref, rs.recl, rs.rec_per_slice);

    update_burst(cfg, rs.mjd_ref);
    if (compute_record_ranges(cfg, rs) < 0) { close_raw_set(rs); return -1; }
    for (int i = 0; i < cfg.nfile; ++i) {
        std::fprintf(stderr, "  file[%d] %s : burst records [%d,%d)\n",
                     i, rs.files[i].path.c_str(),
                     rs.files[i].start_rec,
                     rs.files[i].start_rec + rs.files[i].n_rec);
    }

    const double ref_freq = resolve_ref_freq(cfg);
    const double ch_w     = (cfg.freq1_hz - cfg.freq0_hz) / cfg.nchan;
    std::fprintf(stderr, "MFS ref_freq = %.3f MHz; ch_width = %.3f kHz\n",
                 ref_freq * 1e-6, ch_w * 1e-3);

    // Final MFS sample bag (one per channel × baseline × record).
    GridSamples dirty_samples;
    GridSamples psf_samples;
    std::size_t est_samples = std::size_t(as.nbase) * cfg.nchan * 64;
    dirty_samples.reserve(est_samples);
    psf_samples  .reserve(est_samples);

    // Per-slice raw buffer (thread-local-able later)
    std::vector<float> rbuf(rs.rec_per_slice *
                            static_cast<std::size_t>(rs.recl / sizeof(float)));
    const int nch = rs.channels;

    const double dt_rec = rs.t_slice / rs.rec_per_slice;

    for (int i = 0; i < cfg.nfile; ++i) {
        auto& f = rs.files[i];
        if (f.n_rec <= 0) continue;
        const int slice0 = f.start_rec / rs.rec_per_slice;
        const int slice1 = (f.start_rec + f.n_rec + rs.rec_per_slice - 1) / rs.rec_per_slice;

        std::fprintf(stderr, "file[%d] iterating slices [%d,%d)\n", i, slice0, slice1);
        Bandpass bp;
        for (int sl = slice0; sl < slice1; ++sl) {
            std::fprintf(stderr, "  reading file=%d slice=%d ...\n", i, sl);
            if (read_slice(rs, i, sl, rbuf.data()) < 0) {
                std::fprintf(stderr, "  read_slice FAILED file=%d slice=%d (EOF or short read)\n", i, sl);
                break;
            }
            std::fprintf(stderr, "  read_slice ok file=%d slice=%d\n", i, sl);

            // Build bandpass once per slice (svfits behavior).
            if (cfg.do_band || cfg.do_base)
                make_bandpass(cfg, rs, as, rbuf.data(), i, sl, bp);

            const double t_slice_start = i * rs.t_slice + sl * rs.slice_interval;
            for (int r = 0; r < rs.rec_per_slice; ++r) {
                const int global_rec = sl * rs.rec_per_slice + r;
                if (global_rec < f.start_rec ||
                    global_rec >= f.start_rec + f.n_rec) continue;
                const double trec = t_slice_start + r * dt_rec;
                int cs = 0, ce = nch;
                burst_chans_for_record(cfg, trec, dt_rec, nch, &cs, &ce);
                std::fprintf(stderr,
                    "rec file=%d slice=%d r=%d trec=%.6f cs=%d ce=%d\n",
                    i, sl, r, trec, cs, ce);
                if (cs >= ce) continue;

                // uvw at record midpoint (apparent), then rotate to J2000
                const double mjd_utc = rs.mjd_ref + trec / 86400.0;
                const double ha = lmst_rad(mjd_utc) - cfg.ra_app;
                std::vector<UvwSec> uvw_app;
                compute_uvw_seconds(as, ha, cfg.dec_app, uvw_app);

                static bool _diag_once = false;
                if (!_diag_once) {
                    _diag_once = true;
                    std::fprintf(stderr,
                        "uvw-diag: mjd_ref=%.6f trec=%.6f mjd_utc=%.6f "
                        "lmst=%.6f ra_app=%.6f ha=%.6f dec_app=%.6f epoch=%.1f\n",
                        rs.mjd_ref, trec, mjd_utc, lmst_rad(mjd_utc),
                        cfg.ra_app, ha, cfg.dec_app, cfg.epoch);
                    for (int k = 0; k < 6 && k < (int)as.antenna.size(); ++k) {
                        std::fprintf(stderr,
                            "  ant[%d] name=%s bx=%.3f by=%.3f bz=%.3f\n",
                            k, as.antenna[k].name.c_str(),
                            as.antenna[k].bx, as.antenna[k].by, as.antenna[k].bz);
                    }
                    for (int bb = 0; bb < 6 && bb < as.nbase; ++bb) {
                        std::fprintf(stderr,
                            "  base[%d] a0=%d a1=%d uvw_sec=(%.6e, %.6e, %.6e)\n",
                            bb, as.baseline[bb].s0.ant_id,
                            as.baseline[bb].s1.ant_id,
                            uvw_app[bb].u, uvw_app[bb].v, uvw_app[bb].w);
                    }

                    // ---- antenna-table sanity -------------------------------
                    int n_ant_zero = 0, n_ant_used = 0;
                    for (int k = 0; k < (int)as.antenna.size(); ++k) {
                        const auto& a = as.antenna[k];
                        if (a.bx == 0.0 && a.by == 0.0 && a.bz == 0.0) ++n_ant_zero;
                        else ++n_ant_used;
                    }
                    std::fprintf(stderr,
                        "uvw-diag: antennas total=%zu all-zero(bx=by=bz=0)=%d nonzero=%d\n",
                        as.antenna.size(), n_ant_zero, n_ant_used);

                    // ---- baseline ant_id mapping sanity ---------------------
                    int a0min = 1<<30, a0max = -1, a1min = 1<<30, a1max = -1;
                    int n_oob = 0, n_self = 0;
                    for (int bb = 0; bb < as.nbase; ++bb) {
                        const int a0 = as.baseline[bb].s0.ant_id;
                        const int a1 = as.baseline[bb].s1.ant_id;
                        a0min = std::min(a0min, a0); a0max = std::max(a0max, a0);
                        a1min = std::min(a1min, a1); a1max = std::max(a1max, a1);
                        if (a0 < 0 || a0 >= kMaxAnts || a1 < 0 || a1 >= kMaxAnts) ++n_oob;
                        if (a0 == a1) ++n_self;
                    }
                    std::fprintf(stderr,
                        "uvw-diag: nbase=%d a0=[%d,%d] a1=[%d,%d] "
                        "out-of-range(>=%d or <0)=%d self/autocorr(a0==a1)=%d\n",
                        as.nbase, a0min, a0max, a1min, a1max,
                        kMaxAnts, n_oob, n_self);

                    // ---- aggregate uvw magnitudes (seconds) -----------------
                    double umax = 0, vmax = 0, wmax = 0;
                    double umean = 0, vmean = 0;
                    std::size_t n_zero_bl = 0;
                    for (int bb = 0; bb < as.nbase; ++bb) {
                        const double au = std::fabs(uvw_app[bb].u);
                        const double av = std::fabs(uvw_app[bb].v);
                        const double aw = std::fabs(uvw_app[bb].w);
                        if (au == 0.0 && av == 0.0 && aw == 0.0) ++n_zero_bl;
                        umax = std::max(umax, au); vmax = std::max(vmax, av);
                        wmax = std::max(wmax, aw);
                        umean += au; vmean += av;
                    }
                    if (as.nbase > 0) { umean /= as.nbase; vmean /= as.nbase; }
                    // Convert longest baseline to wavelengths at both band edges
                    // so we can see whether u_lambda will actually spread the grid.
                    const double f_lo = std::min(cfg.freq0_hz, cfg.freq1_hz);
                    const double f_hi = std::max(cfg.freq0_hz, cfg.freq1_hz);
                    std::fprintf(stderr,
                        "uvw-diag: |u|sec max=%.6e mean=%.6e  |v|sec max=%.6e mean=%.6e  "
                        "|w|sec max=%.6e  zero-baselines=%zu/%d\n",
                        umax, umean, vmax, vmean, wmax, n_zero_bl, as.nbase);
                    std::fprintf(stderr,
                        "uvw-diag: longest baseline in lambda: u_max=%.3f..%.3f  "
                        "v_max=%.3f..%.3f  (freq %.4g..%.4g Hz)\n",
                        umax * f_lo, umax * f_hi, vmax * f_lo, vmax * f_hi,
                        f_lo, f_hi);
                }

                if (cfg.epoch >= 1999.0) {
                    for (auto& u : uvw_app) {
                        double in[3]  = { u.u, u.v, u.w };
                        double out[3] = {0,0,0};
                        rotate_uvw_to_j2000(mjd_utc, cfg.iatutc, in, out);
                        u.u = out[0]; u.v = out[1]; u.w = out[2];
                    }
                    static bool _diag_once2 = false;
                    if (!_diag_once2) {
                        _diag_once2 = true;
                        for (int bb = 0; bb < 6 && bb < as.nbase; ++bb) {
                            std::fprintf(stderr,
                                "  base[%d] post-j2000 uvw_sec=(%.6e, %.6e, %.6e)\n",
                                bb, uvw_app[bb].u, uvw_app[bb].v, uvw_app[bb].w);
                        }
                    }
                }

                const float* rec_base = rbuf.data() +
                    static_cast<std::ptrdiff_t>(r) * (nch * as.nbase);

                #pragma omp parallel for num_threads(cfg.num_threads)
                for (int b = 0; b < as.nbase; ++b) {
                    const float* row = rec_base +
                        static_cast<std::ptrdiff_t>(b) * nch;
                    // Decode + calibrate + clip
                    std::vector<Vis> ch(ce - cs);
                    for (int c = cs; c < ce; ++c) {
                        Vis v; decode_packed(row + c, v);
                        if (cfg.do_band || cfg.do_base)
                            apply_calib(v, bp, b, c);
                        ch[c - cs] = v;
                    }
                    if (cfg.do_flag) clip_record(ch, cfg.thresh);

                    GridSamples local_d, local_p;
                    local_d.reserve(ch.size());
                    local_p.reserve(ch.size());
                    for (int c = cs; c < ce; ++c) {
                        const double freq_c = cfg.freq0_hz + c * ch_w;
                        push_sample(local_d, uvw_app[b], freq_c, ch[c - cs]);
                        // PSF samples: same uvw, value = (1, 0)*wt
                        Vis w = ch[c - cs]; w.r = 1.0f; w.i = 0.0f;
                        push_sample(local_p, uvw_app[b], freq_c, w);
                    }
                    #pragma omp critical
                    {
                        dirty_samples.u_lambda.insert(dirty_samples.u_lambda.end(),
                            local_d.u_lambda.begin(), local_d.u_lambda.end());
                        dirty_samples.v_lambda.insert(dirty_samples.v_lambda.end(),
                            local_d.v_lambda.begin(), local_d.v_lambda.end());
                        dirty_samples.c.insert(dirty_samples.c.end(),
                            local_d.c.begin(), local_d.c.end());
                        dirty_samples.wt.insert(dirty_samples.wt.end(),
                            local_d.wt.begin(), local_d.wt.end());
                        psf_samples.u_lambda.insert(psf_samples.u_lambda.end(),
                            local_p.u_lambda.begin(), local_p.u_lambda.end());
                        psf_samples.v_lambda.insert(psf_samples.v_lambda.end(),
                            local_p.v_lambda.begin(), local_p.v_lambda.end());
                        psf_samples.c.insert(psf_samples.c.end(),
                            local_p.c.begin(), local_p.c.end());
                        psf_samples.wt.insert(psf_samples.wt.end(),
                            local_p.wt.begin(), local_p.wt.end());
                    }
                }
            }
        }
    }
    close_raw_set(rs);
    std::fprintf(stderr, "gathered %zu visibility samples\n", dirty_samples.size());

    // Map u_λ, v_λ → FINUFFT coords in [-π, π) for the requested cellsize.
    const double cellsize_rad = cfg.cellsize_asec * (M_PI / 180.0 / 3600.0);
    std::fprintf(stderr,
        "imaging: npix=%d cellsize_asec=%.6g cellsize_rad=%.6g  "
        "pre-coords u_lambda sample=[0]%.6g [N/2]%.6g  v=[0]%.6g [N/2]%.6g\n",
        cfg.npix, cfg.cellsize_asec, cellsize_rad,
        dirty_samples.u_lambda.empty()?0.0:dirty_samples.u_lambda[0],
        dirty_samples.u_lambda.empty()?0.0:dirty_samples.u_lambda[dirty_samples.u_lambda.size()/2],
        dirty_samples.v_lambda.empty()?0.0:dirty_samples.v_lambda[0],
        dirty_samples.v_lambda.empty()?0.0:dirty_samples.v_lambda[dirty_samples.v_lambda.size()/2]);
    to_finufft_coords(dirty_samples.u_lambda, dirty_samples.v_lambda, cellsize_rad);
    to_finufft_coords(psf_samples.u_lambda,   psf_samples.v_lambda,   cellsize_rad);
    std::fprintf(stderr,
        "post-coords u sample=[0]%.6g [N/2]%.6g  v=[0]%.6g [N/2]%.6g\n",
        dirty_samples.u_lambda.empty()?0.0:dirty_samples.u_lambda[0],
        dirty_samples.u_lambda.empty()?0.0:dirty_samples.u_lambda[dirty_samples.u_lambda.size()/2],
        dirty_samples.v_lambda.empty()?0.0:dirty_samples.v_lambda[0],
        dirty_samples.v_lambda.empty()?0.0:dirty_samples.v_lambda[dirty_samples.v_lambda.size()/2]);

    DirtyOut dirty, psf;
    if (run_finufft_dirty(dirty_samples, cfg.npix, cfg.npix, 1e-6, dirty) < 0) return -1;
    if (run_finufft_dirty(psf_samples,   cfg.npix, cfg.npix, 1e-6, psf)   < 0) return -1;

    // Scrub any remaining non-finite pixels to 0 (defence in depth — push_sample
    // already drops NaN samples). Count them so we know if something leaked.
    auto scrub = [](std::vector<double>& a, const char* name) {
        std::size_t nbad = 0;
        for (auto& v : a) if (!std::isfinite(v)) { v = 0.0; ++nbad; }
        std::fprintf(stderr, "  scrubbed %zu/%zu non-finite px from %s\n",
                     nbad, a.size(), name);
    };
    scrub(dirty.img, "dirty");
    scrub(psf.img,   "psf");
    auto audit = [](const std::vector<double>& a, const char* name, double sumw) {
        const std::size_t N = a.size();
        const double* p = a.data();
        std::fprintf(stderr,
            "%s ptr=%p size=%zu sample=[0]%.6g [1]%.6g [N/4]%.6g [N/2]%.6g "
            "[3N/4]%.6g [N-1]%.6g\n",
            name, (const void*)p, N,
            N>0?p[0]:0.0, N>1?p[1]:0.0, N>4?p[N/4]:0.0, N>2?p[N/2]:0.0,
            N>4?p[3*N/4]:0.0, N>0?p[N-1]:0.0);
        double dmin =  1e300, dmax = -1e300;
        std::size_t n_nan = 0, n_inf = 0, n_zero = 0, n_neg = 0, n_pos = 0, n_iter = 0;
        for (std::size_t i = 0; i < N; ++i) {
            ++n_iter;
            double v = p[i];
            if (std::isnan(v))      { ++n_nan; continue; }
            if (std::isinf(v))      { ++n_inf; continue; }
            if (v == 0.0)             ++n_zero;
            else if (v < 0.0)         ++n_neg;
            else                      ++n_pos;
            if (v < dmin) dmin = v;
            if (v > dmax) dmax = v;
        }
        std::fprintf(stderr,
            "%s audit: n_iter=%zu nan=%zu inf=%zu zero=%zu neg=%zu pos=%zu "
            "finite_min=%.6g finite_max=%.6g sumw=%.6g\n",
            name, n_iter, n_nan, n_inf, n_zero, n_neg, n_pos,
            dmin, dmax, sumw);
    };
    audit(dirty.img, "dirty", dirty.sumw);
    audit(psf.img,   "psf",   psf.sumw);

    // Normalise PSF to peak=1 (Hogbom assumption).
    double pk = 0.0;
    for (double v : psf.img) if (std::fabs(v) > pk) pk = std::fabs(v);
    if (pk > 0) for (auto& v : psf.img) v /= pk;
    else        std::fprintf(stderr, "WARNING: PSF peak is zero — image will be invalid\n");

    CleanParams cp;
    cp.niter           = cfg.clean_iter;
    cp.gain            = cfg.clean_gain;
    cp.threshold_sigma = cfg.clean_thresh;
    CleanOut co;
    if (hogbom_clean(dirty.img, psf.img, cfg.npix, cfg.npix, cp, co) < 0) return -1;
    std::fprintf(stderr, "clean: %zu components, beam %.2f x %.2f pix\n",
                 co.components.size(), co.bmaj_pix, co.bmin_pix);

    ImageMeta meta;
    meta.n1 = cfg.npix; meta.n2 = cfg.npix;
    double ra_j = cfg.ra_mean, dec_j = cfg.dec_mean;
    if (cfg.ra_mean == 0.0 && cfg.dec_mean == 0.0) {
        apparent_to_j2000(cfg.ra_app, cfg.dec_app, rs.mjd_ref, cfg.iatutc, ra_j, dec_j);
    }
    meta.ra_rad       = ra_j;
    meta.dec_rad      = dec_j;
    meta.cellsize_rad = cellsize_rad;
    meta.freq_hz      = ref_freq;
    meta.object       = cfg.burst.name;
    meta.bmaj_rad     = co.bmaj_pix * cellsize_rad;
    meta.bmin_rad     = co.bmin_pix * cellsize_rad;
    meta.bpa_rad      = co.bpa_rad;

    // Model image: clean components as Jy/pixel deltas (un-convolved), the
    // WSClean "-model" product. Same row-major idx = iy*npix + ix as clean.
    std::vector<double> model(static_cast<std::size_t>(cfg.npix) * cfg.npix, 0.0);
    for (const auto& cc : co.components)
        model[static_cast<std::size_t>(cc.iy) * cfg.npix + cc.ix] += cc.flux;

    return write_fits_images(cfg.image_fits, meta,
                             dirty.img, psf.img, co.residual, co.restored, model);
}

} // namespace pico

// pico_image/src/dump_vis.cpp
// Write pico's calibrated visibilities to a random-groups UVFITS file matching
// svfits's primary HDU layout exactly:
//   PCOUNT=8 params  : UU,VV,WW (sec), BASELINE(=a0*256+a1+257), DATE,DATE,
//                      SOURCE(=1), FREQSEL(=1)
//   data per group   : COMPLEX(re,im,wt) x STOKES(RR,LL)  -> 6 floats
//   one channel per group (one-chan SPOTLIGHT mode); uvw pre-scaled by
//   freq_c/freq0 exactly as svfits make_onechan_group (svsubs.c:1687).
//
// The calibration path mirrors run_pipeline 1:1 (make_bandpass -> apply_calib ->
// conj), so this file IS what pico would image. Dump pre-clip by default so the
// calibrated values are directly comparable to svfits's (whose flagged samples
// carry the same re,im but wt=-1).

#include "pico/dump_vis.hpp"
#include "pico/calib.hpp"
#include "pico/dm.hpp"
#include "pico/uvw.hpp"
#include "pico/half.hpp"
#include "pico/grid.hpp"

#include <fitsio.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <complex>
#include <vector>
#include <unordered_map>
#include <algorithm>

namespace pico {

namespace {

inline void decode_packed(const float* slot, Vis& v) {
    std::uint32_t bits; std::memcpy(&bits, slot, sizeof(bits));
    v.r  = half_to_float(static_cast<std::uint16_t>(bits & 0xFFFFu));
    v.i  = half_to_float(static_cast<std::uint16_t>((bits >> 16) & 0xFFFFu));
    v.wt = 1.0f;
}

long env_long(const char* k, long dflt) {
    const char* s = std::getenv(k);
    return s ? std::strtol(s, nullptr, 10) : dflt;
}

} // namespace

int dump_uvfits(const Config& cfg, const AntSamp& as, RawSet& rs,
                const std::string& path) {
    const long max_groups = env_long("PICO_DUMP_MAXGROUPS", 5000000);
    const int  only_file  = (int)env_long("PICO_DUMP_FILE", -1);
    const int  nch  = rs.channels;
    const double freq0 = cfg.freq0_hz;
    const double ch_w  = (cfg.freq1_hz - cfg.freq0_hz) / cfg.nchan;
    const double dt_rec = rs.t_slice / rs.rec_per_slice;
    const bool conj_vis = (cfg.freq1_hz >= cfg.freq0_hz);

    // Pair each (a0,a1) RR baseline (band 0) with its LL twin (band 1).
    std::unordered_map<int,int> ll_of;     // key=(a0<<8|a1) -> ll baseline idx
    std::vector<int> rr_list;              // rr baseline indices
    for (int b = 0; b < as.nbase; ++b) {
        const auto& e = as.baseline[b];
        if (e.s0.ant_id == e.s1.ant_id) continue;        // skip autocorr
        if (e.drop) continue;                            // user ANTMASK excluded
        const int key = (e.s0.ant_id << 8) | e.s1.ant_id;
        if (e.s0.band == 0 && e.s1.band == 0) rr_list.push_back(b);
        else if (e.s0.band == 1 && e.s1.band == 1) ll_of[key] = b;
    }

    // Collected groups (write GCOUNT-aware header afterwards).
    std::vector<float> P;  // 8 params per group
    std::vector<float> D;  // 6 data floats per group
    P.reserve(2000000); D.reserve(2000000);
    long ng = 0; bool capped = false;
    // aggregate diagnostics (RR): mean raw, mean off_src, mean abp, mean calib
    double Sraw=0, Soff=0, Sabp=0, Scal=0; long Sn=0; long n_offzero=0, n_abpsmall=0;

    std::vector<float> rbuf(rs.rec_per_slice *
                            (std::size_t)(rs.recl / sizeof(float)));

    for (int i = 0; i < cfg.nfile && !capped; ++i) {
        if (only_file >= 0 && i != only_file) continue;
        auto& f = rs.files[i];
        if (f.n_rec <= 0) continue;
        const int slice0 = f.start_rec / rs.rec_per_slice;
        const int slice1 = (f.start_rec + f.n_rec + rs.rec_per_slice - 1) / rs.rec_per_slice;
        Bandpass bp;
        for (int sl = slice0; sl < slice1 && !capped; ++sl) {
            if (read_slice(rs, i, sl, rbuf.data()) < 0) break;
            if (cfg.do_band || cfg.do_base)
                make_bandpass(cfg, rs, as, rbuf.data(), i, sl, bp);
            const double t_slice_start = f.t_start + sl * rs.slice_interval;
            for (int r = 0; r < rs.rec_per_slice && !capped; ++r) {
                const int grec = sl * rs.rec_per_slice + r;
                if (grec < f.start_rec || grec >= f.start_rec + f.n_rec) continue;
                const double trec = t_slice_start + r * dt_rec;
                int cs = 0, ce = nch;
                burst_chans_for_record(cfg, trec, dt_rec, nch, &cs, &ce);
                if (cs >= ce) continue;

                const double mjd_utc = rs.mjd_ref + trec / 86400.0;
                const double ha = lmst_rad(mjd_utc) - cfg.ra_app;
                std::vector<UvwSec> uvw;
                compute_uvw_seconds(as, ha, cfg.dec_app, uvw);
                if (cfg.epoch >= 1999.0)
                    for (auto& u : uvw) {
                        double in[3]={u.u,u.v,u.w}, out[3]={0,0,0};
                        rotate_uvw_to_j2000(mjd_utc, cfg.iatutc, in, out);
                        u.u=out[0]; u.v=out[1]; u.w=out[2];
                    }
                const double JD = rs.mjd_ref + 2400000.5 + trec / 86400.0;
                const float date1 = (float)std::floor(JD);
                const float date2 = (float)(JD - std::floor(JD) + cfg.iatutc / 86400.0);

                const float* rec_base = rbuf.data() + (std::ptrdiff_t)r * (nch * as.nbase);

                for (int rr : rr_list) {
                    const auto& e = as.baseline[rr];
                    const int key = (e.s0.ant_id << 8) | e.s1.ant_id;
                    auto it = ll_of.find(key);
                    if (it == ll_of.end()) continue;
                    const int ll = it->second;
                    const float* rowRR = rec_base + (std::ptrdiff_t)rr * nch;
                    const float* rowLL = rec_base + (std::ptrdiff_t)ll * nch;
                    const float blcode = (float)(e.s0.ant_id * 256 + e.s1.ant_id + 257);

                    for (int c = cs; c < ce; ++c) {
                        Vis vrr, vll;
                        decode_packed(rowRR + c, vrr);
                        decode_packed(rowLL + c, vll);
                        // --- diagnostic: capture raw, off_src, abp for RR ---
                        const float raw_re = vrr.r;
                        float off_re = 0.f, abp_v = 1.f;
                        if ((cfg.do_band || cfg.do_base) && rr < bp.n_base
                            && !bp.off_src.empty()) {
                            off_re = bp.off_src[rr][c].real();
                            abp_v  = bp.abp[rr][c];
                        }
                        if (cfg.do_band || cfg.do_base) {
                            apply_calib(vrr, bp, rr, c);
                            apply_calib(vll, bp, ll, c);
                        }
                        Sraw += raw_re; Soff += off_re; Sabp += abp_v; Scal += vrr.r; ++Sn;
                        if (off_re == 0.f) ++n_offzero;
                        if (abp_v > 0.f && abp_v < 0.1f) ++n_abpsmall;
                        if (conj_vis) { vrr.i = -vrr.i; vll.i = -vll.i; }
                        const double fc = freq0 + c * ch_w;
                        const double sc = fc / freq0;
                        P.push_back((float)(uvw[rr].u * sc));
                        P.push_back((float)(uvw[rr].v * sc));
                        P.push_back((float)(uvw[rr].w * sc));
                        P.push_back(blcode);
                        P.push_back(date1);
                        P.push_back(date2);
                        P.push_back(1.0f);          // SOURCE
                        P.push_back(1.0f);          // FREQSEL
                        D.push_back(vrr.r); D.push_back(vrr.i); D.push_back(vrr.wt);
                        D.push_back(vll.r); D.push_back(vll.i); D.push_back(vll.wt);
                        if (++ng >= max_groups) { capped = true; break; }
                    }
                    if (capped) break;
                }
            }
        }
    }

    if (ng == 0) { std::fprintf(stderr, "dump_uvfits: no groups collected\n"); return -1; }

    if (Sn > 0)
        std::fprintf(stderr,
            "dump_uvfits DIAG (RR, n=%ld): <raw_re>=%.4f  <off_src_re>=%.4f  "
            "<abp>=%.4f  <calib_re>=%.4f   (raw-off)=%.4f  off==0: %.2f%%  "
            "abp<0.1: %.2f%%\n",
            Sn, Sraw/Sn, Soff/Sn, Sabp/Sn, Scal/Sn, (Sraw-Soff)/Sn,
            100.0*n_offzero/Sn, 100.0*n_abpsmall/Sn);

    // ---- write the UVFITS primary HDU (random groups) -----------------------
    fitsfile* fp = nullptr; int st = 0;
    std::string p = "!" + path;                 // '!' = overwrite
    ffinit(&fp, p.c_str(), &st);
    int   bitpix = -32, naxis = 7;
    long  naxes[7] = {0, 3, 2, 1, 1, 1, 1};      // 0, COMPLEX, STOKES, FREQ, IF, RA, DEC
    long  pcount = 8, gcount = ng;
    ffphpr(fp, 1 /*simple*/, bitpix, naxis, naxes, pcount, gcount, 1 /*extend*/, &st);

    const char* ptype[8] = {"UU---SIN","VV---SIN","WW---SIN","BASELINE",
                            "DATE","DATE","SOURCE","FREQSEL"};
    for (int j = 0; j < 8; ++j) {
        char k[16]; float one = 1.0f, zero = 0.0f;
        std::snprintf(k, sizeof(k), "PTYPE%d", j+1); ffpkys(fp, k, (char*)ptype[j], "", &st);
        std::snprintf(k, sizeof(k), "PSCAL%d", j+1); ffpkye(fp, k, one,  6, "", &st);
        std::snprintf(k, sizeof(k), "PZERO%d", j+1); ffpkye(fp, k, zero, 6, "", &st);
    }
    // axis descriptors (CTYPE2..7); FREQ axis carries the reference frequency.
    const char* ctype[6] = {"COMPLEX","STOKES","FREQ","IF","RA","DEC"};
    float crval[6] = {1.0f, -1.0f /*RR*/, (float)freq0, 1.0f, (float)cfg.ra_mean, (float)cfg.dec_mean};
    float cdelt[6] = {1.0f, -1.0f /*RR,LL*/, (float)ch_w, 1.0f, 1.0f, 1.0f};
    float crpix[6] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
    for (int j = 0; j < 6; ++j) {
        char k[16];
        std::snprintf(k, sizeof(k), "CTYPE%d", j+2); ffpkys(fp, k, (char*)ctype[j], "", &st);
        std::snprintf(k, sizeof(k), "CRVAL%d", j+2); ffpkye(fp, k, crval[j], 8, "", &st);
        std::snprintf(k, sizeof(k), "CDELT%d", j+2); ffpkye(fp, k, cdelt[j], 8, "", &st);
        std::snprintf(k, sizeof(k), "CRPIX%d", j+2); ffpkye(fp, k, crpix[j], 8, "", &st);
    }
    ffpkys(fp, "OBJECT",   (char*)cfg.burst.name.c_str(), "", &st);
    ffpkys(fp, "TELESCOP", (char*)"GMRT", "", &st);
    ffpkys(fp, "INSTRUME", (char*)"GMRT", "", &st);
    ffpkys(fp, "BUNIT",    (char*)"UNCALIB", "", &st);
    { double v=1.0; ffpkyd(fp,"BSCALE",v,11,"",&st); v=0.0; ffpkyd(fp,"BZERO",v,11,"",&st); }
    { double ep=cfg.epoch>0?cfg.epoch:2000.0; ffpkyd(fp,"EPOCH",ep,9,"",&st); }

    // write the groups: params via ffpgpe, data via ffppre (both swap to BE).
    for (long g = 0; g < ng && !st; ++g) {
        ffpgpe(fp, g+1, 1, 8, &P[(std::size_t)g*8], &st);
        ffppre(fp, g+1, 1, 6, &D[(std::size_t)g*6], &st);
    }
    ffclos(fp, &st);
    if (st) { char e[80]; ffgerr(st, e);
        std::fprintf(stderr, "dump_uvfits: cfitsio status=%d (%s)\n", st, e); return -1; }

    std::fprintf(stderr,
        "dump_uvfits: wrote %ld groups (%s%s) to %s  [stokes=2 RR,LL; %d chan/group]\n",
        ng, capped ? "CAPPED at PICO_DUMP_MAXGROUPS=" : "all",
        capped ? std::to_string(max_groups).c_str() : "", path.c_str(), 1);
    return 0;
}

} // namespace pico

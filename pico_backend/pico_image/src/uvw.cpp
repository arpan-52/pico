// pico_image/src/uvw.cpp
// uvw geometry + J2000 rotation. Ports svfits/svsubs.c svgetUvw().
// J2000 rotation uses SuperNOVAS tod_to_j2000 instead of SLA/NOVAS-legacy.

#include "pico/uvw.hpp"
#include <cmath>

extern "C" {
#include "novas.h"
}

namespace pico {

constexpr double C_LIGHT = 2.99792458e8;     // m/s
constexpr double GMRT_LON_RAD = 1.292439;    // 74.0497° E (svfits canonical)

double lmst_rad(double mjd_utc) {
    // Greenwich Mean Sidereal Time, simple IAU-1982 series.
    // Good to ~1 arcsec — fine for uvw geometry.
    const double T = (mjd_utc - 51544.5) / 36525.0;
    double gmst_sec = 24110.54841 + 8640184.812866 * T
                    + 0.093104 * T * T - 6.2e-6 * T * T * T;
    // add UT1 (≈ UTC) contribution
    const double ut1_day = mjd_utc - std::floor(mjd_utc);
    gmst_sec += ut1_day * 86400.0 * 1.00273790935;
    double gmst_rad = std::fmod(gmst_sec * (2.0 * M_PI / 86400.0), 2.0 * M_PI);
    if (gmst_rad < 0) gmst_rad += 2.0 * M_PI;
    double lmst = gmst_rad + GMRT_LON_RAD;
    lmst = std::fmod(lmst, 2.0 * M_PI);
    if (lmst < 0) lmst += 2.0 * M_PI;
    return lmst;
}

void compute_uvw_seconds(const AntSamp& as,
                         double ha_rad, double dec_rad,
                         std::vector<UvwSec>& out_uvw_sec) {
    out_uvw_sec.assign(as.nbase, UvwSec{0, 0, 0});
    const double ch = std::cos(ha_rad),  sh = std::sin(ha_rad);
    const double cd = std::cos(dec_rad), sd = std::sin(dec_rad);

    // Per-antenna uvw (in seconds), then baseline uvw = ant1 - ant0
    std::vector<UvwSec> per_ant(kMaxAnts, {0, 0, 0});
    for (int k = 0; k < kMaxAnts; ++k) {
        const double bx = as.antenna[k].bx;
        const double by = as.antenna[k].by;
        const double bz = as.antenna[k].bz;
        const double bxch = bx * ch;
        const double bysh = by * sh;
        per_ant[k].u = ( bx * sh + by * ch) / C_LIGHT;
        per_ant[k].v = ( bz * cd - sd * (bxch - bysh)) / C_LIGHT;
        per_ant[k].w = ( cd * (bxch - bysh) + sd * bz) / C_LIGHT;
    }
    for (int b = 0; b < as.nbase; ++b) {
        const int a0 = as.baseline[b].s0.ant_id;
        const int a1 = as.baseline[b].s1.ant_id;
        out_uvw_sec[b].u = per_ant[a1].u - per_ant[a0].u;
        out_uvw_sec[b].v = per_ant[a1].v - per_ant[a0].v;
        out_uvw_sec[b].w = per_ant[a1].w - per_ant[a0].w;
    }
}

static double mjd_utc_to_jd_tt(double mjd_utc, double iatutc /*s*/) {
    // JD = MJD + 2400000.5; TT = TAI + 32.184s ; TAI = UTC + iatutc
    const double tt_offset_sec = iatutc + 32.184;
    return mjd_utc + 2400000.5 + tt_offset_sec / 86400.0;
}

void apparent_to_j2000(double ra_app, double dec_app, double mjd_utc,
                       double iatutc, double& ra_j2000, double& dec_j2000) {
    // Build unit vector in apparent (TOD) frame, rotate to J2000 (mean of
    // J2000.0). NOTE: this ignores aberration & light deflection (~20"); for
    // sub-arcsec astrometry use SuperNOVAS frame-based API.
    const double in[3] = { std::cos(dec_app) * std::cos(ra_app),
                           std::cos(dec_app) * std::sin(ra_app),
                           std::sin(dec_app) };
    double out[3] = {0};
    const double jd_tt = mjd_utc_to_jd_tt(mjd_utc, iatutc);
    tod_to_j2000(jd_tt, NOVAS_REDUCED_ACCURACY, in, out);
    ra_j2000  = std::atan2(out[1], out[0]);
    if (ra_j2000 < 0) ra_j2000 += 2.0 * M_PI;
    dec_j2000 = std::asin(out[2]);
}

void rotate_uvw_to_j2000(double mjd_utc, double iatutc,
                         double uvw_in[3], double uvw_out[3]) {
    const double jd_tt = mjd_utc_to_jd_tt(mjd_utc, iatutc);
    tod_to_j2000(jd_tt, NOVAS_REDUCED_ACCURACY, uvw_in, uvw_out);
}

} // namespace pico

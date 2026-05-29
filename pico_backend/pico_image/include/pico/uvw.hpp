#pragma once
// uvw geometry + J2000 rotation. Ports svfits's svgetUvw() and the J2000
// rotation that was done via sla_prenut_vis/novas_prenut_vis. We use
// SuperNOVAS for precession+nutation.

#include "pico/types.hpp"
#include "pico/config.hpp"
#include "pico/antsamp.hpp"
#include <vector>

namespace pico {

// Compute uvw (seconds) for every baseline at the given MJD (UTC).
// `ha_rad` is local hour angle of phase centre, `dec_rad` is declination.
// On entry the AntSamp baseline list is assumed already built.
void compute_uvw_seconds(const AntSamp& as,
                         double ha_rad, double dec_rad,
                         std::vector<UvwSec>& out_uvw_sec);

// Local Mean Sidereal Time (rad) at mjd (UT1≈UTC for our precision).
// GMRT longitude is baked in for now.
double lmst_rad(double mjd_utc);

// Apparent (ra_app, dec_app) at given MJD → J2000 (rc, dc). Wraps SuperNOVAS
// tod_to_j2000() on the unit direction vector. iatutc lets us go UTC→TT.
void apparent_to_j2000(double ra_app, double dec_app, double mjd_utc,
                       double iatutc, double& ra_j2000, double& dec_j2000);

// Rotate a uvw vector (seconds) from apparent frame to J2000 at given MJD.
// Composes precession+nutation via SuperNOVAS tod_to_j2000.
void rotate_uvw_to_j2000(double mjd_utc, double iatutc,
                         double uvw_in[3], double uvw_out[3]);

} // namespace pico

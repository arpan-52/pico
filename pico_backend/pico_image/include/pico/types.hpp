#pragma once
#include <cstdint>
#include <complex>
#include <vector>
#include <string>

namespace pico {

using Complex = std::complex<float>;

struct UvwSec { double u, v, w; };        // uvw in seconds (geometry only, freq-independent)
struct Vis    { float r, i, wt; };        // one (channel,baseline,stokes) sample

// Per-baseline geometry slot derived from antsamp.hdr (svfits BaseUvwType analogue)
struct BaselineGeom {
    int ant0, ant1;
    int band0, band1;          // 0=USB-130(RR), 1=USB-175(LL)
    double bx, by, bz;         // ant1 - ant0, metres
    double u, v, w;            // computed per integration, in seconds
};

// All visibilities for one record (one integration time), one baseline,
// one stokes, all channels.
struct Record {
    double mjd;                // record midpoint MJD (UTC)
    UvwSec uvw;                // current J2000 (or apparent) uvw, seconds
    int    baseline_idx;       // index into BaselineGeom table
    int    stokes;             // 0=RR, 1=LL
    std::vector<Vis> ch;       // size = n_channels_kept
    int    chan_start;         // first kept channel (for DM tracking)
    int    chan_end;           // one past last kept channel
};

} // namespace pico

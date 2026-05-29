#pragma once
// Parser for the GMRT antsamp.hdr file. Ports svfits's get_antenna / get_bands /
// get_sampler / get_baseline / init_corr (only the parts needed for SPOTLIGHT
// RR/LL raw visibility decode).

#include "pico/types.hpp"
#include <array>
#include <string>
#include <vector>
#include <cstdint>

namespace pico {

constexpr int kMaxAnts  = 32;
constexpr int kMaxBands = 4;
constexpr int kMaxSamps = 64;

struct Antenna {
    std::string name;            // e.g. "C00", "E02"
    double bx = 0, by = 0, bz = 0;        // metres
    std::array<double, kMaxBands> d0{};   // delay offsets (unused for imaging)
    std::array<int, kMaxBands> samp_id{}; // sampler index per band, or kMaxSamps if absent
};

struct Sampler {
    int ant_id = kMaxAnts;
    int band   = kMaxBands;     // 0=USB-130 (RR), 1=USB-175 (LL)
};

// Single baseline as ordered by svfits's get_baseline (all RR baselines first,
// then all LL). Matches the per-record visibility layout in the raw files.
struct BaselineEntry {
    Sampler s0;
    Sampler s1;
};

struct AntSamp {
    std::vector<Antenna>       antenna;     // size 32
    std::vector<Sampler>       sampler;     // size = nsamp (typically 64)
    std::vector<BaselineEntry> baseline;    // size = nsamp*(nsamp+1)/2  per stokes
    uint32_t                   antmask  = 0;
    uint32_t                   bandmask = 0;
    int                        nsamp    = 0;
    int                        nbase    = 0;   // total (== RR + LL)
    int                        stokes   = 2;
};

int load_antsamp(const std::string& path, AntSamp& out);

// Apply user antmask (from svfits_par) on top of the antsamp-derived antmask,
// rebuild the baseline list. Output baselines is the set passed to the imager.
int rebuild_baselines(AntSamp& as, uint32_t user_antmask);

} // namespace pico

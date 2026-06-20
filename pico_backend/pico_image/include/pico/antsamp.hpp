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
//
// IMPORTANT: the baseline list is built from the FULL hardware sampler set
// (every antenna that has a sampler), so its size and ordering match the raw
// file's on-disk record layout exactly (= svfits corr->daspar.baselines, the
// record stride). The per-record byte offset of baseline b is therefore just
// b*channels*sizeof(float) — contiguous. The user ANTMASK (svfits user->antmask)
// does NOT shrink this list; instead baselines whose antennas are excluded by
// the user mask are marked `drop` and skipped at imaging/calibration (this is
// svfits's init_vispar selection, mapped back onto the contiguous input list).
struct BaselineEntry {
    Sampler s0;
    Sampler s1;
    bool    drop = false;   // antenna excluded by user ANTMASK → not imaged
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

// Build the FULL hardware baseline list (size/order = on-disk record layout,
// = record stride) and mark baselines excluded by the user ANTMASK as `drop`.
// nbase is the on-disk baseline count (e.g. 1056), NOT the masked subset.
int rebuild_baselines(AntSamp& as, uint32_t user_antmask);

} // namespace pico

#pragma once
// Per-file raw visibility reader. Replaces svfits's open_file + read_slice
// pattern that fopens every slice — we keep one FD per file open for the
// life of the run and use pread() (thread-safe positional read) so multiple
// slices can be pulled in parallel without locking.

#include "pico/config.hpp"
#include "pico/antsamp.hpp"
#include <string>
#include <vector>
#include <cstdint>

namespace pico {

struct RawFile {
    int    fd        = -1;
    int    file_idx  = -1;      // 0..15
    std::string path;
    // From per-file metadata (currently HAVE_IDX=1 path):
    double t_start   = 0.0;     // first record time (seconds since mjd_ref)
    double b_start   = 0.0;     // first record-with-burst time
    int    start_rec = 0;
    int    n_rec     = 0;
};

struct RawSet {
    std::vector<RawFile>  files;   // size = nfile (≤16)
    int    rec_per_slice = 0;
    double t_slice       = 0.0;   // duration of one slice (s)
    double slice_interval = 0.0;  // gap between slices on the SAME file
    double mjd_ref       = 0.0;
    int    channels      = 4096;
    int    baselines     = 0;     // total (RR+LL) baselines on disk

    // bytes per (one record == one integration, all chans, all baselines):
    // channels * baselines * sizeof(float) — half-float pair packed into float
    std::size_t recl     = 0;
    std::size_t timeval_size = 16;   // sizeof(struct timeval) — 16 on x86_64
};

int open_raw_set(const Config& cfg, const AntSamp& as, RawSet& out);
void close_raw_set(RawSet& rs);

// Read raw slice `slice` from file `idx` into `dst`. dst must have rs.recl *
// rs.rec_per_slice bytes. Thread-safe; uses pread().
int read_slice(const RawSet& rs, int idx, int slice, void* dst);

// Total slices available in a given file (computed from filesize).
int n_slices(const RawSet& rs, int idx);

} // namespace pico

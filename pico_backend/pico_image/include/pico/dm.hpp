#pragma once
// DM-aware time / frequency channel selection per record. Ports svfits's
// get_chan_num() + the burst-slice / start_rec / n_rec setup.

#include "pico/config.hpp"
#include "pico/antsamp.hpp"
#include "pico/raw_io.hpp"

namespace pico {

// Cold-plasma dispersion delay (s) at freq_hz, relative to infinite frequency.
// t = K * DM / freq_GHz^2 ;  K = 4.148808e-3 (s · GHz^2 / pc/cm^3) ← svfits convention
inline double dm_delay_sec(double DM, double freq_hz) {
    const double f_ghz = freq_hz * 1e-9;
    return 4.148808e-3 * DM / (f_ghz * f_ghz);
}

// For a record at time trec (s since mjd_ref) of duration integ_sec, determine
// the [start,end) channel range expected to contain the burst signal.
// Returns 0 on success; if the burst is not in this record, *cs = *ce = nchan_full.
int burst_chans_for_record(const Config& cfg, double trec, double integ_sec,
                           int nchan_full, int* cs, int* ce);

// Compute per-file start_rec / n_rec (which records to read, where the burst
// is) — ports svfits's get_rec_num. nf is files count, t_slice in seconds.
int compute_record_ranges(const Config& cfg, RawSet& rs);

} // namespace pico

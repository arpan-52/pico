#pragma once
// Top-level pico_image pipeline driver. mirrors svfits's copy_burst() flow:
//
//   for file f in files_with_burst:
//     for slice s overlapping burst in f:
//       read_slice(f, s)
//       make_bandpass(off-source records)
//       for record r in slice:
//         (cs, ce) = burst_chans_for_record(...)
//         for each baseline:
//           apply_calib(vis); clip; uvw = svgetUvw(t); rotate→J2000;
//           push_sample → GridSamples
//
//   run_finufft → dirty + psf
//   hogbom_clean
//   write_fits_images → <base>-{dirty,psf,image,residual,model}.fits

#include "pico/config.hpp"
#include "pico/antsamp.hpp"
#include "pico/raw_io.hpp"
#include "pico/grid.hpp"

namespace pico {

int run_pipeline(const Config& cfg);

} // namespace pico

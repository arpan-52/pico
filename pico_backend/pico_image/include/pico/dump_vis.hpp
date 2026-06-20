#pragma once
// Debug: dump pico's calibrated visibilities to a random-groups UVFITS file in
// the SAME format svfits writes (UvwParType 8 params + Cmplx3Type[stokes]), so
// the two can be imaged in CASA identically and diffed visibility-for-visibility.
//
// Enabled from run_pipeline when env PICO_DUMP_UVFITS=<path> is set.
// Optional caps: PICO_DUMP_MAXGROUPS (default 5,000,000), PICO_DUMP_FILE (only
// that file index), PICO_DUMP_CLIP=1 (apply the per-slice MAD clip → wt=-1).

#include "pico/config.hpp"
#include "pico/antsamp.hpp"
#include "pico/raw_io.hpp"
#include <string>

namespace pico {

int dump_uvfits(const Config& cfg, const AntSamp& as, RawSet& rs,
                const std::string& path);

} // namespace pico

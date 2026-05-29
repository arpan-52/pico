#pragma once
// CASA/WSClean-friendly FITS image writer. We write one file per product —
// <base>-dirty.fits, <base>-psf.fits, <base>-image.fits (restored),
// <base>-residual.fits, <base>-model.fits — mirroring WSClean's layout.
// Coordinates are RA---SIN / DEC--SIN tangent-plane at the phase centre
// (J2000); RA axis runs east-to-left (CDELT1 < 0, the CASA/WSClean
// convention). Frequency axis is the MFS reference frequency.

#include "pico/config.hpp"
#include <string>
#include <vector>

namespace pico {

struct ImageMeta {
    int    n1, n2;
    double ra_rad, dec_rad;        // J2000 phase centre
    double cellsize_rad;           // radians / pix (both axes)
    double freq_hz;                // MFS ref freq for FREQ axis
    double bmaj_rad = 0;           // restoring beam (0 ⇒ omit BMAJ/BMIN/BPA)
    double bmin_rad = 0;
    double bpa_rad  = 0;
    std::string telescope = "GMRT";
    std::string object;
    std::string date_obs;          // ISO 8601
};

// Write all five products to separate files derived from `base` (a trailing
// ".fits"/".FITS" is stripped before appending the WSClean suffixes).
int write_fits_images(const std::string& base, const ImageMeta& meta,
                      const std::vector<double>& dirty,
                      const std::vector<double>& psf,
                      const std::vector<double>& residual,
                      const std::vector<double>& restored,
                      const std::vector<double>& model);

} // namespace pico

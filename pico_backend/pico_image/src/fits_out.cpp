// pico_image/src/fits_out.cpp
// WSClean-style FITS writer: one single-image file per product. Coordinates
// are J2000 RA---SIN / DEC--SIN tangent-plane, RA running east-to-left
// (CDELT1 < 0), frequency third axis at the MFS reference frequency. The
// restoring beam (BMAJ/BMIN/BPA) is written when meta carries one and the
// product is the restored/residual image.

#include "pico/fits_out.hpp"
#include <fitsio.h>
#include <cstdio>
#include <vector>
#include <cstring>
#include <cctype>
#include <string>

namespace pico {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kRad2Deg = 180.0 / kPi;

// Write `img` as the primary HDU of a freshly created file at `path`, with
// full WCS. `write_beam` adds BMAJ/BMIN/BPA from meta (restored/residual only).
int write_image_file(const std::string& path, const ImageMeta& m,
                     const std::vector<double>& img, bool write_beam) {
    std::remove(path.c_str());
    fitsfile* fp = nullptr;
    int status = 0;
    fits_create_file(&fp, path.c_str(), &status);
    if (status) { fits_report_error(stderr, status); return -1; }

    long naxes[4] = { m.n1, m.n2, 1, 1 };
    fits_create_img(fp, DOUBLE_IMG, 4, naxes, &status);

    const double cdelt  = m.cellsize_rad * kRad2Deg;
    const double crpix1 = m.n1 / 2.0 + 1.0;
    const double crpix2 = m.n2 / 2.0 + 1.0;
    const double crval1 = m.ra_rad  * kRad2Deg;
    const double crval2 = m.dec_rad * kRad2Deg;
    const double cdelt_neg = -cdelt;   // RA increases east-to-left
    const double one = 1.0;

    fits_update_key(fp, TSTRING, "CTYPE1", const_cast<char*>("RA---SIN"), nullptr, &status);
    fits_update_key(fp, TSTRING, "CTYPE2", const_cast<char*>("DEC--SIN"), nullptr, &status);
    fits_update_key(fp, TSTRING, "CTYPE3", const_cast<char*>("FREQ"),     nullptr, &status);
    fits_update_key(fp, TSTRING, "CTYPE4", const_cast<char*>("STOKES"),   nullptr, &status);
    fits_update_key(fp, TSTRING, "CUNIT1", const_cast<char*>("deg"), nullptr, &status);
    fits_update_key(fp, TSTRING, "CUNIT2", const_cast<char*>("deg"), nullptr, &status);
    fits_update_key(fp, TSTRING, "CUNIT3", const_cast<char*>("Hz"),  nullptr, &status);
    fits_update_key(fp, TDOUBLE, "CRVAL1", const_cast<double*>(&crval1), nullptr, &status);
    fits_update_key(fp, TDOUBLE, "CRVAL2", const_cast<double*>(&crval2), nullptr, &status);
    fits_update_key(fp, TDOUBLE, "CRVAL3", const_cast<double*>(&m.freq_hz), nullptr, &status);
    fits_update_key(fp, TDOUBLE, "CRVAL4", const_cast<double*>(&one), nullptr, &status);
    fits_update_key(fp, TDOUBLE, "CRPIX1", const_cast<double*>(&crpix1), nullptr, &status);
    fits_update_key(fp, TDOUBLE, "CRPIX2", const_cast<double*>(&crpix2), nullptr, &status);
    fits_update_key(fp, TDOUBLE, "CRPIX3", const_cast<double*>(&one), nullptr, &status);
    fits_update_key(fp, TDOUBLE, "CRPIX4", const_cast<double*>(&one), nullptr, &status);
    fits_update_key(fp, TDOUBLE, "CDELT1", const_cast<double*>(&cdelt_neg), nullptr, &status);
    fits_update_key(fp, TDOUBLE, "CDELT2", const_cast<double*>(&cdelt), nullptr, &status);
    fits_update_key(fp, TDOUBLE, "CDELT3", const_cast<double*>(&one), nullptr, &status);
    fits_update_key(fp, TDOUBLE, "CDELT4", const_cast<double*>(&one), nullptr, &status);
    fits_update_key(fp, TSTRING, "TELESCOP", const_cast<char*>(m.telescope.c_str()), nullptr, &status);
    if (!m.object.empty())
        fits_update_key(fp, TSTRING, "OBJECT", const_cast<char*>(m.object.c_str()), nullptr, &status);
    if (!m.date_obs.empty())
        fits_update_key(fp, TSTRING, "DATE-OBS", const_cast<char*>(m.date_obs.c_str()), nullptr, &status);

    if (write_beam && m.bmaj_rad > 0.0) {
        double bmaj = m.bmaj_rad * kRad2Deg;
        double bmin = m.bmin_rad * kRad2Deg;
        double bpa  = m.bpa_rad  * kRad2Deg;
        fits_update_key(fp, TDOUBLE, "BMAJ", &bmaj, nullptr, &status);
        fits_update_key(fp, TDOUBLE, "BMIN", &bmin, nullptr, &status);
        fits_update_key(fp, TDOUBLE, "BPA",  &bpa,  nullptr, &status);
    }

    const long npix = static_cast<long>(m.n1) * m.n2;
    long fpixel = 1;
    fits_write_img(fp, TDOUBLE, fpixel, npix, const_cast<double*>(img.data()), &status);

    fits_close_file(fp, &status);
    if (status) { fits_report_error(stderr, status); return -1; }
    return 0;
}

} // namespace

int write_fits_images(const std::string& base, const ImageMeta& meta,
                      const std::vector<double>& dirty,
                      const std::vector<double>& psf,
                      const std::vector<double>& residual,
                      const std::vector<double>& restored,
                      const std::vector<double>& model) {
    // Strip a trailing .fits/.FITS so suffixes append cleanly.
    std::string stem = base;
    const std::size_t dot = stem.rfind('.');
    if (dot != std::string::npos) {
        std::string ext = stem.substr(dot);
        for (auto& c : ext) c = static_cast<char>(::tolower(c));
        if (ext == ".fits") stem.erase(dot);
    }

    int rc = 0;
    rc |= write_image_file(stem + "-dirty.fits",    meta, dirty,    false);
    rc |= write_image_file(stem + "-psf.fits",      meta, psf,      false);
    rc |= write_image_file(stem + "-image.fits",    meta, restored, true);
    rc |= write_image_file(stem + "-residual.fits", meta, residual, true);
    rc |= write_image_file(stem + "-model.fits",    meta, model,    false);
    return rc ? -1 : 0;
}

} // namespace pico

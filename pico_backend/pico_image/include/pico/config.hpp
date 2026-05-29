#pragma once
// pico_image configuration. Ports svfits's SvSelectionType + svuserInp/init_user,
// extended with imaging keys (NPIX, CELLSIZE_ASEC, ...).
//
// All keys parsed from svfits_par.txt; unknown keys are ignored (svfits behavior).

#include <string>
#include <vector>
#include <cstdint>

namespace pico {

struct BurstPar {
    std::string name;
    double mjd          = 0.0;
    double t            = 0.0;     // arrival time (s, IST)
    double dt           = 0.0;
    double width        = 0.0;     // observed width (s, dispersion-dominated)
    double int_wd       = 0.0;     // intrinsic width (s)
    double DM           = 0.0;
    double dDM          = 0.0;
    double freq         = 0.0;     // ref frequency for t (Hz)
    int    bm_id        = 0;
    double ra_app       = 0.0;
    double dec_app      = 0.0;
};

struct Config {
    // ---- I/O (svfits par keys) ----
    int                       nfile     = 16;
    int                       have_idx  = 1;
    std::string               path;
    std::vector<std::string>  input;          // per-file names
    std::string               fits_uv;        // FITS UV (svfits output) — unused here, kept for compat
    int                       all_chan  = 0;
    int                       all_data  = 0;
    double                    lta       = 1.0;
    int                       n_lta     = -1;
    int                       nchav     = 16;

    // ---- Observing (svfits par keys) ----
    double                    obs_mjd      = 0.0;
    double                    freq0_hz     = 0.0;   // FREQ_SET first
    double                    freq1_hz     = 0.0;
    int                       nchan        = 4096;
    uint32_t                  antmask      = 0xFFFFFFFFu;
    double                    epoch        = 2000.0;
    int                       stokes_type  = 0;     // 0=RL, 1=XY
    double                    iatutc       = 37.0;
    double                    ra_mean      = 0.0;
    double                    dec_mean     = 0.0;
    double                    ra_app       = 0.0;
    double                    dec_app      = 0.0;

    // ---- Burst (svfits par keys) ----
    BurstPar                  burst;
    int                       update_burst = 1;

    // ---- Processing (svfits par keys) ----
    int                       do_band      = 1;
    int                       do_base      = 1;
    int                       do_flag      = 1;
    float                     thresh       = 20.0f;
    int                       post_corr    = 0;
    int                       drop_csq     = 0;
    int                       recentre     = 0;
    int                       num_threads  = 16;

    // ---- pico_image imaging extensions ----
    int                       npix         = 2048;     // image side
    double                    cellsize_asec = 1.0;     // arcsec / pix
    std::string               weighting    = "natural"; // natural | uniform
    int                       clean_iter   = 1000;
    float                     clean_gain   = 0.1f;
    float                     clean_thresh = 3.0f;     // CLEAN_THRESH, in sigma
    std::string               image_fits   = "PICO_IMAGE.FITS";
    double                    ref_freq_hz  = 0.0;      // 0 ⇒ midband

    // Resolved at load time:
    std::string               ant_hdr_path;
};

// Parse svfits_par.txt-style file. Lines starting with '*' are comments,
// '!' starts a trailing comment. KEYWORD VALUE per line.
// Returns 0 on success, <0 on error.
int load_config(const std::string& path, Config& out);

// Compute reference frequency if ref_freq_hz==0 (default: band midpoint).
double resolve_ref_freq(const Config& c);

} // namespace pico

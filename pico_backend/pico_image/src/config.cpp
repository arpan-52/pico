// pico_image/src/config.cpp
// KEYWORD VALUE parser for svfits_par.txt + extended imaging keys.
// Ports the relevant parts of svfits/svsubs.c svuserInp(). Lines starting
// with '*' are comments; '!' starts a trailing comment.

#include "pico/config.hpp"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cctype>
#include <string>
#include <sstream>

namespace pico {

namespace {

std::string strip_trailing_comment(std::string s) {
    auto p = s.find('!');
    if (p != std::string::npos) s.resize(p);
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) s.pop_back();
    std::size_t i = 0;
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    return s.substr(i);
}

std::vector<std::string> split_csv(const std::string& s) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : s) {
        if (c == ',') { out.push_back(cur); cur.clear(); }
        else if (!std::isspace(static_cast<unsigned char>(c))) cur.push_back(c);
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
}

} // namespace

int load_config(const std::string& path, Config& out) {
    FILE* fp = std::fopen(path.c_str(), "r");
    if (!fp) { std::fprintf(stderr, "config: cannot open %s\n", path.c_str()); return -1; }

    char line[4096];
    while (std::fgets(line, sizeof(line), fp)) {
        if (line[0] == '*' || line[0] == '#' || line[0] == '\n') continue;

        // split first whitespace-token from rest
        char key[64] = {0};
        int  ko = 0;
        const char* p = line;
        while (*p && std::isspace(static_cast<unsigned char>(*p))) ++p;
        while (*p && !std::isspace(static_cast<unsigned char>(*p)) && ko < 63) key[ko++] = *p++;
        while (*p && std::isspace(static_cast<unsigned char>(*p))) ++p;
        std::string val = strip_trailing_comment(p ? std::string(p) : std::string());
        if (!*key) continue;

        // svfits keys
        if (!std::strcmp(key, "NFILE"))         out.nfile     = std::atoi(val.c_str());
        else if (!std::strcmp(key, "HAVE_IDX")) out.have_idx  = std::atoi(val.c_str());
        else if (!std::strcmp(key, "PATH"))     out.path      = val;
        else if (!std::strcmp(key, "INPUT"))    out.input     = split_csv(val);
        else if (!std::strcmp(key, "FITS"))     out.fits_uv   = val;
        else if (!std::strcmp(key, "ALL_CHAN")) out.all_chan  = std::atoi(val.c_str());
        else if (!std::strcmp(key, "ALL_DATA")) out.all_data  = std::atoi(val.c_str());
        else if (!std::strcmp(key, "LTA"))      out.lta       = std::atof(val.c_str());
        else if (!std::strcmp(key, "N_LTA"))    out.n_lta     = std::atoi(val.c_str());
        else if (!std::strcmp(key, "NCHAV"))    out.nchav     = std::atoi(val.c_str());
        else if (!std::strcmp(key, "OBS_MJD"))  out.obs_mjd   = std::atof(val.c_str());
        else if (!std::strcmp(key, "FREQ_SET")) {
            // "freq0:freq1:nchan"
            double f0 = 0, f1 = 0; int nc = 4096;
            std::sscanf(val.c_str(), "%lf:%lf:%d", &f0, &f1, &nc);
            out.freq0_hz = f0; out.freq1_hz = f1; out.nchan = nc;
        }
        else if (!std::strcmp(key, "ANTMASK"))     out.antmask     = static_cast<uint32_t>(std::strtoul(val.c_str(), nullptr, 0));
        else if (!std::strcmp(key, "EPOCH"))       out.epoch       = std::atof(val.c_str());
        else if (!std::strcmp(key, "STOKES_TYPE")) out.stokes_type = std::atoi(val.c_str());
        else if (!std::strcmp(key, "IATUTC"))      out.iatutc      = std::atof(val.c_str());
        else if (!std::strcmp(key, "RA_MEAN"))     out.ra_mean     = std::atof(val.c_str());
        else if (!std::strcmp(key, "DEC_MEAN"))    out.dec_mean    = std::atof(val.c_str());
        else if (!std::strcmp(key, "RA_APP"))      out.ra_app      = std::atof(val.c_str());
        else if (!std::strcmp(key, "DEC_APP"))     out.dec_app     = std::atof(val.c_str());
        else if (!std::strcmp(key, "BURST_NAME"))  out.burst.name   = val;
        else if (!std::strcmp(key, "BURST_MJD"))   out.burst.mjd    = std::atof(val.c_str());
        else if (!std::strcmp(key, "BURST_TIME"))  out.burst.t      = std::atof(val.c_str());
        else if (!std::strcmp(key, "BURST_DT"))    out.burst.dt     = std::atof(val.c_str());
        else if (!std::strcmp(key, "BURST_INTWD")) out.burst.int_wd = std::atof(val.c_str());
        else if (!std::strcmp(key, "BURST_WIDTH")) out.burst.width  = std::atof(val.c_str());
        else if (!std::strcmp(key, "BURST_DM"))    out.burst.DM     = std::atof(val.c_str());
        else if (!std::strcmp(key, "BURST_DDM"))   out.burst.dDM    = std::atof(val.c_str());
        else if (!std::strcmp(key, "BURST_FREQ"))  out.burst.freq   = std::atof(val.c_str());
        else if (!std::strcmp(key, "BURST_BM_ID")) out.burst.bm_id  = std::atoi(val.c_str());
        else if (!std::strcmp(key, "BURST_RA"))    out.burst.ra_app = std::atof(val.c_str());
        else if (!std::strcmp(key, "BURST_DEC"))   out.burst.dec_app= std::atof(val.c_str());
        else if (!std::strcmp(key, "UPDATE_BURST")) out.update_burst= std::atoi(val.c_str());
        else if (!std::strcmp(key, "DO_BAND"))     out.do_band     = std::atoi(val.c_str());
        else if (!std::strcmp(key, "DO_BASE"))     out.do_base     = std::atoi(val.c_str());
        else if (!std::strcmp(key, "DO_FLAG"))     out.do_flag     = std::atoi(val.c_str());
        else if (!std::strcmp(key, "THRESH"))      out.thresh      = static_cast<float>(std::atof(val.c_str()));
        else if (!std::strcmp(key, "POST_CORR"))   out.post_corr   = std::atoi(val.c_str());
        else if (!std::strcmp(key, "DROP_CSQ"))    out.drop_csq    = std::atoi(val.c_str());
        else if (!std::strcmp(key, "RECENTRE"))    out.recentre    = std::atoi(val.c_str());
        else if (!std::strcmp(key, "NUM_THREADS")) out.num_threads = std::atoi(val.c_str());
        // pico_image extensions
        else if (!std::strcmp(key, "NPIX"))          out.npix          = std::atoi(val.c_str());
        else if (!std::strcmp(key, "CELLSIZE_ASEC")) out.cellsize_asec = std::atof(val.c_str());
        else if (!std::strcmp(key, "WEIGHTING"))     out.weighting     = val;
        else if (!std::strcmp(key, "CLEAN_ITER"))    out.clean_iter    = std::atoi(val.c_str());
        else if (!std::strcmp(key, "CLEAN_GAIN"))    out.clean_gain    = static_cast<float>(std::atof(val.c_str()));
        else if (!std::strcmp(key, "CLEAN_THRESH"))  out.clean_thresh  = static_cast<float>(std::atof(val.c_str()));
        else if (!std::strcmp(key, "IMAGE_FITS"))    out.image_fits    = val;
        else if (!std::strcmp(key, "REF_FREQ_HZ"))   out.ref_freq_hz   = std::atof(val.c_str());
        else if (!std::strcmp(key, "ANT_HDR"))       out.ant_hdr_path  = val;
        // silently ignore other keys (svfits behavior)
    }
    std::fclose(fp);

    if (out.input.size() < static_cast<std::size_t>(out.nfile)) {
        std::fprintf(stderr, "config: INPUT lists %zu files but NFILE=%d\n",
                     out.input.size(), out.nfile);
        return -1;
    }
    if (out.ant_hdr_path.empty()) out.ant_hdr_path = "antsamp.hdr";
    return 0;
}

double resolve_ref_freq(const Config& c) {
    if (c.ref_freq_hz > 0.0) return c.ref_freq_hz;
    return 0.5 * (c.freq0_hz + c.freq1_hz);
}

} // namespace pico

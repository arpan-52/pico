// pico_image/src/main.cpp
//   pico_image -c svfits_par.txt   [-a antsamp.hdr]
//
// Mirrors svfits's CLI surface so existing PICO job scripts can swap binaries
// with minimal change.

#include "pico/config.hpp"
#include "pico/pipeline.hpp"
#include <cstdio>
#include <cstring>
#include <unistd.h>

static void usage() {
    std::fprintf(stderr,
        "Usage: pico_image -c <par.txt> [-a antsamp.hdr] [-o IMAGE.FITS]\n"
        "  Imaging keys read from <par.txt> (extends svfits_par.txt):\n"
        "    NPIX, CELLSIZE_ASEC, WEIGHTING, CLEAN_ITER, CLEAN_GAIN,\n"
        "    CLEAN_THRESH, IMAGE_FITS, REF_FREQ_HZ, ANT_HDR\n");
}

int main(int argc, char** argv) {
    std::string parfile = "svfits_par.txt";
    std::string anthdr;
    std::string outfits;
    int opt;
    while ((opt = getopt(argc, argv, "c:a:o:h")) != -1) {
        switch (opt) {
            case 'c': parfile = optarg; break;
            case 'a': anthdr  = optarg; break;
            case 'o': outfits = optarg; break;
            case 'h':
            default:  usage(); return (opt == 'h') ? 0 : 2;
        }
    }
    pico::Config cfg;
    if (pico::load_config(parfile, cfg) < 0) return 1;
    if (!anthdr.empty())  cfg.ant_hdr_path = anthdr;
    if (!outfits.empty()) cfg.image_fits   = outfits;
    return (pico::run_pipeline(cfg) == 0) ? 0 : 1;
}

// pico_image/src/antsamp.cpp
// Parse the GMRT antsamp.hdr file. Two sections (Antenna.def, Sampler.def)
// each delimited by '{' / '}'.  Comments start with '*'.

#include "pico/antsamp.hpp"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cctype>
#include <string>

namespace pico {

namespace {

bool eat_section(FILE* fp, const char* tag, char* buf, std::size_t bufsz) {
    while (std::fgets(buf, static_cast<int>(bufsz), fp)) {
        if (buf[0] == '*') continue;
        if (buf[0] == '{' && std::strstr(buf, tag)) return true;
    }
    return false;
}

int parse_antennas(FILE* fp, AntSamp& out) {
    out.antenna.assign(kMaxAnts, Antenna{});
    for (int i = 0; i < kMaxAnts; ++i)
        for (int b = 0; b < kMaxBands; ++b) out.antenna[i].samp_id[b] = kMaxSamps;

    char line[1024];
    while (std::fgets(line, sizeof(line), fp)) {
        if (line[0] == '*') continue;
        if (line[0] == '}') break;
        if (std::strncmp(line, "ANT", 3) != 0) continue;
        int id = -1;
        std::sscanf(line + 3, "%d", &id);
        if (id < 0 || id >= kMaxAnts) continue;
        char name[16] = {0};
        double bx, by, bz, fixed = 0;
        if (std::sscanf(line + 10, "%15s %lf %lf %lf %lf",
                        name, &bx, &by, &bz, &fixed) < 4) continue;
        out.antenna[id].name = name;
        out.antenna[id].bx = bx;
        out.antenna[id].by = by;
        out.antenna[id].bz = bz;
        for (int b = 0; b < kMaxBands; ++b) out.antenna[id].d0[b] = fixed;
    }
    return 0;
}

int parse_samplers(FILE* fp, AntSamp& out) {
    // SPOTLIGHT: USB-130 → band 0 (RR), USB-175 → band 1 (LL)
    auto band_id = [](const char* s) -> int {
        if (!std::strncmp(s, "USB-130", 7)) return 0;
        if (!std::strncmp(s, "USB-175", 7)) return 1;
        if (!std::strncmp(s, "LSB-130", 7)) return 2;
        if (!std::strncmp(s, "LSB-175", 7)) return 3;
        return -1;
    };

    out.sampler.clear();
    out.antmask = 0;
    out.bandmask = 0;

    char line[1024];
    while (std::fgets(line, sizeof(line), fp)) {
        if (line[0] == '*') continue;
        if (line[0] == '}') break;
        if (std::strncmp(line, "SMP", 3) != 0) continue;

        int id = -1;
        std::sscanf(line + 3, "%d", &id);
        if (id < 0 || id >= kMaxSamps) continue;

        char ant_name[16] = {0}, band_name[32] = {0};
        if (std::sscanf(line + 10, "%15s %31s", ant_name, band_name) < 2) continue;

        int ant_id = -1;
        for (int a = 0; a < kMaxAnts; ++a) {
            if (out.antenna[a].name == ant_name) { ant_id = a; break; }
        }
        if (ant_id < 0) continue;
        int b = band_id(band_name);
        if (b < 0) continue;

        Sampler s;
        s.ant_id = ant_id;
        s.band   = b;
        out.antenna[ant_id].samp_id[b] = static_cast<int>(out.sampler.size());
        out.sampler.push_back(s);
        out.antmask  |= (1u << ant_id);
        out.bandmask |= (1u << b);
    }
    out.nsamp = static_cast<int>(out.sampler.size());
    return 0;
}

} // namespace

int load_antsamp(const std::string& path, AntSamp& out) {
    FILE* fp = std::fopen(path.c_str(), "r");
    if (!fp) {
        std::fprintf(stderr, "antsamp: cannot open %s\n", path.c_str());
        return -1;
    }

    char line[1024];
    if (!eat_section(fp, "Antenna", line, sizeof(line))) {
        std::fprintf(stderr, "antsamp: no Antenna section\n"); std::fclose(fp); return -1;
    }
    parse_antennas(fp, out);

    if (!eat_section(fp, "Sampler", line, sizeof(line))) {
        std::fprintf(stderr, "antsamp: no Sampler section\n"); std::fclose(fp); return -1;
    }
    parse_samplers(fp, out);

    std::fclose(fp);

    // Build baselines using default (full) antmask. Caller can rebuild later.
    rebuild_baselines(out, out.antmask);
    return 0;
}

int rebuild_baselines(AntSamp& as, uint32_t user_antmask) {
    as.baseline.clear();
    const uint32_t am = as.antmask & user_antmask;
    // svfits get_baseline order: per band, s0 then s1 ≥ s0, both samplers
    // must have band == band and belong to enabled antennas.
    for (int band = 0; band < as.stokes; ++band) {
        for (int s0 = 0; s0 < as.nsamp; ++s0) {
            const Sampler& a = as.sampler[s0];
            if (a.band != band) continue;
            if (!((1u << a.ant_id) & am)) continue;
            for (int s1 = s0; s1 < as.nsamp; ++s1) {
                const Sampler& b = as.sampler[s1];
                if (b.band != band) continue;
                if (!((1u << b.ant_id) & am)) continue;
                BaselineEntry be; be.s0 = a; be.s1 = b;
                as.baseline.push_back(be);
            }
        }
    }
    as.nbase = static_cast<int>(as.baseline.size());
    return 0;
}

} // namespace pico

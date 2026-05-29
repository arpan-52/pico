// pico_image/src/raw_io.cpp
// Persistent-FD reader for the SPOTLIGHT raw visibility files. Each file's
// FD stays open for the life of the run, and reads use pread() so threads
// can pull different slices in parallel without locking.
//
// Per-slice layout (HAVE_IDX=1):
//   16 bytes  struct timeval  { time_t tv_sec; suseconds_t tv_usec; }
//    4 bytes  int file_index  (1-based on disk)
//   rec_per_slice * (channels * baselines * 4) bytes  visibility data
//
// HAVE_IDX=0 drops the 4-byte index field.

#include "pico/raw_io.hpp"
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include <cstdio>
#include <cstring>
#include <ctime>

namespace pico {

namespace {

std::size_t per_slice_hdr_size(int have_idx) {
    // sizeof(struct timeval) is 16 on x86_64 Linux, 8 on 32-bit. We hard-code
    // 16 because the GMRT cluster is x86_64 (svfits assumes this too).
    return have_idx ? (16 + sizeof(int)) : 16;
}

} // namespace

int open_raw_set(const Config& cfg, const AntSamp& as, RawSet& out) {
    out.files.clear();
    out.files.resize(cfg.nfile);
    out.channels  = cfg.nchan;
    out.baselines = as.nbase;
    out.recl      = static_cast<std::size_t>(out.channels) *
                    static_cast<std::size_t>(out.baselines) * sizeof(float);

    // svfits defaults: rec_per_slice=50, t_slice = NCHAN*sta/clock * 50 ≈ 65 ms ?
    // svfits inits rfile->t_slice from corr params; we read the first file's
    // record count to infer per-slice timing.
    out.rec_per_slice  = 50;   // MAX_REC_PER_SLICE in svfits ; SPOTLIGHT default
    out.t_slice        = 50 * 0.001310720; // 1.31072 ms per integration × 50
    out.slice_interval = cfg.nfile * out.t_slice;
    out.timeval_size   = 16;

    for (int i = 0; i < cfg.nfile; ++i) {
        RawFile& f = out.files[i];
        f.file_idx = i;
        f.path = cfg.path;
        if (!f.path.empty() && f.path.back() != '/') f.path.push_back('/');
        f.path += cfg.input[i];
        // strip trailing newline / whitespace
        while (!f.path.empty() && (f.path.back() == '\n' || f.path.back() == ' '
                                   || f.path.back() == '\t'))
            f.path.pop_back();
        f.fd = ::open(f.path.c_str(), O_RDONLY);
        if (f.fd < 0) {
            std::fprintf(stderr, "raw_io: cannot open %s: %s\n",
                         f.path.c_str(), std::strerror(errno));
            return -1;
        }
    }

    // Read slice-0 timestamp from file 0 to establish mjd_ref.
    if (out.files.empty()) return -1;
    struct timeval tv{};
    const std::size_t hdrsz = per_slice_hdr_size(cfg.have_idx);
    if (::pread(out.files[0].fd, &tv, sizeof(tv), 0) != static_cast<ssize_t>(sizeof(tv))) {
        std::fprintf(stderr, "raw_io: cannot read timestamp from %s\n",
                     out.files[0].path.c_str());
        return -1;
    }
    // Unix epoch → MJD: MJD = unix_sec / 86400 + 40587
    const double unix_sec = static_cast<double>(tv.tv_sec) + tv.tv_usec * 1e-6;
    out.mjd_ref = unix_sec / 86400.0 + 40587.0;

    (void)hdrsz;
    return 0;
}

void close_raw_set(RawSet& rs) {
    for (auto& f : rs.files) {
        if (f.fd >= 0) { ::close(f.fd); f.fd = -1; }
    }
}

int read_slice(const RawSet& rs, int idx, int slice, void* dst) {
    if (idx < 0 || idx >= static_cast<int>(rs.files.size())) return -1;
    const RawFile& f = rs.files[idx];

    // svfits get_slice_time computes offset as:
    //   off = slice * (hdrsize + rec_per_slice * recl) + timesize
    // Note svfits uses `timesize` (16) as the post-header offset even when
    // HAVE_IDX=1; the int index is read from the slot AFTER timeval but data
    // still starts at timesize+sizeof(int). Stay consistent: skip whole header.
    // We replicate svfits's layout: data start = full per-slice header.
    const std::size_t per_slice = 16 /*timeval*/ + sizeof(int) /*idx, if HAVE_IDX*/
                                 + rs.rec_per_slice * rs.recl;
    const off_t base = static_cast<off_t>(slice) * static_cast<off_t>(per_slice);
    const off_t data_off = base + 16 + static_cast<off_t>(sizeof(int));

    const std::size_t want = static_cast<std::size_t>(rs.rec_per_slice) * rs.recl;
    std::size_t got = 0;
    auto* p = static_cast<char*>(dst);
    while (got < want) {
        ssize_t n = ::pread(f.fd, p + got, want - got, data_off + static_cast<off_t>(got));
        if (n <= 0) {
            if (n == 0) {
                std::fprintf(stderr,
                    "raw_io: EOF on %s slice %d at byte %zu/%zu (offset %lld)\n",
                    f.path.c_str(), slice, got, want,
                    (long long)(data_off + (off_t)got));
                return -1;
            }
            if (errno == EINTR) continue;
            std::fprintf(stderr, "raw_io: pread failed on %s slice %d: %s\n",
                         f.path.c_str(), slice, std::strerror(errno));
            return -1;
        }
        got += n;
    }
    return 0;
}

int n_slices(const RawSet& rs, int idx) {
    if (idx < 0 || idx >= static_cast<int>(rs.files.size())) return 0;
    struct stat st{};
    if (::fstat(rs.files[idx].fd, &st) < 0) return 0;
    const std::size_t per_slice = 16 + sizeof(int) + rs.rec_per_slice * rs.recl;
    return static_cast<int>(st.st_size / per_slice);
}

} // namespace pico

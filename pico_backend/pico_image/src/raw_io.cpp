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
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>

namespace pico {

int open_raw_set(const Config& cfg, const AntSamp& as, RawSet& out) {
    out.files.clear();
    out.files.resize(cfg.nfile);
    out.channels  = cfg.nchan;
    out.baselines = as.nbase;
    out.recl      = static_cast<std::size_t>(out.channels) *
                    static_cast<std::size_t>(out.baselines) * sizeof(float);

    // Slice timing derived as in svfits init_corr (svsubs.c:460-477,670):
    //   statime = 2*channels*sta/clock ; integ = lta*statime
    //   t_slice = rec_per_slice * integ
    // SPOTLIGHT correlator constants: clock=400 MHz, sta=1, lta=64.
    const double clock_hz = 400e6;
    const double sta      = 1.0;
    const double lta      = 64.0;
    const double statime  = 2.0 * cfg.nchan * sta / clock_hz;
    out.rec_per_slice  = 50;   // 50 lta records per slice (svsubs.c:668)
    out.t_slice        = out.rec_per_slice * lta * statime;
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

    // Per-file time anchoring, ported from svfits get_slice_time
    // (svsubs.c:813-887): read each file's slice-0 timestamp AND the embedded
    // file index (1-based int after the timeval), and IGNORE the order the
    // files were supplied in. The timestamp marks the start of file[0]'s
    // slice (svsubs.c:860), so
    //   t_start[file] = (unix_file - unix_ref) + embedded_idx * t_slice
    // with embedded_idx 0-based after the -1. svfits explicitly distrusts the
    // supplied order ("read the index ... from the file itself", svsubs.c:846).
    if (out.files.empty()) return -1;
    double unix_ref = 0.0;
    for (int i = 0; i < cfg.nfile; ++i) {
        RawFile& f = out.files[i];
        struct timeval tv{};
        if (::pread(f.fd, &tv, sizeof(tv), 0) != static_cast<ssize_t>(sizeof(tv))) {
            std::fprintf(stderr, "raw_io: cannot read timestamp from %s\n",
                         f.path.c_str());
            return -1;
        }
        const double unix_sec = static_cast<double>(tv.tv_sec) + tv.tv_usec * 1e-6;

        int k = i;  // HAVE_IDX=0: fall back to list position
        if (cfg.have_idx) {
            int idx_raw = 0;
            if (::pread(f.fd, &idx_raw, sizeof(idx_raw), 16) !=
                static_cast<ssize_t>(sizeof(idx_raw))) {
                std::fprintf(stderr, "raw_io: cannot read file index from %s\n",
                             f.path.c_str());
                return -1;
            }
            k = idx_raw - 1;  // 1-based on disk (svsubs.c:852)
            if (k < 0 || k >= cfg.nfile) {
                std::fprintf(stderr, "raw_io: idx %d out of range in %s\n",
                             idx_raw, f.path.c_str());
                return -1;
            }
        }
        f.file_idx = k;

        if (i == 0) {
            unix_ref = unix_sec;
            // Unix epoch → MJD: MJD = unix_sec / 86400 + 40587
            out.mjd_ref = unix_sec / 86400.0 + 40587.0;
        }
        f.t_start = (unix_sec - unix_ref) + k * out.t_slice;
        if (k != i)
            std::fprintf(stderr,
                "raw_io: WARNING list_pos=%d but embedded idx=%d for %s "
                "(supplied order wrong; using embedded idx)\n",
                i, k, f.path.c_str());
        std::fprintf(stderr,
            "raw_io: file[%2d] idx=%2d tv=%ld.%06ld t_start=%.6f s  %s\n",
            i, k, static_cast<long>(tv.tv_sec), static_cast<long>(tv.tv_usec),
            f.t_start, f.path.c_str());
    }

    // Cross-check the clock-derived t_slice against the data — but ONLY when
    // file[0] genuinely stacks a second timestamped slice. svfits supports
    // multiple slices per file (svsubs.c:874-881), where consecutive slices on
    // the same file are slice_interval = nfile*t_slice apart and each carries
    // its own timestamp. SPOTLIGHT burst dumps, however, are one slice per file:
    // all nfile files carry file[0]'s identical timestamp and are ordered solely
    // by the embedded index (t_start = idx*t_slice — see the loop above). There
    // is then no second timestamp inside file[0]; a read at +per_slice lands in
    // visibility (half-float) data, which decoded as time_t is absurd. Treat an
    // implausible second timestamp as "single slice per file" and skip — exactly
    // like the file-too-short case — instead of falsely rejecting good data.
    {
        const std::size_t per_slice = 16 + (cfg.have_idx ? sizeof(int) : 0)
                                    + out.rec_per_slice * out.recl;
        struct timeval tv0{}, tv1{};
        if (::pread(out.files[0].fd, &tv0, sizeof(tv0), 0) == (ssize_t)sizeof(tv0) &&
            ::pread(out.files[0].fd, &tv1, sizeof(tv1),
                    static_cast<off_t>(per_slice)) == (ssize_t)sizeof(tv1)) {
            const double dt = (tv1.tv_sec - tv0.tv_sec)
                            + (tv1.tv_usec - tv0.tv_usec) * 1e-6;
            // A real second slice is ~slice_interval (≈1 s) after the first.
            // Anything outside (0, 1 day] or with an out-of-range usec is data,
            // not a timestamp ⇒ one slice per file, nothing to cross-check.
            const bool tv1_is_time = tv1.tv_usec >= 0 && tv1.tv_usec < 1000000
                                   && dt > 0.0 && dt < 86400.0;
            if (!tv1_is_time) {
                std::fprintf(stderr,
                    "raw_io: single slice per file (no 2nd timestamp in file[0]); "
                    "t_slice derived=%.9f s [clock-derived, unchecked]\n",
                    out.t_slice);
            } else {
                const double t_slice_meas = dt / cfg.nfile;
                const double tiny = (out.t_slice / out.rec_per_slice) * 1.0e-3;
                std::fprintf(stderr,
                    "raw_io: t_slice derived=%.9f s, measured from timestamps=%.9f s\n",
                    out.t_slice, t_slice_meas);
                if (std::fabs(t_slice_meas - out.t_slice) > tiny) {
                    std::fprintf(stderr,
                        "raw_io: ERROR t_slice mismatch — correlator constants "
                        "(clock/sta/lta) wrong for these data\n");
                    return -1;
                }
            }
        }
        // (file too short for 2 slices: skip the check)
    }
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

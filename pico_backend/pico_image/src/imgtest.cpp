// pico_image/src/imgtest.cpp
// Imaging-path debugger. Drives the REAL pico imaging functions
//   push_sample -> global_clip_samples -> to_finufft_coords -> run_finufft_dirty
//   -> hogbom_clean
// with synthetic visibilities, to isolate gridder / weighting / flagger behaviour
// from the raw-data calibration. Each scenario also dumps its dirty image to
//   imgtest_out/scenNN.f64   (int32 npix, then npix*npix float64, row-major)
// for plotting (see plot_imgtest.py).

#include "pico/grid.hpp"
#include "pico/clean.hpp"
#include "pico/types.hpp"

#include <cstdio>
#include <cstdint>
#include <cmath>
#include <random>
#include <vector>
#include <algorithm>
#include <complex>
#include <string>
#include <sys/stat.h>

using namespace pico;
using cd = std::complex<double>;

namespace {

struct UV { double u, v; };

struct Array {
    std::vector<UV> uv;       // size Nbase (same uv every record = snapshot)
    int Nbase, Nrec;
    double Umax;
};

Array make_array(int Nbase, int Nrec, double Umax, unsigned seed = 1) {
    Array a; a.Nbase = Nbase; a.Nrec = Nrec; a.Umax = Umax;
    std::mt19937 rng(seed);
    std::normal_distribution<double> g(0.0, Umax / 3.0);
    a.uv.resize(Nbase);
    for (int b = 0; b < Nbase; ++b) {
        double u = std::max(-Umax, std::min(Umax, g(rng)));
        double v = std::max(-Umax, std::min(Umax, g(rng)));
        a.uv[b] = {u, v};
    }
    return a;
}

struct ImgStat { double centre, vmax, vmin, noise; int px, py; };

ImgStat stat_image(const std::vector<double>& img, int n) {
    ImgStat s{};
    const int cx = n / 2, cy = n / 2;
    s.centre = img[(std::size_t)cy * n + cx];
    s.vmax = -1e300; s.vmin = 1e300; s.px = cx; s.py = cy;
    for (int y = 0; y < n; ++y)
        for (int x = 0; x < n; ++x) {
            double v = img[(std::size_t)y * n + x];
            if (v > s.vmax) { s.vmax = v; s.px = x; s.py = y; }
            if (v < s.vmin) s.vmin = v;
        }
    std::vector<double> tmp(img);
    auto med = [](std::vector<double>& a) {
        std::nth_element(a.begin(), a.begin() + a.size() / 2, a.end());
        return a[a.size() / 2];
    };
    double m = med(tmp);
    for (double& v : tmp) v = std::fabs(v - m);
    s.noise = 1.4826 * med(tmp);
    return s;
}

void dump_image(const std::string& path, const std::vector<double>& img, int n) {
    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) { std::fprintf(stderr, "cannot write %s\n", path.c_str()); return; }
    std::int32_t npix = n;
    std::fwrite(&npix, sizeof(npix), 1, f);
    std::fwrite(img.data(), sizeof(double), img.size(), f);
    std::fclose(f);
}

// Per-baseline DC ("zero-spacing") subtraction: for each baseline subtract the
// mean visibility over the records flagged off-source. off_mask empty -> use all
// records (pure stationary-RFI removal). Mirrors svfits off_src in spirit.
void dc_subtract(std::vector<cd>& vis, const Array& A,
                 const std::vector<char>& off_mask) {
    for (int b = 0; b < A.Nbase; ++b) {
        cd sum(0, 0); int n = 0;
        for (int r = 0; r < A.Nrec; ++r) {
            if (!off_mask.empty() && !off_mask[r]) continue;
            sum += vis[(std::size_t)r * A.Nbase + b]; ++n;
        }
        cd dc = (n > 0) ? sum / double(n) : cd(0, 0);
        for (int r = 0; r < A.Nrec; ++r)
            vis[(std::size_t)r * A.Nbase + b] -= dc;
    }
}

// Run one scenario. image_mask: which records to grid (empty -> all).
struct Result { ImgStat ds; std::size_t flagged, kept; double cflux, cmaxf; int cmx, cmy; };

Result run_scenario(const std::string& name, const std::string& outfile,
                    const Array& A, const std::vector<cd>& vis,
                    const std::vector<char>& image_mask,
                    bool do_clip, double thresh, int npix, double cellsize_rad) {
    GridSamples dirty, psf;
    for (int r = 0; r < A.Nrec; ++r) {
        if (!image_mask.empty() && !image_mask[r]) continue;
        for (int b = 0; b < A.Nbase; ++b) {
            const UV p = A.uv[b];
            cd V = vis[(std::size_t)r * A.Nbase + b];
            Vis vd; vd.r = (float)V.real(); vd.i = (float)V.imag(); vd.wt = 1.0f;
            UvwSec g{p.u, p.v, 0.0};
            push_sample(dirty, g, 1.0, vd);
            Vis vp; vp.r = 1.0f; vp.i = 0.0f; vp.wt = 1.0f;
            push_sample(psf, g, 1.0, vp);
        }
    }
    const std::size_t n_in = dirty.size();
    std::size_t flagged = 0;
    if (do_clip) flagged = global_clip_samples(dirty, psf, thresh, 0);

    to_finufft_coords(dirty.u_lambda, dirty.v_lambda, cellsize_rad);
    to_finufft_coords(psf.u_lambda,   psf.v_lambda,   cellsize_rad);

    DirtyOut d, p;
    run_finufft_dirty(dirty, npix, npix, 1e-7, d);
    run_finufft_dirty(psf,   npix, npix, 1e-7, p);
    double pk = 0; for (double v : p.img) pk = std::max(pk, std::fabs(v));
    if (pk > 0) for (double& v : p.img) v /= pk;

    Result R{};
    R.ds = stat_image(d.img, npix);
    R.flagged = flagged; R.kept = dirty.size();

    CleanParams cp; cp.niter = 2000; cp.gain = 0.1f; cp.threshold_sigma = 5.0f;
    CleanOut co;
    hogbom_clean(d.img, p.img, npix, npix, cp, co);
    R.cflux = 0; R.cmaxf = 0; R.cmx = npix/2; R.cmy = npix/2;
    for (auto& c : co.components) {
        R.cflux += c.flux;
        if (std::fabs(c.flux) > std::fabs(R.cmaxf)) { R.cmaxf = c.flux; R.cmx = c.ix; R.cmy = c.iy; }
    }
    dump_image(outfile, d.img, npix);

    std::printf("\n=== %s ===\n", name.c_str());
    std::printf("  in=%zu flagged=%zu(%.2f%%) kept=%zu | DIRTY centre=%.4g max=%.4g@off(%+d,%+d) "
                "min=%.4g noise=%.4g pk/no=%.1f | CLEAN n=%zu flux=%.4g bright=%.4g@off(%+d,%+d)\n",
                n_in, flagged, 100.0*flagged/std::max<std::size_t>(1,n_in), R.kept,
                R.ds.centre, R.ds.vmax, R.ds.px-npix/2, R.ds.py-npix/2,
                R.ds.vmin, R.ds.noise, R.ds.noise>0?R.ds.vmax/R.ds.noise:0.0,
                co.components.size(), R.cflux, R.cmaxf, R.cmx-npix/2, R.cmy-npix/2);
    return R;
}

} // namespace

int main() {
    const int    npix = 256;
    const int    Nbase = 400, Nrec = 100;
    const double Umax  = 4000.0;
    const double cellsize_rad = 0.12 / Umax;
    const double rfi_A = 130.0;
    const double noise_sigma = 1.0;

    ::mkdir("imgtest_out", 0755);

    Array A = make_array(Nbase, Nrec, Umax);
    const std::size_t Nsamp = (std::size_t)Nbase * Nrec;

    std::mt19937 rng(12345);
    std::normal_distribution<double> gn(0.0, noise_sigma);
    auto noise = [&]() { return cd(gn(rng), gn(rng)); };

    const double l0 = 30 * cellsize_rad, m0 = 20 * cellsize_rad;  // off-centre src
    auto src_phase = [&](int b) {
        return 2 * M_PI * (A.uv[b].u * l0 + A.uv[b].v * m0);
    };
    std::vector<char> none;  // empty mask = all records

    // 1. point source @ centre
    {
        std::vector<cd> v(Nsamp, cd(1, 0));
        run_scenario("1. point src @ centre (S=1)", "imgtest_out/scen01.f64",
                     A, v, none, false, 20, npix, cellsize_rad);
    }
    // 2. point source off-centre
    {
        std::vector<cd> v(Nsamp);
        for (int r=0;r<Nrec;++r) for (int b=0;b<Nbase;++b)
            v[(std::size_t)r*Nbase+b] = std::polar(1.0, src_phase(b));
        run_scenario("2. point src off-centre (S=1)", "imgtest_out/scen02.f64",
                     A, v, none, false, 20, npix, cellsize_rad);
    }
    // 3. pure noise
    {
        std::vector<cd> v(Nsamp); for (auto& x:v) x=noise();
        run_scenario("3. pure noise (blank)", "imgtest_out/scen03.f64",
                     A, v, none, false, 20, npix, cellsize_rad);
    }
    // 4. noise + SPIKY rfi + clip
    {
        std::vector<cd> v(Nsamp); for (auto& x:v) x=noise();
        std::uniform_real_distribution<double> uni(0,1), ph(0,2*M_PI);
        for (auto& x:v) if (uni(rng)<0.02) x += std::polar(2000.0, ph(rng));
        run_scenario("4. noise + SPIKY rfi + CLIP", "imgtest_out/scen04.f64",
                     A, v, none, true, 20, npix, cellsize_rad);
    }
    // 5. noise + COHERENT-DC rfi + clip  (FAILS)
    {
        std::vector<cd> v(Nsamp); for (auto& x:v) x=noise()+cd(rfi_A,0);
        run_scenario("5. noise + COHERENT-DC rfi + CLIP (fails)", "imgtest_out/scen05.f64",
                     A, v, none, true, 20, npix, cellsize_rad);
    }
    // 6. noise + COHERENT-DC rfi + DC-SUBTRACTION  (FIX)
    {
        std::vector<cd> v(Nsamp); for (auto& x:v) x=noise()+cd(rfi_A,0);
        dc_subtract(v, A, none);   // per-baseline mean over all records
        run_scenario("6. noise + COHERENT-DC rfi + DC-SUB (fix)", "imgtest_out/scen06.f64",
                     A, v, none, false, 20, npix, cellsize_rad);
    }
    // transient burst: source present only in records [45,55)
    std::vector<char> burst(Nrec,0), offsrc(Nrec,1);
    for (int r=45;r<55;++r){ burst[r]=1; offsrc[r]=0; }
    const double S_burst = 8.0;
    // 7. burst + COHERENT-DC rfi, CLIP only, image burst recs  (source buried)
    {
        std::vector<cd> v(Nsamp);
        for (int r=0;r<Nrec;++r) for (int b=0;b<Nbase;++b) {
            cd x = noise() + cd(rfi_A,0);
            if (burst[r]) x += std::polar(S_burst, src_phase(b));
            v[(std::size_t)r*Nbase+b]=x;
        }
        run_scenario("7. burst(S=8) + DC rfi, CLIP only (buried)", "imgtest_out/scen07.f64",
                     A, v, burst, true, 20, npix, cellsize_rad);
    }
    // 8. burst + COHERENT-DC rfi, OFF-SOURCE DC-SUB, image burst recs  (recovered)
    {
        std::vector<cd> v(Nsamp);
        for (int r=0;r<Nrec;++r) for (int b=0;b<Nbase;++b) {
            cd x = noise() + cd(rfi_A,0);
            if (burst[r]) x += std::polar(S_burst, src_phase(b));
            v[(std::size_t)r*Nbase+b]=x;
        }
        dc_subtract(v, A, offsrc);   // subtract per-baseline mean over OFF-burst recs
        run_scenario("8. burst(S=8) + DC rfi, OFF-SRC DC-SUB (recovered)", "imgtest_out/scen08.f64",
                     A, v, burst, true, 20, npix, cellsize_rad);
    }

    std::printf("\nwrote 8 dirty images to imgtest_out/scen0[1-8].f64\n");
    return 0;
}

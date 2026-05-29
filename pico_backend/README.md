# pico_backend

C++ imaging backend for the SPOTLIGHT/GMRT transient pipeline. Replaces the
SVFITS → UVFITS → CASA/wsclean detour with a single binary that reads the 16
raw 1.31 ms visibility dumps, applies DM-aware time/freq selection, computes
uvw + J2000 rotation, runs an MFS FINUFFT gridder, Hogbom CLEAN, and writes a
FITS image.

## Layout

```
pico_backend/
  pico_image/    # the imaging binary (sources + CMake)
  finufft/       # FINUFFT 2.6 (NUFFT library; CPU)
  SuperNOVAS/    # SuperNOVAS 1.6 (precession/nutation for J2000)
```

`pico_image/CMakeLists.txt` consumes `finufft/` and `SuperNOVAS/` as
sibling `add_subdirectory()` builds. cfitsio is expected from the system
(install `libcfitsio-dev` on Debian/Ubuntu, `cfitsio-devel` on Fedora).

## Build

```bash
cd pico_image
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Produces `pico_image/build/pico_image` (statically links FINUFFT and
SuperNOVAS; dynamically links cfitsio + OpenMP).

## Run

```bash
pico_image -c svfits_par.txt -a antsamp.hdr
```

Imaging-specific keys (extend the standard `svfits_par.txt`):

| key             | default       | meaning                                 |
|-----------------|---------------|------------------------------------------|
| `NPIX`          | 2048          | image side in pixels                     |
| `CELLSIZE_ASEC` | 1.0           | pixel size in arcsec                     |
| `WEIGHTING`     | natural       | natural \| uniform (uniform NYI)         |
| `CLEAN_ITER`    | 1000          | Hogbom minor cycles                      |
| `CLEAN_GAIN`    | 0.1           | loop gain                                |
| `CLEAN_THRESH`  | 0.0           | stop when peak < threshold (Jy)          |
| `IMAGE_FITS`    | PICO_IMAGE.FITS | output filename                        |
| `REF_FREQ_HZ`   | midband       | MFS reference freq                       |
| `ANT_HDR`       | antsamp.hdr   | path to antenna/sampler header           |

Standard svfits keys (NFILE, PATH, INPUT, ANTMASK, FREQ_SET, BURST_*,
DO_BAND, DO_BASE, DO_FLAG, THRESH, RECENTRE, NUM_THREADS, …) are honoured
unchanged. Unknown keys are silently ignored — same behavior as svfits.

## Output

A single FITS file with four image HDUs: `DIRTY`, `PSF`, `RESIDUAL`,
`RESTORED`. WCS is J2000 RA---SIN / DEC--SIN at the configured phase
centre, frequency axis at the MFS reference frequency.

## Antenna selection (ANTMASK)

Identical semantics to svfits: bit `i` in `ANTMASK` selects ANT`i` from the
order in `antsamp.hdr`. GMRT's spare/duplicate slots (ANT30 = C07, ANT31 =
S05 in the canonical antsamp.hdr) are excluded by the default mask
`ANTMASK 1073741823` (= `0x3FFFFFFF`, bits 0..29).

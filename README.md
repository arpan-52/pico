# PICO v1 — Prompt Imaging of Cosmic Outbursts

PICO is the imaging stage of the **SPOTLIGHT/GMRT** transient pipeline. When the
real-time search flags a candidate burst, PICO turns its raw visibility dump
into a science-ready radio image — automatically, on whichever compute node is
least busy, with no human in the loop.

This directory is the **distributed orchestration layer**: a trigger, a daemon,
and a node-side runner that fan out [`pico_backend/`](pico_backend/) (the C++
imaging engine) across the cluster. One command in, a FITS image and a
publication-quality PNG out.

---

## Why it exists

The legacy path imaged a burst via `SVFITS → UVFITS → CASA/wsclean` — several
processes, several intermediate files, minutes of wall-clock per burst. PICO
collapses that into a **single binary** (`pico_image`) that reads the 16 raw
1.31 ms visibility dumps, applies DM-aware time/frequency selection, computes
the `uvw` coordinates with full J2000 rotation, grids with FINUFFT, runs Högbom
CLEAN, and writes FITS — and then wraps it in lightweight scheduling so a node
is chosen, the job is tracked, and the result is plotted without operator input.

---

## Architecture

```
            pico.sh                pico_monitor.py                pico_run_job.sh
        (trigger / submit)      (daemon / scheduler)            (node-side worker)
                │                        │                              │
   svfits_par.txt + outdir              │                              │
                │                        │                              │
                ├─ parse burst params    │                              │
                ├─ write jobs_pending/<id>.job                          │
                └─ ensure daemon up ────▶ scan .job files               │
                                         ├─ record in pico_jobs.db       │
                                         ├─ pick freest of 8 nodes       │
                                         └─ ssh node ─────────────────▶  run pico_image
                                                                         ├─ STAGE 1: image  → burst_<mjd>-*.fits
                                                                         └─ STAGE 2: plot   → burst_<mjd>.png
```

### Components

| File | Role |
|------|------|
| `pico.sh` | **Trigger.** Parses the burst config, drops a `jobs_pending/<job_id>.job` descriptor, and makes sure the monitor daemon is running. Sets no imaging parameters itself. |
| `pico_monitor.py` | **Scheduler daemon.** Watches `jobs_pending/`, records state in `pico_jobs.db`, and dispatches each job over SSH to the least-loaded node (`rggpur00`–`rggpur07`). |
| `pico_run_job.sh` | **Node-side worker.** Stages the config + `antsamp.hdr`, derives imaging parameters from the observing band, runs `pico_image`, then renders the PNG. |
| `pico_status.py` | **Operator CLI.** Live view of node pressure, the job queue, and per-job logs. |
| `pico_backend/` | The C++ imaging engine (`pico_image`) plus its vendored deps (FINUFFT, SuperNOVAS). See [`pico_backend/README.md`](pico_backend/README.md). |
| `jobs_pending/` | Spool directory — one `*.job` file per submitted burst. |

### Why a custom scheduler?

Earlier versions always sent every job to `rggpur00`: node selection was
decoupled from accounting, so the "active jobs" counter never moved and the
scheduler fell back to alphabetical order. PICO v1 fixes this with a single
locked SQLite transaction — `claim_freest_server()` selects the minimum-load
node **and** increments its counter atomically; `release_server()` decrements on
completion. Load is therefore real, and bursts spread across all eight nodes.

---

## Quick start

```bash
# Submit a burst for imaging.
#   -c   the svfits-style parameter file (raw-data PATH lives inside it)
#   arg  the output directory for products + per-job log
./pico.sh -c svfits_par.txt /lustre_scratch/spotlight/data/<obs>/imaging_out
```

That is the whole interface. The trigger ensures the daemon is up, the daemon
picks a node, and within ~minutes the output directory contains the FITS
products, the PNG, and a full `pico_<job_id>.log` trace.

```bash
# Watch a specific job.
tail -f /lustre_scratch/.../imaging_out/pico_<job_id>.log

# Monitor the whole system.
python3 pico_status.py pressure --watch     # per-node load + queue summary
python3 pico_status.py jobs --status running
python3 pico_status.py detail <job_id>
python3 pico_status.py logs   <job_id>
```

---

## Configuration

PICO consumes a standard **`svfits_par.txt`** — the same format the rest of the
SPOTLIGHT chain uses. The raw-data location is **not** a command-line argument;
it is the `PATH` key inside the par file (the burst's `DT_…` directory), and
`antsamp.hdr` is expected alongside those raw files.

> **svfits par quirk:** `!` is an inline comment marker glued to the value,
> e.g. `PATH /data/DT_xxx/! directory containing raw dumps`. Parsing stops at
> the `!`.

`pico_run_job.sh` appends a band-derived imaging block to a **copy** of the par
in the output directory (the original you point `-c` at is left untouched).
Current defaults, uniform across GMRT bands:

| key | value | meaning |
|-----|-------|---------|
| `NPIX` | 4096 | image side, pixels |
| `CELLSIZE_ASEC` | 1.0 | pixel size, arcsec |
| `WEIGHTING` | natural | natural \| uniform (uniform NYI) |
| `CLEAN_ITER` | 1000 | Högbom minor cycles |
| `CLEAN_GAIN` | 0.1 | loop gain |
| `CLEAN_THRESH` | 10.0 | stop at N·σ (MAD noise) |
| `REF_FREQ_HZ` | 0 | 0 → midband, else explicit MFS reference |
| `ANT_HDR` | resolved | absolute path to `antsamp.hdr` |

See [`pico_backend/README.md`](pico_backend/README.md) for the full key
reference honoured by the imaging binary.

---

## Outputs

Written to the submission's output directory:

| file | description |
|------|-------------|
| `burst_<mjd>-image.fits` | restored (CLEANed) image — the science product |
| `burst_<mjd>-dirty.fits` | dirty image |
| `burst_<mjd>-psf.fits` | point-spread function |
| `burst_<mjd>-residual.fits` | post-CLEAN residual |
| `burst_<mjd>-model.fits` | CLEAN component model |
| `burst_<mjd>.png` | publication-quality 5′ cutout around the peak |
| `pico_<job_id>.log` | complete per-burst trace |

**Output naming is driven by the `IMAGE_FITS` par key:** `pico_image` strips the
extension and writes `<stem>-{dirty,psf,image,residual,model}.fits`. The runner
forces `IMAGE_FITS = burst_<mjd>.FITS` so products line up with what the plotter
and the success check expect.

The PNG is rendered for publication: a 5-arcmin field centred on the positive
peak, J2000 WCS axes, the synthesised clean-beam ellipse, a 1′ scale bar, an
`afmhot` colour map with a ZScale stretch, and the peak flux / signal-to-noise
in the title.

---

## Deployment

PICO runs from `/lustre_archive/apps/astrosoft/test_imaging/picov1` on the
cluster. To deploy:

1. **Build the backend** (once, or after engine changes):
   ```bash
   cd pico_backend && ./build.sh        # → pico_backend/pico_image/build/pico_image
   ```
2. **Sync this directory** to the deploy path. The scripts resolve their own
   location, so they work wherever they are dropped — but the deployed copy and
   this source tree must not drift (a stale `pico_run_job.sh` is the classic
   cause of "produced no `burst_*-image.fits`" while `PICO_IMAGE-*.fits` files
   sit in the output directory).

All three entry points activate the `pico` micromamba environment before any
Python runs, so `pico_image`, astropy, and matplotlib resolve to one consistent
interpreter over non-login SSH.

---

## Glossary

- **MFS** — multi-frequency synthesis; imaging across the full band at once.
- **DM** — dispersion measure; sets the time/frequency selection for the burst.
- **MAD** — median absolute deviation; the robust noise estimate behind the
  CLEAN threshold and the PNG signal-to-noise.
- **antsamp.hdr** — antenna/sampler header defining array geometry and the
  `ANTMASK` antenna selection.

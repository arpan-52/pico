#!/bin/bash -i
#
# PICO v1 node-side job. Runs on the compute node selected by the monitor.
#   pico_run_job.sh JOB_ID OUTPUT_DIR CONFIG_FILE BURST_MJD JOB_LOG PICO_IMAGE_BIN
#
# Steps:
#   1. copy svfits_par.txt + antsamp.hdr into the working dir
#   2. append band-derived imaging params to the par (pico_image extends svfits_par.txt)
#   3. run pico_image -> burst_<mjd>-{dirty,psf,image,residual,model}.fits
#   4. plot burst_<mjd>.png from the restored -image.fits

# ===== deploy paths + environment (hardcoded for the cluster) =====
TEST_IMAGING="/lustre_archive/apps/astrosoft/test_imaging"
PICO_HOME="$TEST_IMAGING/picov1"
ANTSAMP_SRC="$TEST_IMAGING/antsamp.hdr"

# Bring up the same micromamba 'pico' env used by the rest of the pipeline so
# pico_image, astropy and matplotlib all resolve to one interpreter. start.sh
# initialises micromamba; the shell hook makes 'activate' work over non-login ssh.
source "$TEST_IMAGING/start.sh"
eval "$(micromamba shell hook --shell bash 2>/dev/null)" || true
micromamba activate pico
export CASACORE_DATA="$TEST_IMAGING/data"
PICO_PY="$(command -v python3)"     # python3 from the activated 'pico' env

JOB_ID="$1"
OUTPUT_DIR="$2"
CONFIG_FILE="$3"
BURST_MJD="$4"
JOB_LOG="$5"
PICO_IMAGE_BIN="$6"

log_msg() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$JOB_LOG"; }
trap 'log_msg "[FAILED] Job error on $(hostname)"; exit 1' ERR

cd "$OUTPUT_DIR"

CONFIG_BASENAME=$(basename "$CONFIG_FILE")
[ -f "$CONFIG_BASENAME" ] || cp "$CONFIG_FILE" ./

# Resolve antsamp.hdr to an absolute path and copy a local provenance copy.
# Search order: deploy dir, OUTPUT_DIR, the raw-data PATH from the config.
# Stop at whitespace or '!' (inline comment marker in svfits par format).
DATA_PATH=$(grep -Po '^PATH\s+\K[^!\s]+' "$CONFIG_BASENAME" | head -1)
ANTSAMP_ABS=""
for cand in "$ANTSAMP_SRC" "$OUTPUT_DIR/antsamp.hdr" "${DATA_PATH%/}/antsamp.hdr"; do
    if [ -n "$cand" ] && [ -f "$cand" ]; then ANTSAMP_ABS="$(readlink -f "$cand")"; break; fi
done
[ -n "$ANTSAMP_ABS" ] || { log_msg "[ERROR] antsamp.hdr not found (looked in: $ANTSAMP_SRC, $OUTPUT_DIR, $DATA_PATH)"; exit 1; }
[ "$ANTSAMP_ABS" = "$(readlink -f ./antsamp.hdr 2>/dev/null)" ] || cp "$ANTSAMP_ABS" ./antsamp.hdr
log_msg "[ANTSAMP] using $ANTSAMP_ABS"

MJD_TAG=$(echo "$BURST_MJD" | tr . _)
BASE="burst_${MJD_TAG}"

# ===== imaging params (band-derived; currently uniform across all GMRT bands) =====
# Band centre from FREQ_SET freq0:freq1:nchan -> used for logging / future tuning.
FREQ_SET=$(grep -Po 'FREQ_SET\s+\K[\d.:eE+-]+' "$CONFIG_BASENAME" | head -1)
if [ -n "$FREQ_SET" ]; then
    F0=${FREQ_SET%%:*}; F1=$(echo "$FREQ_SET" | cut -d: -f2)
    BAND_MHZ=$("$PICO_PY" -c "print(int((float('$F0')+float('$F1'))/2/1e6))")
else
    BAND_MHZ=1400
fi

NPIX=4096
CELLSIZE_ASEC=1.0
WEIGHTING=natural
CLEAN_ITER=1000
CLEAN_GAIN=0.1
CLEAN_THRESH=10.0      # sigma (pico_image stops at nsigma * MAD noise)
REF_FREQ_HZ=0          # 0 -> midband
log_msg "[PARAMS] band_centre=${BAND_MHZ}MHz NPIX=$NPIX CELL=${CELLSIZE_ASEC}\" thresh=${CLEAN_THRESH}sigma"

# Append imaging block once (idempotent on re-run).
PAR="$CONFIG_BASENAME"
if ! grep -q '^NPIX' "$PAR"; then
    cat >> "$PAR" << EOF

* ===== pico_image imaging params (appended by pico_run_job.sh) =====
NPIX            $NPIX
CELLSIZE_ASEC   $CELLSIZE_ASEC
WEIGHTING       $WEIGHTING
CLEAN_ITER      $CLEAN_ITER
CLEAN_GAIN      $CLEAN_GAIN
CLEAN_THRESH    $CLEAN_THRESH
REF_FREQ_HZ     $REF_FREQ_HZ
ANT_HDR         $ANTSAMP_ABS
EOF
    log_msg "[PARAMS] appended imaging block to $PAR"
fi

# Force exactly one IMAGE_FITS = ours. The original svfits_par.txt ships its own
# (e.g. PICO_IMAGE.FITS); pico_image takes the FIRST IMAGE_FITS, so a leftover
# line would name the output wrong. Strip all, then write ours.
sed -i '/^IMAGE_FITS/d' "$PAR"
echo "IMAGE_FITS      ${BASE}.FITS" >> "$PAR"

# ===== STAGE 1: pico_image =====
log_msg "[STAGE: pico_image] Starting imaging on $(hostname)"
[ -x "$PICO_IMAGE_BIN" ] || { log_msg "[ERROR] pico_image binary not found/executable: $PICO_IMAGE_BIN"; exit 1; }

start=$(date +%s)
"$PICO_IMAGE_BIN" -c "$PAR" -a "$ANTSAMP_ABS" >> "$JOB_LOG" 2>&1
dur=$(( $(date +%s) - start ))

if [ -f "${BASE}-image.fits" ]; then
    log_msg "[STAGE: pico_image] Completed in ${dur}s ($(ls -lh ${BASE}-image.fits | awk '{print $5}'))"
else
    log_msg "[ERROR] pico_image produced no ${BASE}-image.fits"
    exit 1
fi

# ===== STAGE 2: PNG =====
# Publication-quality 5 arcmin cutout around the peak: WCS axes, clean-beam
# ellipse, 1 arcmin scalebar, afmhot/zscale. NOTE: this heredoc is unquoted so
# the shell expands $OUTPUT_DIR / ${BASE} -- keep '$' out of the Python (no
# matplotlib mathtext); use plain/unicode labels instead.
FIELD_ARCMIN=5.0
log_msg "[STAGE: png] Locating positive peak, plotting ${FIELD_ARCMIN} arcmin zoom"
"$PICO_PY" << PNGEOF
import sys, os, glob
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.visualization import ImageNormalize, LinearStretch, ZScaleInterval
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

FIELD_ARCMIN = ${FIELD_ARCMIN}

matplotlib.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "dejavuserif",
    "axes.linewidth": 1.1, "font.size": 13,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.major.size": 6, "ytick.major.size": 6,
})

output_dir = "$OUTPUT_DIR"
fits_files = glob.glob(os.path.join(output_dir, "${BASE}-image.fits")) or glob.glob(os.path.join(output_dir, "*-image.fits"))
if not fits_files:
    print("[ERROR] No -image.fits found"); sys.exit(1)

fits_file = fits_files[0]
burst_name = os.path.basename(fits_file).split("-")[0]   # burst_<mjd>

with fits.open(fits_file) as hdul:
    hdr = hdul[0].header
    data = hdul[0].data
    data = data[0, 0] if data.ndim == 4 else data
    wcs = WCS(hdr)
    while wcs.naxis > 2:
        wcs = wcs.dropaxis(2)

pix_asec = abs(hdr["CDELT2"]) * 3600.0
bmaj = hdr.get("BMAJ"); bmin = hdr.get("BMIN"); bpa = hdr.get("BPA")
have_beam = None not in (bmaj, bmin, bpa)
if have_beam:
    bmaj *= 3600.0; bmin *= 3600.0

# Brightest positive pixel in the full image -> zoom target.
peak_idx = np.nanargmax(data)
py, px = np.unravel_index(peak_idx, data.shape)
peak_val = float(data[py, px])
peak_world = wcs.pixel_to_world(px, py)

half = max(1, int(round(FIELD_ARCMIN * 60.0 / pix_asec / 2)))
y1, y2 = max(0, py - half), min(data.shape[0], py + half)
x1, x2 = max(0, px - half), min(data.shape[1], px + half)
cut = data[y1:y2, x1:x2]
cut_wcs = wcs[y1:y2, x1:x2]

# Robust noise (MAD) for the colour scale / S-N.
v = cut[np.isfinite(cut)]
med = float(np.median(v)) if v.size else 0.0
rms = float(np.median(np.abs(v - med)) * 1.4826) if v.size else 1.0
snr = peak_val / rms if rms > 0 else 0.0
print("[PNG] peak=%.4g Jy/beam rms~%.4g SNR~%.0f field=%g arcmin" % (peak_val, rms, snr, FIELD_ARCMIN))

# Beam PA in the DISPLAY frame. FITS BPA is deg East-of-North; matplotlib's
# Ellipse angle is CCW from +x screen axis. Derive from the WCS (step one
# beam-major along PA, map back to pixels) so it is correct for any RA-flip /
# rotation, then apply a 90 deg flip (CASA vs astropy beam convention).
def beam_angle_display(w, x0, y0):
    c0 = w.pixel_to_world(x0, y0)
    pa = np.deg2rad(bpa)
    dN = np.cos(pa) * (bmaj / 3600.0)
    dE = np.sin(pa) * (bmaj / 3600.0)
    c1 = SkyCoord((c0.ra.deg + dE / np.cos(c0.dec.radian)) * u.deg,
                  (c0.dec.deg + dN) * u.deg)
    x1p, y1p = w.world_to_pixel(c1)
    return np.degrees(np.arctan2(y1p - y0, x1p - x0)) - 90.0

lo, hi = ZScaleInterval(contrast=0.15).get_limits(cut)
norm = ImageNormalize(vmin=lo, vmax=peak_val, stretch=LinearStretch())

fig = plt.figure(figsize=(8.5, 7.6), facecolor="white")
ax = fig.add_subplot(111, projection=cut_wcs)
im = ax.imshow(cut, origin="lower", cmap="afmhot", norm=norm, interpolation="nearest")

ext = cut.shape[0]
if have_beam:
    ang = float(beam_angle_display(wcs, px, py)) + 90.0
    bx, by = 0.085 * ext, 0.085 * ext
    ax.add_patch(Ellipse((bx, by), bmin / pix_asec, bmaj / pix_asec, angle=ang,
                         fc="white", ec="black", lw=1.0, alpha=0.9, zorder=6))
    ax.text(bx, by + 0.5 * bmaj / pix_asec + 6, "beam", color="white",
            ha="center", va="bottom", fontsize=9, zorder=6)

# 1 arcmin scalebar.
L = 60.0 / pix_asec
x0 = 0.07 * ext; yb = 0.93 * ext
ax.plot([x0, x0 + L], [yb, yb], color="white", lw=2.5, zorder=6, solid_capstyle="butt")
ax.text(x0 + L / 2, yb - 7, u"1'", color="white", ha="center", va="top", fontsize=10, zorder=6)

ax.coords[0].set_axislabel("Right Ascension (J2000)", fontsize=13)
ax.coords[1].set_axislabel("Declination (J2000)", fontsize=13)
ax.coords[0].set_major_formatter("hh:mm:ss")
ax.coords[1].set_major_formatter("dd:mm:ss")
ax.coords[0].set_ticks(number=5); ax.coords[1].set_ticks(number=5)
ax.set_title("%s  |  peak %.3g Jy/beam,  S/N ~ %.0f  |  %g arcmin" %
             (burst_name, peak_val, snr, FIELD_ARCMIN), fontsize=13, pad=12)

cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
cb.set_label(u"Flux density (Jy beam⁻¹)", fontsize=12)
cb.ax.tick_params(labelsize=10)

fig.tight_layout()
png = os.path.join(output_dir, burst_name + ".png")
fig.savefig(png, dpi=200, bbox_inches="tight", facecolor="white")
plt.close(fig)
print("[PNG] saved " + png)
PNGEOF

log_msg "[STAGE: png] Done"
log_msg "[SUCCESS] Pipeline completed on $(hostname)"

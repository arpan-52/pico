#!/usr/bin/env python
# Plot the 8 dirty images dumped by pico_imgtest into a single PNG.
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os, sys

OUT = os.path.join(os.path.dirname(__file__), "build", "imgtest_out")
if not os.path.isdir(OUT):
    OUT = "imgtest_out"

titles = {
    1: "1. point src @ centre\n(imaging sanity)",
    2: "2. point src off-centre\n(geometry sanity)",
    3: "3. pure noise\n(blank)",
    4: "4. noise + SPIKY rfi + CLIP\n(blank recovered)",
    5: "5. noise + COHERENT-DC rfi + CLIP\n(FAILS: central source)",
    6: "6. coherent-DC rfi + DC-SUB\n(FIX: blank)",
    7: "7. burst(S=8) + DC rfi, CLIP only\n(source buried)",
    8: "8. burst(S=8) + DC rfi, DC-SUB\n(source recovered)",
}

def load(i):
    p = os.path.join(OUT, f"scen{i:02d}.f64")
    with open(p, "rb") as f:
        n = int(np.fromfile(f, dtype=np.int32, count=1)[0])
        img = np.fromfile(f, dtype=np.float64, count=n*n).reshape(n, n)
    return img

fig, axes = plt.subplots(2, 4, figsize=(20, 10))
for i in range(1, 9):
    ax = axes[(i-1)//4, (i-1)%4]
    img = load(i)
    n = img.shape[0]
    # robust scale so noise scenarios aren't washed out by a bright RFI panel
    med = np.median(img)
    mad = 1.4826*np.median(np.abs(img-med))
    vmax = max(img.max(), med + 5*mad)
    vmin = min(img.min(), med - 5*mad)
    im = ax.imshow(img, origin="lower", cmap="inferno", vmin=vmin, vmax=vmax)
    ax.set_title(titles[i], fontsize=11)
    ax.axhline(n//2, color="cyan", lw=0.4, alpha=0.4)
    ax.axvline(n//2, color="cyan", lw=0.4, alpha=0.4)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.ax.tick_params(labelsize=8)
    ax.text(0.03, 0.97, f"peak={img.max():.3g}\nmin={img.min():.3g}\ncentre={img[n//2,n//2]:.3g}",
            transform=ax.transAxes, va="top", ha="left", color="white", fontsize=8,
            bbox=dict(facecolor="black", alpha=0.5, pad=2))

fig.suptitle("pico imaging-path test: clip handles spiky RFI, but only DC-subtraction "
             "removes coherent RFI (and preserves a transient source)", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.97])
outpng = os.path.join(os.path.dirname(__file__), "imgtest.png")
fig.savefig(outpng, dpi=110)
print("wrote", outpng)

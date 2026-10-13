# -*- coding: utf-8 -*-
r"""Automated retrieval check over every released record.

Each released 90-degree radiograph is compared with parallel-beam projections of its own
recorded crop window and of 19 distractor crops drawn from the same subset, scoring
gradient correlation maximised over translation, rotation and magnification, with every
candidate scored identically.

The matched score of a record depends only on that record, so a record's own score is the
same whether the distractors come from a sample or from the whole release; only its rank
depends on the distractor pool.

Method, in the order it is applied:

  radiograph   DICOM -> 1-99.5 percentile window -> 1/6 scale (139 x 114) -> gradient magnitude
  crop         NIfTI -> mean along the depth axis -> a 25 x 30 mm window
                       (matching the 24.1 x 29.2 mm radiograph field) -> transpose
                       -> 8 (rotation, magnification) variants -> 1/6 scale -> gradient magnitude
  similarity   normalised cross-correlation by FFT, zero-padded to twice the image size,
               maximised over integer shifts inside a +/-27 x +/-20 mask (1/6 scale,
               equivalent to +/-162 x +/-120 full-resolution pixels, or +/-5.67 x +/-4.20 mm)
  rank         1 + the number of distractors scoring higher than the record's own crop

Memory: the naive form memoises every candidate variant transform, which is bounded for a
sample but reaches roughly 9 GB over 2,229 records. This version inverts the loop instead:
each record is turned into its eight transforms once, scored against every radiograph that
needs it, then released. Same work, bounded peak memory.

Output, in --out: correspondence_check.csv, correspondence_check.npz, correspondence_check.txt
"""
import argparse
import csv
import io
import os
import time

import numpy as np
import nibabel as nib
import pydicom
from PIL import Image

_ap = argparse.ArgumentParser(
    description="Automated retrieval check over the released records.")
_ap.add_argument("--dataset", default="physionet-paired-radiograph-cbct",
                 help="released dataset directory, relative to the current working "
                      "directory (default: %(default)s)")
_ap.add_argument("--out", default=".", help="output directory (default: %(default)s)")
_ap.add_argument("--distractors", type=int, default=19,
                 help="distractors per record (default: %(default)s)")
_ap.add_argument("--seed", type=int, default=20260920,
                 help="distractor sampling seed (default: %(default)s)")
_args = _ap.parse_args()

D = _args.dataset
OUTDIR = _args.out
os.makedirs(OUTDIR, exist_ok=True)
OUTP = os.path.join(OUTDIR, "correspondence_check.txt")
NPZ = os.path.join(OUTDIR, "correspondence_check.npz")
CSVP = os.path.join(OUTDIR, "correspondence_check.csv")

# ── protocol constants ──────────────────────────────────────────────────────
SC = 6                                   # working scale: 1/6
SHIFT = (27, 20)                         # translation search, in 1/6-scale pixels
ND = _args.distractors
SEED_DISTRACT = _args.seed
SEED_CI = 20260918
N_BOOT = 4000
VARIANTS = [(0, 1.00), (-6, 1.00), (6, 1.00), (-3, 1.00), (3, 1.00),
            (0, 1.08), (-6, 1.08), (6, 1.08)]
CROP_WINDOW = (slice(10, 110), slice(20, 140))   # 25 x 30 mm inside a 30 x 40 mm crop

out = io.open(OUTP, "w", encoding="utf-8")


def P(*a):
    s = " ".join(str(x) for x in a)
    out.write(s + "\n")
    out.flush()
    print(s, flush=True)


def win(a, lo=1.0, hi=99.5):
    """Percentile window to [0, 1]; robust to hot and cold pixels."""
    l, h = np.percentile(a, lo), np.percentile(a, hi)
    return np.clip((a - l) / max(h - l, 1e-6), 0, 1)


def grad_mag(x):
    """Central-difference gradient magnitude, normalised. Drops absolute intensity, so
    the two modalities become comparable."""
    gx = np.zeros_like(x)
    gy = np.zeros_like(x)
    gx[:, 1:-1] = x[:, 2:] - x[:, :-2]
    gy[1:-1, :] = x[2:, :] - x[:-2, :]
    g = np.sqrt(gx * gx + gy * gy)
    return g / (g.max() + 1e-9)


# ── the release, labelled with the subset names used for distractor pooling ──
meta = {r["sample_id"]: r for r in csv.DictReader(
    io.open(os.path.join(D, "metadata", "metadata.csv"), encoding="utf-8-sig"))}
SUBSET_OF_REGION = {"limb_unspecified": "limb", "maxillary": "maxilla",
                    "mandibular": "mandible"}


def subset_of(m):
    return "mouse" if m["species"] == "mouse" else SUBSET_OF_REGION[m["anatomical_region"]]


ids = sorted(meta)
subset = {s: subset_of(meta[s]) for s in ids}
P("records: %d" % len(ids))
P("subsets: %s" % {g: sum(1 for s in ids if subset[s] == g)
                   for g in ("limb", "maxilla", "mandible", "mouse")})

t0 = time.time()
PROJ, FR = {}, {}
RSIZE = (688 // SC, 834 // SC)
for k, sid in enumerate(ids):
    m = meta[sid]
    base = os.path.join(D, "data", m["species"], m["data_group_id"],
                        m["source_volume_id"], sid)
    vol = np.asanyarray(nib.load(
        os.path.join(base, sid + "__cbct_crop.nii.gz")).dataobj).astype(np.float32)
    PROJ[sid] = Image.fromarray(
        (win(vol.mean(axis=2)[CROP_WINDOW].T) * 255).astype(np.uint8))
    im = pydicom.dcmread(
        os.path.join(base, sid + "__xray__a090.dcm")).pixel_array.astype(np.float32)
    FR[sid] = grad_mag(np.asarray(
        Image.fromarray((win(im) * 255).astype(np.uint8)).resize(
            (RSIZE[0], RSIZE[1]), Image.BILINEAR), np.float32) / 255.0)
    if (k + 1) % 250 == 0:
        P("  loaded %d/%d  %.0fs" % (k + 1, len(ids), time.time() - t0))

shape = (RSIZE[1] * 2, RSIZE[0] * 2)
f0 = np.fft.fftfreq(shape[0], 1.0 / shape[0]).astype(int)
f1 = np.fft.fftfreq(shape[1], 1.0 / shape[1]).astype(int)
MASK = (np.abs(f0)[:, None] <= SHIFT[0]) & (np.abs(f1)[None, :] <= SHIFT[1])


def fft_e(g):
    a = g - g.mean()
    return np.fft.rfft2(a, s=shape), float((a * a).sum())


P("transforming radiographs")
FRF = {}
for k, s in enumerate(ids):
    FRF[s] = fft_e(FR[s])
    if (k + 1) % 500 == 0:
        P("  fft %d/%d  %.0fs" % (k + 1, len(ids), time.time() - t0))
del FR


def variant_ffts(sid):
    """The eight (rotation, magnification) transforms of one crop projection."""
    img0 = PROJ[sid]
    out = []
    for rot, mag in VARIANTS:
        img = img0
        if mag > 1.0:
            w, h = img.size
            cw, ch = int(round(w / mag)), int(round(h / mag))
            img = img.crop(((w - cw) // 2, (h - ch) // 2,
                            (w + cw) // 2, (h + ch) // 2))
        if rot:
            img = img.rotate(rot, resample=Image.BILINEAR, fillcolor=0)
        g = grad_mag(np.asarray(
            img.resize((RSIZE[0], RSIZE[1]), Image.BILINEAR), np.float32) / 255.0)
        out.append(fft_e(g))
    return out


rngd = np.random.default_rng(SEED_DISTRACT)
distr_of = {}
for sid in ids:
    pool = [o for o in ids if o != sid and subset[o] == subset[sid]]
    assert len(pool) >= ND, (sid, len(pool))
    distr_of[sid] = list(rngd.choice(pool, ND, replace=False))

needed_by = {s: [] for s in ids}
for sid in ids:
    needed_by[sid].append(sid)
    for c in distr_of[sid]:
        needed_by[c].append(sid)

P("scoring: %d records x %d candidates" % (len(ids), ND + 1))
best = {s: [] for s in ids}
own = {}
for k, c in enumerate(ids):
    vf = variant_ffts(c)
    for sid in needed_by[c]:
        fr, nr = FRF[sid]
        b = -9.0
        for fd, nd in vf:
            cc = np.fft.irfft2(fr * np.conj(fd), s=shape) / np.sqrt(nr * nd)
            v = float(np.where(MASK, cc, -np.inf).max())
            if v > b:
                b = v
        if sid == c:
            own[sid] = b
        else:
            best[sid].append((b, c))
    if (k + 1) % 100 == 0:
        P("  scored %d/%d  %.0fs" % (k + 1, len(ids), time.time() - t0))

P("scoring done %.0fs" % (time.time() - t0))

per_record = []
for sid in ids:
    ms = [v for v, _c in best[sid]]
    rk = 1 + sum(1 for v in ms if v > own[sid])
    top = max(range(len(ms)), key=lambda i: ms[i])
    per_record.append((sid, subset[sid], rk, own[sid], ms[top], best[sid][top][1]))

matched = np.array([r[3] for r in per_record])
mismatched = np.array([v for sid in ids for v, _c in best[sid]])
ranks = np.array([r[2] for r in per_record])


def auc(pos, neg):
    a = np.concatenate([pos, neg])
    r = np.argsort(np.argsort(a)) + 1
    return (r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


A = auc(matched, mismatched)
boot = ranks[np.random.default_rng(SEED_CI).integers(0, len(ranks), (N_BOOT, len(ranks)))]
t1lo, t1hi = np.percentile(100.0 * (boot <= 1).mean(axis=1), [2.5, 97.5])
t3lo, t3hi = np.percentile(100.0 * (boot <= 3).mean(axis=1), [2.5, 97.5])

P("")
P("=" * 84)
P("automated retrieval check, all %d records; distractors from the whole release" % len(ids))
P("=" * 84)
P("  first   %d (%.1f%%)  [%.1f-%.1f]   top-3 %.1f%% [%.1f-%.1f]   mean rank %.2f"
  % (int((ranks == 1).sum()), 100 * (ranks == 1).mean(), t1lo, t1hi,
     100 * (ranks <= 3).mean(), t3lo, t3hi, ranks.mean()))
P("  chance  5.0%% first, mean rank %.1f" % ((1 + ND + 1) / 2))
P("  matched %d   mismatched %d   AUC %.3f" % (len(matched), len(mismatched), A))
for g in ("limb", "maxilla", "mandible", "mouse"):
    m = np.array([subset[s] for s in ids]) == g
    P("    %-9s n=%4d  first %5.1f%%  top-3 %5.1f%%  mean rank %5.2f"
      % (g, int(m.sum()), 100 * (ranks[m] == 1).mean(),
         100 * (ranks[m] <= 3).mean(), ranks[m].mean()))
P("")
P("  This is a description of the release, not a held-out estimate: the search variants")
P("  and the shift window are fixed constants, so nothing is fitted to these records.")

np.savez(NPZ, matched=matched, mismatched=mismatched, ranks=ranks,
         ids=np.array(ids), subset=np.array([subset[s] for s in ids]))
with io.open(CSVP, "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f)
    w.writerow(["sample_id", "subset", "rank", "n_candidates", "score_matched",
                "score_best_other", "margin", "best_other_id"])
    for sid, sp, rk, o, b, cid in per_record:
        w.writerow([sid, sp, rk, ND + 1, "%.4f" % o, "%.4f" % b,
                    "%.4f" % (o - b), cid])
P("\nwrote %s\nwrote %s\ntotal %.0fs" % (CSVP, NPZ, time.time() - t0))
out.close()

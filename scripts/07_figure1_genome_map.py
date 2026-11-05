#!/usr/bin/env python3
"""07_figure1_genome_map.py — Figure 1: circular map of the three CJM-SD1 replicons.

Input : a PGAP GenBank flat file. Works with EITHER
          - the local submission file (record ids "chromosome", "pCJMSD1_1", "pCJMSD1_2"), or
          - the NCBI download for GCA_062071506.1 (ids CM205908.1, JCFRYO0100000xx.1)
        Replicons are resolved by length, so no id mapping is needed.
Output: figs/Fig1_genome_map.{png,svg}

Rings, outside in: forward CDS, reverse CDS, RNA features, GC skew, GC content.
The eps cluster, the plasmid polysaccharide genes and the two proviruses are
called out by name. Radii are scaled per replicon so the two 32 kb plasmids are
legible beside a 2.45 Mb chromosome; each panel therefore has its own scale bar.

Updated 2026-11-05: added the two proviruses found by geNomad (see
01_paper1_genome_and_EPS.md §4h). Previously the map showed no prophage content.
"""
import sys, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
from Bio import SeqIO

GB  = sys.argv[1] if len(sys.argv) > 1 else "CJM-SD1.bgpipe.output_2765522.gb"
OUT = sys.argv[2] if len(sys.argv) > 2 else "figs/Fig1_genome_map"

# features to label, keyed by CANONICAL replicon name: (start, end, label, colour)
HIGHLIGHT = {
    "chromosome": [
        (713383,  728136,  "eps cluster\nADATZ4_03440–03505", "#c2410c"),
        # geNomad 2026-11-05 — both complete, with holin + GH25 endolysin
        (1015613, 1058411, "provirus 1\n42.8 kb",             "#86198f"),
        (1110395, 1159112, "provirus 2\n48.7 kb",             "#86198f"),
    ],
    "pCJMSD1_2":  [(17754, 20250,  "polysaccharide genes\n12305–12310",  "#c2410c"),
                   (16156, 17331,  "IS256", "#6b7280"),
                   (20931, 21999,  "IS3",   "#6b7280"),
                   (22410, 23339,  "ISLpl1","#6b7280")],
    "pCJMSD1_1":  [],
}
FWD, REV, RNA, PHAGE = "#2563eb", "#7c3aed", "#059669", "#86198f"


def resolve_replicons(recs):
    """Map canonical replicon names -> record ids present in THIS file.

    Resolved by length (chromosome 2.45 Mb > pCJMSD1_1 32,797 > pCJMSD1_2 32,245),
    so the same script works on the local submission file and on the NCBI download,
    whose record ids differ.
    """
    by_len = sorted(recs.values(), key=lambda r: -len(r.seq))
    if len(by_len) != 3:
        sys.exit(f"expected 3 replicons, found {len(by_len)}: {[r.id for r in by_len]}")
    canon = {"chromosome": by_len[0].id,
             "pCJMSD1_1":  by_len[1].id,
             "pCJMSD1_2":  by_len[2].id}
    exp = {"chromosome": 2449864, "pCJMSD1_1": 32797, "pCJMSD1_2": 32245}
    for name, rid in canon.items():
        got = len(recs[rid].seq)
        if got != exp[name]:
            print(f"  WARNING {name} ({rid}) is {got:,} bp, expected {exp[name]:,}")
        else:
            print(f"  {name:12s} -> {rid}  ({got:,} bp)")
    return canon


def gc_tracks(seq, nwin=720):
    L = len(seq); step = max(1, L // nwin)
    gc, skew, mid = [], [], []
    for i in range(0, L, step):
        w = seq[i:i+step].upper()
        g, c = w.count("G"), w.count("C")
        gc.append(100*(g+c)/max(1,len(w)))
        skew.append((g-c)/(g+c) if (g+c) else 0.0)
        mid.append(i + len(w)/2)
    return np.array(mid), np.array(gc), np.array(skew)


def draw(ax, rec, canon_name):
    L = len(rec.seq)
    ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
    ax.set_ylim(0, 10); ax.axis("off"); ax.set_rorigin(0)
    deg = lambda p: 360.0 * p / L          # position -> POLAR degrees, clockwise from north

    # Wedge() is drawn through ax.transData._b, the plain Cartesian transform under the
    # polar axes. It therefore ignores set_theta_zero_location("N") and
    # set_theta_direction(-1), which DO apply to annotate() and to the tick marks below.
    # Without this conversion every wedge is mirrored and rotated 90 deg away from its
    # own label and from the kb scale. wdeg() maps polar degrees to Wedge's convention
    # (counterclockwise from east); Wedge sweeps theta1 -> theta2 counterclockwise, so
    # the END of a feature supplies theta1 and the START supplies theta2.
    wdeg = lambda p: 90.0 - deg(p)

    # CDS rings
    for f in rec.features:
        if f.type not in ("CDS", "rRNA", "tRNA", "tmRNA", "ncRNA"):
            continue
        s, e = int(f.location.start), int(f.location.end)
        if f.type == "CDS":
            r0, col = (8.3, FWD) if f.location.strand == 1 else (7.5, REV)
            h = 0.7
        else:
            r0, col, h = 6.8, RNA, 0.55
        span = max(deg(e) - deg(s), 0.06)
        t2 = wdeg(s)
        ax.add_patch(Wedge((0, 0), r0 + h, t2 - span, t2,
                           width=h, transform=ax.transData._b,
                           facecolor=col, edgecolor="none"))

    mid, gc, skew = gc_tracks(str(rec.seq))
    th = 2*np.pi*mid/L
    g = (gc - gc.mean()) / max(1e-9, np.abs(gc - gc.mean()).max())
    ax.fill_between(th, 5.6, 5.6 + g*0.7, color="#334155", lw=0, step="mid")
    sk = skew / max(1e-9, np.abs(skew).max())
    ax.fill_between(th, 4.4, 4.4 + sk*0.7, where=sk >= 0, color="#0891b2", lw=0, step="mid")
    ax.fill_between(th, 4.4, 4.4 + sk*0.7, where=sk < 0,  color="#be123c", lw=0, step="mid")
    for r in (5.6, 4.4):
        ax.plot(np.linspace(0, 2*np.pi, 400), [r]*400, color="#94a3b8", lw=0.4)

    # highlights. Small features get a minimum angular width so they are visible
    # at chromosome scale; the label carries the true coordinates.
    MINDEG = 2.2
    for k, (s, e, lab, col) in enumerate(HIGHLIGHT.get(canon_name, [])):
        d0, d1 = deg(s), deg(e)
        if d1 - d0 < MINDEG:
            mid_d = (d0 + d1) / 2; d0, d1 = mid_d - MINDEG/2, mid_d + MINDEG/2
        ax.add_patch(Wedge((0, 0), 9.5, 90.0 - d1, 90.0 - d0, width=1.0,
                           transform=ax.transData._b, facecolor=col, edgecolor="none"))
        a = math.radians((d0 + d1) / 2)
        rr = 13.2 + 1.5 * (k % 2)          # stagger so adjacent callouts clear
        ax.annotate(lab, xy=(a, 9.6), xytext=(a, rr),
                    ha="center", va="center", fontsize=7.2, color=col, weight="bold",
                    annotation_clip=False,
                    arrowprops=dict(arrowstyle="-", color=col, lw=0.8,
                                    shrinkA=0, shrinkB=2))
    # ticks
    span = 500_000 if L > 1e6 else 10_000
    for p in range(0, L, span):
        a = math.radians(deg(p))
        ax.plot([a, a], [9.5, 9.9], color="#475569", lw=0.7)
        ax.text(a, 10.6, f"{p//1000:,} kb" if p else "0",
                ha="center", va="center", fontsize=6.2, color="#475569")

    seq = str(rec.seq).upper()
    gcpc = 100*(seq.count("G")+seq.count("C"))/len(seq)
    ax.text(0.5, 0.5, f"{canon_name}\n{L:,} bp\nGC {gcpc:.1f}%",
            ha="center", va="center", fontsize=9.2, weight="bold",
            color="#0f172a", transform=ax.transAxes, zorder=5)


recs = {r.id: r for r in SeqIO.parse(GB, "genbank")}
print(f"parsed {len(recs)} records from {GB}")
canon = resolve_replicons(recs)

order = ["chromosome", "pCJMSD1_1", "pCJMSD1_2"]
fig = plt.figure(figsize=(14.5, 6.4), dpi=200)
widths = [1.5, 1, 1]
for i, (name, w) in enumerate(zip(order, widths)):
    ax = fig.add_subplot(1, 3, i+1, projection="polar")
    draw(ax, recs[canon[name]], name)

handles = [plt.Line2D([], [], color=c, lw=6, label=l) for c, l in
           [(FWD, "CDS, forward"), (REV, "CDS, reverse"), (RNA, "rRNA / tRNA / ncRNA"),
            ("#334155", "GC content"), ("#0891b2", "GC skew +"), ("#be123c", "GC skew −"),
            ("#c2410c", "polysaccharide loci"), (PHAGE, "provirus"),
            ("#6b7280", "IS elements")]]
fig.legend(handles=handles, loc="lower center", ncol=9, frameon=False, fontsize=7.4,
           bbox_to_anchor=(0.5, -0.01))
fig.suptitle("Figure 1 — $\\it{Levilactobacillus\\ brevis}$ CJM-SD1, three closed replicons "
             "(PGAP annotation, GenBank JCFRYO010000000)", fontsize=10.5, y=1.0)
fig.subplots_adjust(left=0.04, right=0.96, top=0.86, bottom=0.12, wspace=0.42)
for ext in ("png", "svg"):
    fig.savefig(f"{OUT}.{ext}", bbox_inches="tight", facecolor="white")
print(f"wrote {OUT}.png and {OUT}.svg")

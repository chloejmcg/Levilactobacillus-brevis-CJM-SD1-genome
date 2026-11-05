#!/usr/bin/env python3
"""07b_figure2_eps_cluster.py — Figure 2: the eps locus, its precursor route,
and the plasmid-borne polysaccharide genes.

A  chromosomal eps cluster, ADATZ4_03440-03505, 14 CDS, Wzx/Wzy assignments
B  nucleotide-sugar precursor route feeding the cluster, with the alpha-PGM chokepoint
C  pCJMSD1_2 polysaccharide genes in their IS context

Input : a PGAP GenBank flat file. Works with EITHER the local submission file
        (record ids "chromosome"/"pCJMSD1_2") or the NCBI download for
        GCA_062071506.1 (CM205908.1 / JCFRYO0100000xx.1) — replicons are
        resolved by length.
Output: figs/Fig2_eps_cluster.{png,svg}

Updated 2026-11-05 (see 01_paper1_genome_and_EPS.md §3a, §4h):
  - Panel A no longer claims "single operon". The transcript arrow is now solid
    across the compact core (03445-03495, gaps <=59 bp, two overlapping ORFs)
    and dashed over the loosely-attached flanks (252 / 73 / 114 bp gaps).
    Transcript structure was not determined; this is a CONTENT fix, not styling.
  - Panel B: galactose confirmed and now IN the carbon panel (complete Leloir
    operon, GalT 1.2e-34); fructose DROPPED (mtlD present, ADATZ4_11390).
"""
import sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, FancyBboxPatch, FancyArrowPatch
from Bio import SeqIO

GB  = sys.argv[1] if len(sys.argv) > 1 else "CJM-SD1.bgpipe.output_2765522.gb"
OUT = sys.argv[2] if len(sys.argv) > 2 else "figs/Fig2_eps_cluster"

ROLE = {  # locus_tag -> (short label, functional class)
 "ADATZ4_03440": ("lcp",      "export"),
 "ADATZ4_03445": ("yveK",     "regulation"),
 "ADATZ4_03450": ("cpsD",     "regulation"),
 "ADATZ4_03455": ("cpsB",     "regulation"),
 "ADATZ4_03460": ("priming\nGT", "priming"),
 "ADATZ4_03465": ("gt",       "gt"),
 "ADATZ4_03470": ("gt2",      "gt"),
 "ADATZ4_03475": ("gt",       "gt"),
 "ADATZ4_03480": ("gro-P\ntransf.", "substituent"),
 "ADATZ4_03485": ("glf",      "precursor"),
 "ADATZ4_03490": ("wzx",      "transport"),
 "ADATZ4_03495": ("wzy",      "transport"),
 "ADATZ4_03500": ("lcp",      "export"),
 "ADATZ4_03505": ("rmlB",     "orphan"),
}
COL = {"regulation":"#7c3aed","priming":"#c2410c","gt":"#2563eb","substituent":"#059669",
       "precursor":"#db2777","transport":"#0891b2","export":"#64748b","orphan":"#a16207"}
PLASMID = {"ADATZ4_12295":("IS256","is"),"ADATZ4_12300":("hyp","other"),
           "ADATZ4_12305":("LPS biosynth.\nprotein","cargo"),
           "ADATZ4_12310":("pyruvyl\ntransferase","cargo"),
           "ADATZ4_12315":("tnp","is"),"ADATZ4_12320":("IS3","is"),
           "ADATZ4_12325":("ISLpl1","is")}
COL_P = {"is":"#6b7280","cargo":"#c2410c","other":"#cbd5e1"}

# Panel A transcript annotation. Coordinates are PGAP-rotated chromosome.
CORE_START, CORE_END = 714691, 725864      # 03445 start -> 03495 end, gaps <= 59 bp
FLANK_L_START        = 713383              # 03440 start, 252 bp gap to the core
FLANK_R_END          = 728136              # 03505 end, 73 and 114 bp gaps

recs = {r.id: r for r in SeqIO.parse(GB, "genbank")}
_by_len = sorted(recs.values(), key=lambda r: -len(r.seq))
if len(_by_len) != 3:
    sys.exit(f"expected 3 replicons, found {len(_by_len)}")
CANON = {"chromosome": _by_len[0].id, "pCJMSD1_1": _by_len[1].id, "pCJMSD1_2": _by_len[2].id}
print("  " + "  ".join(f"{k}->{v}" for k, v in CANON.items()))

def grab(canon_name, keys):
    out=[]
    for f in recs[CANON[canon_name]].features:
        if f.type!="CDS": continue
        lt=f.qualifiers.get("locus_tag",[""])[0]
        if lt in keys:
            out.append((lt,int(f.location.start)+1,int(f.location.end),f.location.strand,
                        f.qualifiers.get("product",[""])[0]))
    return sorted(out,key=lambda x:x[1])

cluster=grab("chromosome",ROLE); plas=grab("pCJMSD1_2",PLASMID)
if len(cluster)!=14: print(f"  WARNING: {len(cluster)} eps CDS found, expected 14")

fig=plt.figure(figsize=(14,8.8),dpi=200)
gs=fig.add_gridspec(3,1,height_ratios=[1.2,1.5,0.95],hspace=0.48)

def draw_genes(ax, genes, roles, colmap, y=0.0, h=0.30, fs=7.0):
    x0=min(g[1] for g in genes); x1=max(g[2] for g in genes)
    for lt,s,e,strand,prod in genes:
        lab,cls=roles[lt]; col=colmap[cls]
        L=e-s; head=min(L*0.42,(x1-x0)*0.016)
        xs,dx=(s,L-head) if strand==1 else (e,-(L-head))
        ax.add_patch(FancyArrow(xs,y,dx,0,width=h,head_width=h*1.75,head_length=head,
                                length_includes_head=True,facecolor=col,edgecolor="white",lw=0.7))
        ax.text((s+e)/2,y-h*2.3,lab,ha="center",va="top",fontsize=fs,color=col,weight="bold")
        ax.text((s+e)/2,y+h*2.0,lt.split("_")[1],ha="center",va="bottom",fontsize=5.8,color="#475569")
    ax.set_xlim(x0-(x1-x0)*0.035,x1+(x1-x0)*0.035); ax.set_ylim(-1.5,1.55)
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.set_yticks([]); ax.tick_params(axis="x",labelsize=7,colors="#475569",length=3)
    return x0,x1

# ---- A ----
axA=fig.add_subplot(gs[0])
x0,x1=draw_genes(axA,cluster,ROLE,COL)
axA.set_xlabel("position on chromosome (bp)",fontsize=7.5,color="#475569",labelpad=1)
axA.set_title("A   Chromosomal $\\it{eps}$ cluster — 14 co-oriented CDS, ADATZ4_03440–03505, "
              f"{x1-x0+1:,} bp",loc="left",fontsize=10,weight="bold",pad=16)

# Transcript annotation — see module docstring and paper §3a.
# Solid over the compact core; dashed over the loosely-attached flanks.
YT = 0.95
axA.annotate("", xy=(CORE_END, YT), xytext=(CORE_START, YT),
             arrowprops=dict(arrowstyle="-|>", color="#64748b", lw=1.4))
for xa, xb in ((FLANK_L_START, CORE_START), (CORE_END, FLANK_R_END)):
    axA.annotate("", xy=(xb, YT), xytext=(xa, YT),
                 arrowprops=dict(arrowstyle="-", color="#cbd5e1", lw=1.1,
                                 linestyle=(0, (3, 2))))
axA.text((CORE_START+CORE_END)/2, YT+0.12,
         "compact core — gaps ≤ 59 bp, two overlapping ORFs",
         ha="center", fontsize=6.8, color="#64748b")
# above the arrow, not below: at YT-0.38 this collided with the locus-tag row (y=+0.60)
axA.text(x1, YT+0.34, "transcript structure not determined",
         ha="right", fontsize=6.6, color="#94a3b8", style="italic")

h=[plt.Line2D([],[],color=c,lw=6,label=l) for l,c in
   [("chain-length regulation",COL["regulation"]),("priming GT",COL["priming"]),
    ("glycosyltransferase",COL["gt"]),("Gro-P substituent",COL["substituent"]),
    ("Gal$\\it{f}$ precursor",COL["precursor"]),("flippase / polymerase",COL["transport"]),
    ("LCP export",COL["export"]),("orphan $\\it{rmlB}$",COL["orphan"])]]
_LEG=h  # drawn as a figure-level legend after all panels exist

# ---- B ----
axB=fig.add_subplot(gs[1]); axB.axis("off")
axB.set_xlim(0,100); axB.set_ylim(0,100)
axB.set_title("B   Nucleotide-sugar precursor route — every sugar funnels through one α-phosphoglucomutase",
              loc="left",fontsize=10,weight="bold",pad=4)
feeds=[("glucose","Glk 07275",             "#0f172a"),
       ("maltose","MalP 10295 → β-PGM 04105","#0f172a"),
       ("galactose","GalK 02220 → GalT 02230","#059669"),
       ("pentose","XylA/B 11175/11170 → Tkt/Tal — reverse PPP","#0f172a")]
for i,(sug,route,tc) in enumerate(feeds):
    yy=88-i*11
    axB.text(1,yy,sug,fontsize=8.2,weight="bold",color=tc,va="center")
    axB.text(13.5,yy,route,fontsize=6.9,color="#64748b" if tc=="#0f172a" else tc,va="center")
    axB.add_patch(FancyArrowPatch((33.5,yy),(34.5,70),arrowstyle="-",color="#cbd5e1",lw=0.9,
                                  connectionstyle="arc3,rad=0.18"))
axB.text(1,88-4*11,"sucrose",fontsize=8.2,weight="bold",color="#b91c1c",va="center")
axB.text(13.5,88-4*11,"no route — no hydrolase, no GH32, no PTS",fontsize=6.9,color="#b91c1c",va="center")

def box(x,y,w,hh,txt,fc,ec,fs=7.6,tc="white",bold=True):
    axB.add_patch(FancyBboxPatch((x,y-hh/2),w,hh,boxstyle="round,pad=0.6,rounding_size=1.6",
                                 facecolor=fc,edgecolor=ec,lw=1.1))
    axB.text(x+w/2,y,txt,ha="center",va="center",fontsize=fs,color=tc,
             weight="bold" if bold else "normal")
def arrow(x1,y1,x2,y2,lab="",col="#334155",dx=0,dy=5.5,ha="center"):
    axB.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=11,
                                  color=col,lw=1.3))
    if lab: axB.text((x1+x2)/2+dx,(y1+y2)/2+dy,lab,ha=ha,fontsize=6.6,color=col,style="italic")

box(36,70,13,11,"G6P","#475569","#475569")
arrow(49,70,57,70,"ADATZ4_02205\nα-PGM","#b91c1c",dy=7.0)
box(57,70,14,11,"α-G1P","#475569","#475569")
arrow(71,70,79,70,"GalU 08865")
box(79,70,16,11,"UDP-Glc","#2563eb","#2563eb")
arrow(87,64.5,87,50,"GalE 02225",dx=2.2,dy=0,ha="left")
box(79,44,16,11,"UDP-Gal","#db2777","#db2777")
arrow(79,44,66,44,"Glf 03485")
box(50,44,16,11,"UDP-Gal$\\it{f}$","#db2777","#db2777")
axB.add_patch(FancyArrowPatch((58,38.5),(58,26),arrowstyle="-|>",mutation_scale=11,color="#c2410c",lw=1.6))
box(42,20,32,11,"eps repeat unit  (ADATZ4_03460 →)","#c2410c","#c2410c",fs=7.8)
axB.text(53,94,"CHOKEPOINT — the only enzyme in the genome that can make α-G1P",
         ha="center",fontsize=7.1,color="#b91c1c",weight="bold")

# Galactose bypass: Leloir GalT releases Glc-1-P DOWNSTREAM of the mutase, so the
# arrow must pass UNDER the G6P box and the alpha-PGM step and land on alpha-G1P.
# Routed through the box before, which drew it crossing the two things it bypasses.
axB.add_patch(FancyArrowPatch((34.0,64.0),(60.0,64.3),arrowstyle="-|>",mutation_scale=11,
                              color="#059669",lw=1.5,linestyle=(0,(4,2)),
                              connectionstyle="arc3,rad=0.25"))
axB.text(72,55,"galactose bypass (Leloir):\nGlc-1-P downstream of α-PGM",
         fontsize=6.6,color="#059669",style="italic",ha="center")

axB.text(1,20,"no PFK, no FBA, no FBPase, no $\\it{glpX}$\n→ no gluconeogenic route from\n"
             "   triose-P back to hexose-P",
         fontsize=6.9,color="#b91c1c",va="center")
axB.text(1,6,"fructose dropped from the panel — $\\it{mtlD}$ (ADATZ4_11390) makes it an "
             "electron acceptor, confounding redox and energy state",
         fontsize=7.0,color="#a16207",va="center")

# ---- C ----
axC=fig.add_subplot(gs[2])
draw_genes(axC,plas,PLASMID,COL_P,h=0.26,fs=6.6)
axC.set_xlabel("position on pCJMSD1_2 (bp)",fontsize=7.5,color="#475569",labelpad=1)
axC.set_title("C   pCJMSD1_2 polysaccharide cargo, bracketed by three IS families",
              loc="left",fontsize=10,weight="bold",pad=10)
axC.text(0.995,0.04,"RepA (12270) is a pseudogene — p2 likely replicates $\\it{in\\ trans}$ off p1",
         transform=axC.transAxes,ha="right",fontsize=6.8,color="#64748b",style="italic")

fig.legend(handles=_LEG,loc="lower center",ncol=8,frameon=False,fontsize=6.9,
           bbox_to_anchor=(0.5,0.002),columnspacing=1.2,handlelength=1.3)
fig.suptitle("Figure 2 — $\\it{eps}$ locus architecture and the precursor route that feeds it",
             fontsize=11.5,y=0.985)
fig.subplots_adjust(bottom=0.075)
for ext in ("png","svg"):
    fig.savefig(f"{OUT}.{ext}",bbox_inches="tight",facecolor="white")
print(f"wrote {OUT}.png and {OUT}.svg")

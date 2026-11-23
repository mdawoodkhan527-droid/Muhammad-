from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "carbon_closure_lhs_sensitivity_corrected.csv"
OUT = Path(__file__).resolve().parent
rows = list(csv.DictReader(DATA.open(encoding="utf-8")))
values = np.array([float(r["predicted_delta13C_calcite_permil"]) for r in rows])
observed, spread = -13.7, 2.1
median, low, high = float(np.median(values)), float(values.min()), float(values.max())
within = int(np.sum(np.abs(values - observed) <= spread))

ink, muted, grid, blue, orange, pale = "#263238", "#56636B", "#D9E1E5", "#0072B2", "#D55E00", "#F4D9B8"
plt.rcParams.update({"font.family":"sans-serif", "font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
                     "font.size":9, "axes.labelcolor":ink, "text.color":ink,
                     "xtick.color":muted, "ytick.color":muted, "pdf.fonttype":42, "ps.fonttype":42})
fig, ax = plt.subplots(figsize=(7.4,3.1))
fig.subplots_adjust(left=.25,right=.98,bottom=.25,top=.78)
ax.axvspan(observed-spread, observed+spread, color=pale, alpha=.62, zorder=0)
ax.axvline(observed, color=orange, lw=1.5, zorder=2)
ax.hlines(1, low, high, color=blue, lw=3, zorder=3)
ax.vlines([low,high], .9, 1.1, color=blue, lw=1.4, zorder=3)
ax.scatter([median],[1], marker="D", s=54, color=blue, edgecolor="white", lw=.8, zorder=4)
ax.axvline(median, color=blue, ls=(0,(4,2.5)), lw=1.15, zorder=1)
ax.set_xlim(-25.2,-9.4); ax.set_ylim(-.48,1.58); ax.set_yticks([1,0])
ax.set_yticklabels(["Model design ensemble","Published concretion"],fontsize=8.4)
ax.set_xlabel("Carbonate δ¹³C (‰ VPDB)",fontsize=9,labelpad=7)
ax.set_title("Conditional carbon-isotope closure without fitting",fontsize=10.5,fontweight="bold",color=ink,pad=12)
ax.text(.985,.91,f"{within} of {len(values)} design points\nwithin published spread",transform=ax.transAxes,ha="right",va="top",fontsize=8.2,
        bbox=dict(boxstyle="round,pad=.35",facecolor="white",edgecolor=grid,lw=.8))
ax.text(.02,.91,f"Median: {median:.2f}‰\nRange: {low:.2f} to {high:.2f}‰",transform=ax.transAxes,ha="left",va="top",fontsize=7.8,
        bbox=dict(boxstyle="round,pad=.35",facecolor="white",edgecolor=grid,lw=.8))
handles=[Patch(facecolor=pale,edgecolor=orange,lw=.7,alpha=.75,label="Published mean ± 2.1‰ spread"),
         Line2D([0],[0],color=orange,lw=1.5,label="Published mean"),
         Line2D([0],[0],color=blue,ls=(0,(4,2.5)),lw=1.15,label="Model median"),
         Line2D([0],[0],color=blue,lw=3,label="Model range")]
ax.legend(handles=handles,loc="lower right",bbox_to_anchor=(.985,.04),ncol=2,frameon=True,facecolor="white",edgecolor=grid,framealpha=.96,fontsize=7.2,columnspacing=1.1,handlelength=1.8,borderpad=.45)
ax.grid(axis="x",color=grid,lw=.7,zorder=0); ax.set_axisbelow(True)
for s in ("top","right"): ax.spines[s].set_visible(False)
for s in ("left","bottom"): ax.spines[s].set_color("#7D8B91"); ax.spines[s].set_linewidth(.8)
ax.tick_params(axis="y",length=0,pad=8); ax.tick_params(axis="x",labelsize=8,width=.7)
fig.savefig(OUT/"carbon_closure_corrected.png",dpi=300,bbox_inches="tight",pad_inches=.08)
fig.savefig(OUT/"carbon_closure_corrected.pdf",bbox_inches="tight",pad_inches=.08)

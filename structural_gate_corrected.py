from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

OUT=Path(__file__).resolve().parent
INK="#263238"; MUTED="#56636B"; GRID="#D9E1E5"; BLUE="#0072B2"; ORANGE="#D55E00"; GREY="#9AA7AD"
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],"font.size":9,
                     "axes.labelcolor":INK,"text.color":INK,"xtick.color":MUTED,"ytick.color":MUTED,
                     "pdf.fonttype":42,"ps.fonttype":42})
fig=plt.figure(figsize=(7.4,3.25))
gs=fig.add_gridspec(1,2,width_ratios=(1.15,1.0),left=.11,right=.98,bottom=.24,top=.80,wspace=.28)
ax0=fig.add_subplot(gs[0,0]); ax1=fig.add_subplot(gs[0,1])
labels=["Elastic moduli /\nhomogenisation","Biot coefficient /\nloads","Mechanical boundary\nconditions"]
x=np.arange(len(labels)); ax0.axhline(0,color=GRID,lw=1,zorder=1)
ax0.scatter(x,np.zeros_like(x),s=48,color=BLUE,edgecolor="white",lw=.9,zorder=3)
ax0.set_xlim(-.45,len(labels)-.55); ax0.set_ylim(-.42,.42); ax0.set_xticks(x,labels,fontsize=7.6); ax0.set_yticks([0],["0"],fontsize=8)
ax0.set_ylabel("Sensitivity of selected observations",fontsize=8.5,labelpad=7)
ax0.set_title("(a) Exact structural zeros",loc="left",fontsize=10.2,fontweight="bold",pad=10,color=INK)
ax0.text(.5,.81,r"$\partial O/\partial\theta_{\mathrm{mech}} = 0$",transform=ax0.transAxes,ha="center",va="center",fontsize=10.2)
for s in ("top","right"): ax0.spines[s].set_visible(False)
for s in ("left","bottom"): ax0.spines[s].set_color("#7D8B91"); ax0.spines[s].set_linewidth(.8)
ax0.tick_params(axis="x",length=0,pad=7); ax0.tick_params(axis="y",length=0,pad=5)

effect=.4258; threshold=1.0; refinement=2.7783
ax1.axhspan(0,refinement,color=GREY,alpha=.15,zorder=0,label="Coarse-to-fine refinement envelope (2.778%)")
ax1.bar([0],[effect],width=.62,color=BLUE,alpha=.88,zorder=2)
ax1.axhline(threshold,color=ORANGE,ls=(0,(4,2.5)),lw=1.35,zorder=3,label="Predeclared threshold (1.0%)")
ax1.text(0,effect+.07,"0.426%",ha="center",va="bottom",fontsize=9.2,fontweight="bold",color=INK)
ax1.set_xlim(-.58,.58); ax1.set_ylim(0,3.2); ax1.set_xticks([0],["Stage 1 → Stage 2"],fontsize=8.2)
ax1.set_ylabel("Relative output difference (%)",fontsize=8.3,labelpad=7)
ax1.set_title("(b) Coupling comparison is unresolved",loc="left",fontsize=10.2,fontweight="bold",pad=10,color=INK)
ax1.legend(loc="upper left",bbox_to_anchor=(.02,.99),frameon=False,fontsize=6.9,handlelength=1.7)
ax1.grid(axis="y",color=GRID,lw=.65,zorder=0); ax1.set_axisbelow(True)
for s in ("top","right"): ax1.spines[s].set_visible(False)
for s in ("left","bottom"): ax1.spines[s].set_color("#7D8B91"); ax1.spines[s].set_linewidth(.8)
ax1.tick_params(axis="x",length=0,pad=7); ax1.tick_params(axis="y",labelsize=7.8,width=.7)
fig.savefig(OUT/"structural_gate_corrected.png",dpi=300,bbox_inches="tight",pad_inches=.08)
fig.savefig(OUT/"structural_gate_corrected.pdf",bbox_inches="tight",pad_inches=.08)

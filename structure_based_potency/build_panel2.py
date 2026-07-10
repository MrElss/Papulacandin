"""Merge ADME-v1 + Inductive Bio LogD for the rational-design panel; make CSV + figure."""
import csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# id, class, linkage, LogD, LogD_flag, ADME_perm, ADME_lipo, sol_conf, rationale
rows = [
 ("papB_native","native","di-ester",1.39,"OOD",0.0344,0.367,"high","native rigid polyene tail"),
 ("A1_caspo_ester","aliphatic","ester",2.81,"",0.318,0.854,"med","pos1 caspofungin branched-sat"),
 ("A2_palmitoyl_ester","aliphatic","ester",2.98,"",0.289,1.113,"med","pos1 palmitoyl straight-sat"),
 ("A3_onekink_ester","aliphatic","ester",2.76,"",0.259,0.723,"med","pos1 C16:1 one-kink control"),
 ("G_anidula_ester","aromatic","ester",3.22,"",0.129,1.229,"med","pos1 anidulafungin terphenyl (rigid)"),
 ("G_mica_ester","aromatic","ester",3.00,"",0.205,0.921,"med","pos1 micafungin isoxazole-diaryl (rigid)"),
 ("M_palmitoyl_amide","aliphatic","amide",2.99,"",0.497,2.031,"med","A2 with ester->amide"),
 ("M_caspo_amide","aliphatic","amide",2.82,"",0.486,1.726,"med","A1 with ester->amide"),
 ("M_anidula_amide","aromatic","amide",3.14,"",0.153,1.962,"med","G_anidula with ester->amide"),
 ("D_caspo_both_ester","dual","di-ester",3.96,"OOD",0.749,3.772,"med","caspo at both positions, di-ester"),
 ("D_caspo_both_amide","dual","di-amide",3.92,"",0.617,3.236,"med","caspo at both positions, di-amide (hardened)"),
 ("ref_ibrexafungerp","reference","-",2.96,"OOD",0.0634,2.025,"med","approved oral GS inhibitor"),
 ("ref_fluconazole","reference","-",0.73,"",1.605,0.650,"high","small oral anchor"),
]

with open("/home/claude/papulacandin_work/rational_design_panel.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["id","class","linkage","LogD_pH7.4","LogD_flag","ADME_permeability","ADME_lipophilicity","solubility_conf","rationale"])
    for r in rows: w.writerow(r)

colors={"native":"#C0392B","aliphatic":"#1F6FB2","aromatic":"#7D3C98","dual":"#E67E22","reference":"#7F8C8D"}
fig,ax=plt.subplots(figsize=(11,7))
pts={}
for r in rows:
    cid,cls,link,logd,flag,perm=r[0],r[1],r[2],r[3],r[4],r[5]
    marker="^" if "amide" in link else ("*" if cls=="reference" else "o")
    ax.scatter(logd,perm,s=200 if cls=="reference" else 150,c=colors[cls],marker=marker,
               edgecolor="black",linewidth=0.8,zorder=3)
    pts[cid]=(logd,perm)
    lab=cid.replace("_"," ").replace("ref ","").replace(" ester","").replace(" amide","")
    ax.annotate(lab,(logd,perm),xytext=(6,6),textcoords="offset points",fontsize=8.5,
                fontweight="bold" if cls!="reference" else "normal")

# ester -> amide connector arrows (matched pairs): shows the effect of hardening the linkage
for est,ami in [("A2_palmitoyl_ester","M_palmitoyl_amide"),
                ("A1_caspo_ester","M_caspo_amide"),
                ("G_anidula_ester","M_anidula_amide")]:
    x0,y0=pts[est]; x1,y1=pts[ami]
    ax.annotate("",xy=(x1,y1),xytext=(x0,y0),
                arrowprops=dict(arrowstyle="->",color="#27AE60",lw=1.4,alpha=0.8),zorder=2)

ax.set_yscale("log")
ax.set_xlabel("Predicted LogD @ pH 7.4  (Inductive Bio, real log units)",fontsize=11)
ax.set_ylabel("Predicted passive permeability  (Boltz ADME-v1, model units, log)",fontsize=11)
ax.set_title("Rational-design property panel: echinocandin grafts + ester→amide bioisosteres",
             fontsize=13,fontweight="bold")
ax.axvspan(2.6,3.1,color="#2ECC71",alpha=0.07,zorder=0)
ax.text(2.85,ax.get_ylim()[0]*1.3,"ibrexafungerp\nLogD band",ha="center",fontsize=8.5,color="#1E8449")

# legend
from matplotlib.lines import Line2D
leg=[Line2D([0],[0],marker='o',color='w',markerfacecolor=colors[k],markeredgecolor='k',markersize=10,label=k)
     for k in ["native","aliphatic","aromatic","dual","reference"]]
leg+=[Line2D([0],[0],marker='o',color='w',markerfacecolor='#BBB',markeredgecolor='k',markersize=10,label='ester (circle)'),
      Line2D([0],[0],marker='^',color='w',markerfacecolor='#BBB',markeredgecolor='k',markersize=10,label='amide (triangle)'),
      Line2D([0],[0],color='#27AE60',lw=1.6,label='ester→amide swap')]
ax.legend(handles=leg,loc="upper left",fontsize=8.5,framealpha=0.9,ncol=2)

note=("Native = permeability floor (relies on tail for delivery).\n"
      "Green arrows: ester→amide hardening RAISES predicted permeability (no penalty).\n"
      "Rigid aromatic grafts (anidula/mica) sit LOWER than flexible aliphatic tails.\n"
      "Orientation only: model-internal scales, out-of-domain scaffold; NO f_u / albumin / potency here.")
ax.text(0.985,0.03,note,transform=ax.transAxes,fontsize=8,va="bottom",ha="right",
        bbox=dict(boxstyle="round,pad=0.4",fc="#FBFCFC",ec="#BDC3C7"))
ax.grid(True,which="both",alpha=0.22); ax.margins(x=0.16)
plt.tight_layout()
plt.savefig("/home/claude/papulacandin_work/rational_design_panel.png",dpi=160)
print("done")

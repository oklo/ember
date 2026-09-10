#!/usr/bin/env python3
"""Plot completed 0.1 Msun tracks with identical interior physics/settings."""
import argparse
import json
from pathlib import Path


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("nongrey",type=Path)
    parser.add_argument("corrected_cond",type=Path)
    parser.add_argument("--output",type=Path,default=Path("docs/results/evolution_nongrey_m010_1tyr"))
    args=parser.parse_args()
    runs=[json.loads(path.read_text()) for path in [args.nongrey,args.corrected_cond]]
    for run in runs:
        if not run["converged"] or run["history"][-1][0]!=1e12 or abs(run["mass_Msun"]-.1)>1e-12:
            raise ValueError("two completed 0.1 Msun trillion-year tracks are required")
    if not runs[0]["atmosphere_model"].startswith("nongrey:") or runs[1]["atmosphere_model"]!="cond-corrected":
        raise ValueError("expected non-grey and corrected-COND tracks, in that order")
    for key in ["points","nuclear_model","transport_model","mass_basis","step_error_tolerances"]:
        if runs[0][key]!=runs[1][key]:
            raise ValueError(f"comparison differs in {key}")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False})
    figure,axes=plt.subplots(2,2,figsize=(9,6.5),sharex=True)
    for run,style,color in zip(runs,["-","--"],["#246e91","#797979"]):
        columns={name:i for i,name in enumerate(run["columns"])}
        def values(name):return [row[columns[name]] for row in run["history"]]
        age=[v/1e12 for v in values("age_yr")]
        for name,label,hue in [("central_X","Hydrogen-1","#b34b24"),
                               ("central_Y3","Helium-3","#246e91")]:
            axes[0,0].plot(age,values(name),style,color=hue,label=label if style=="-" else None)
        he4=[.98-x-y for x,y in zip(values("central_X"),values("central_Y3"))]
        axes[0,0].plot(age,he4,style,color="#657344",label="Helium-4" if style=="-" else None)
        for ax,name,label,scale in [(axes[0,1],"L_Lsun",r"Luminosity ($10^{-3} L_\odot$)",1000),
                                    (axes[1,0],"R_Rsun",r"Radius ($R_\odot$)",1),
                                    (axes[1,1],"Teff_K","Effective temperature (K)",1)]:
            ax.plot(age,[v*scale for v in values(name)],style,color=color)
            ax.set_ylabel(label)
    axes[0,0].set_ylabel("Baryonic mass fraction")
    axes[0,0].legend(frameon=False,fontsize=9,loc="center left",bbox_to_anchor=(.02,.63))
    for ax in axes.flat:
        ax.grid(alpha=.18)
        ax.set_xlim(0,1)
    for ax in axes[1]:ax.set_xlabel("Elapsed time (trillion years)")
    figure.suptitle(r"$0.1\,M_\odot$: atmosphere comparison",y=.985)
    figure.legend([Line2D([],[],color="#246e91"),Line2D([],[],color="#797979",linestyle="--")],
                  ["Non-grey atmosphere","COND + composition correction"],
                  loc="upper center",bbox_to_anchor=(.5,.953),ncol=2,frameon=False)
    figure.text(.5,.012,f"{runs[0]['points']} mass points; identical interior physics and timestep tolerances. "
                "Time starts at the prescribed initial composition.",ha="center",fontsize=8)
    figure.tight_layout(rect=(0,.035,1,.9))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    for suffix in [".png",".pdf"]:figure.savefig(args.output.with_suffix(suffix),dpi=180)
    plt.close(figure)


if __name__=="__main__":
    main()

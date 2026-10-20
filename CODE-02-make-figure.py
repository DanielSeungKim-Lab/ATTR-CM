"""Rebuild the two-panel figure from the deposited data.

Usage: python CODE-02-make-figure.py   ->  FIGURE1.png
Benchmark values and intervals are recomputed here; see CODE-01 for the printed results.
"""
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

SEED, N_BOOT = 20260929, 20000
POOLED = ["ATTR-ACT", "APOLLO-B", "ATTRibute-CM", "CARDIO-TTRansform", "HELIOS-B"]
BLUE, ORANGE, GREY = "#1b3b6f", "#c2571a", "#8a8f98"

race = pd.read_csv("DATA-01-trial-race-counts.csv")
inp = pd.read_csv("DATA-02-benchmark-inputs.csv").set_index("key")["value"]
cty = pd.read_csv("DATA-03-site-counties.csv", dtype={"county_fips": str})

share = lambda b, w: 100.0 * b / (b + w)
def disease_weight(p_pct, ratio):
    p = p_pct / 100.0
    return 100.0 * (ratio * p) / ((1 - p) + ratio * p)

pool = race[race.trial.isin(POOLED)]
obs = share(pool.black.sum(), pool.white.sum())
rng = np.random.default_rng(SEED)
b_arr, w_arr = pool.black.to_numpy(), pool.white.to_numpy()
draws = rng.integers(0, len(pool), size=(N_BOOT, len(pool)))
lo, hi = np.percentile(share(b_arr[draws].sum(1), w_arr[draws].sum(1)), [2.5, 97.5])

r = inp["prevalence_black_per_100000"] / inp["prevalence_white_per_100000"]
nat75 = share(inp["acs_black_75plus"], inp["acs_white_75plus"])
p5 = cty[cty.trial.isin(POOLED)]
cat5 = share(p5.black_alone_75plus.sum(), p5.white_alone_75plus.sum())

# (label, value, general-population benchmark?) ordered as plotted, least to most stringent
rows = [
    ("US adults \u226520 y", share(inp["acs_black_20plus"], inp["acs_white_20plus"]), True),
    ("US adults \u226565 y", share(inp["acs_black_65plus"], inp["acs_white_65plus"]), True),
    ("US adults \u226575 y", nat75, True),
    ("Population in US site counties \u226575 y", cat5, True),
    ("Disease-referenced US adults \u226575 y", disease_weight(nat75, r), False),
    ("Registry, normalized", share(inp["registry_black"], inp["registry_white"]), False),
    ("Diagnosed claims cohort", share(inp["claims_cohort_black_pct"], inp["claims_cohort_white_pct"]), False),
    ("Disease-referenced US site counties \u226575 y", disease_weight(cat5, r), False),
]

per = (cty.groupby("trial")[["black_alone_75plus", "white_alone_75plus"]].sum()
          .join(race.set_index("trial")[["black", "white"]]))
per["county_pct"] = share(per.black_alone_75plus, per.white_alone_75plus)
per["enrolled_pct"] = share(per.black, per.white)
per["ratio"] = per.enrolled_pct / per.county_pct
per = per.sort_values("ratio")

mpl.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
                     "xtick.labelsize": 6, "ytick.labelsize": 6, "figure.dpi": 150,
                     "axes.linewidth": 0.8, "savefig.bbox": "tight"})
fig, (axA, axB) = plt.subplots(2, 1, figsize=(6.9, 6.8),
                               gridspec_kw=dict(height_ratios=[1, 0.9], hspace=0.45))

for k, (lab, val, general) in enumerate(rows):
    y = len(rows) - 1 - k
    col = BLUE if general else ORANGE
    axA.plot([lo / val, hi / val], [y, y], color=col, lw=1.3, solid_capstyle="round")
    axA.plot(obs / val, y, "o", color=col, ms=6.5, zorder=4)
axA.axvline(1, color="0.6", lw=0.9)
axA.set_yticks(range(len(rows)))
axA.set_yticklabels([lab for lab, _, _ in reversed(rows)])
axA.set_xscale("log"); axA.set_xlim(0.14, 1.7)
axA.set_xticks([0.2, 0.3, 0.5, 1.0, 1.5])
axA.get_xaxis().set_minor_locator(mpl.ticker.NullLocator())
axA.get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
axA.set_ylim(-0.75, len(rows) - 0.25)
axA.set_xlabel("Pooled Black proportion / reference proportion")
axA.set_title("Participation ratios against 8 US reference populations", loc="left")
axA.text(1.03, len(rows) - 0.45, "parity", fontsize=6, color="0.45", va="center")

allc = per[per.index.isin(POOLED)]
end = per.loc["ENDEAVOUR"]
y_end = len(allc) + 0.7
for i, (name, row) in enumerate(allc.iterrows()):
    axB.annotate("", xy=(row.enrolled_pct, i), xytext=(row.county_pct, i),
                 arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=1.3, shrinkA=3.5, shrinkB=2.5))
    axB.plot(row.county_pct, i, "o", mfc="white", mec=GREY, mew=1.3, ms=6.5, zorder=3)
    axB.plot(row.enrolled_pct, i, "o", color=BLUE, ms=6.5, zorder=4)
    axB.text(row.enrolled_pct - 0.9, i, f"{row.ratio:.2f}", ha="right", va="center",
             fontsize=6, color=BLUE)
axB.annotate("", xy=(end.enrolled_pct, y_end), xytext=(end.county_pct, y_end),
             arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=1.3, shrinkA=3.5, shrinkB=2.5))
axB.plot(end.county_pct, y_end, "o", mfc="white", mec=GREY, mew=1.3, ms=6.5, zorder=3)
axB.plot(end.enrolled_pct, y_end, "o", color=ORANGE, ms=6.5, zorder=4)
axB.text(end.enrolled_pct + 1.1, y_end, f"{end.ratio:.2f}", ha="left", va="center",
         fontsize=6, color=ORANGE)
axB.axhline(len(allc) - 0.15, color="0.82", lw=0.8, ls=(0, (4, 3)))
axB.set_yticks(list(range(len(allc))) + [y_end])
axB.set_yticklabels(list(allc.index) + ["ENDEAVOUR"])
axB.set_xlim(0, 58); axB.set_ylim(-0.8, y_end + 0.9)
axB.set_xlabel("Black proportion among Black or White individuals, %")
axB.set_title("Global trial enrollment vs population in each trial's US site counties", loc="left")
axB.text(end.county_pct + 1.0, y_end + 0.5,
         "enrollment restricted to pathogenic TTR variant carriers",
         fontsize=6, color=ORANGE, style="italic")
axB.legend(handles=[
    Line2D([], [], marker="o", mfc="white", mec=GREY, mew=1.3, ls="", ms=6.5,
           label="Population aged \u226575 y in that trial's US site counties"),
    Line2D([], [], marker="o", color=BLUE, ls="", ms=6.5,
           label="Global trial enrollment, 5 pooled trials"),
    Line2D([], [], marker="o", color=ORANGE, ls="", ms=6.5,
           label="Global trial enrollment, carrier-restricted trial")],
    loc="lower right", frameon=False, fontsize=6, handletextpad=0.4, borderpad=0.2)
for ax in (axA, axB):
    ax.spines[["top", "right"]].set_visible(False)
for ax, letter in [(axA, "A"), (axB, "B")]:
    ax.text(-0.08, 1.02, letter, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="bottom", ha="left")
fig.savefig("FIGURE1.png", dpi=300)

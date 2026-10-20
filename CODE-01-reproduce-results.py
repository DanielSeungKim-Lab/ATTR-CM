"""Reproduce every number reported in the ATTR-CM representation letter.

Usage: python CODE-01-reproduce-results.py
Requires pandas and numpy. Writes nothing; prints the reported values.
"""
import numpy as np
import pandas as pd

SEED = 20260929
N_BOOT = 20000
POOLED = ["ATTR-ACT", "APOLLO-B", "ATTRibute-CM", "CARDIO-TTRansform", "HELIOS-B"]

race = pd.read_csv("DATA-01-trial-race-counts.csv")
inp = pd.read_csv("DATA-02-benchmark-inputs.csv").set_index("key")["value"]
cty = pd.read_csv("DATA-03-site-counties.csv", dtype={"county_fips": str})

share = lambda b, w: 100.0 * b / (b + w)


def disease_weight(p_pct, ratio):
    """Black share p (percent) rescaled by a Black-to-White prevalence ratio, on the odds scale."""
    p = p_pct / 100.0
    return 100.0 * (ratio * p) / ((1 - p) + ratio * p)


# ---- pooled Black proportion and its cluster bootstrap interval ----
pool = race[race.trial.isin(POOLED)]
black, bw = pool.black.sum(), (pool.black + pool.white).sum()
obs = share(black, pool.white.sum())

rng = np.random.default_rng(SEED)
b_arr, w_arr = pool.black.to_numpy(), pool.white.to_numpy()
draws = rng.integers(0, len(pool), size=(N_BOOT, len(pool)))
reps = share(b_arr[draws].sum(1), w_arr[draws].sum(1))
lo, hi = np.percentile(reps, [2.5, 97.5])

print(f"Pooled: {black} of {bw} Black or White participants were Black "
      f"({obs:.2f}%; 95% CI, {lo:.2f}%-{hi:.2f}%)\n")

# ---- the eight benchmarks ----
r = inp["prevalence_black_per_100000"] / inp["prevalence_white_per_100000"]
nat75 = share(inp["acs_black_75plus"], inp["acs_white_75plus"])

# Counties are summed within a trial, then the five trials' county sets are added.
# A county hosting sites in more than one pooled trial therefore contributes once per trial.
p5 = cty[cty.trial.isin(POOLED)]
cat5 = share(p5.black_alone_75plus.sum(), p5.white_alone_75plus.sum())

bench = {
    "US adults 20 y and over": share(inp["acs_black_20plus"], inp["acs_white_20plus"]),
    "US adults 65 y and over": share(inp["acs_black_65plus"], inp["acs_white_65plus"]),
    "US adults 75 y and over": nat75,
    "US site counties, 75 y and over": cat5,
    "Disease-referenced, US adults 75 y and over": disease_weight(nat75, r),
    "Disease registry, transthyretin stratum": share(inp["registry_black"], inp["registry_white"]),
    "Diagnosed claims cohort": share(inp["claims_cohort_black_pct"], inp["claims_cohort_white_pct"]),
    "Disease-referenced, US site counties 75 y and over": disease_weight(cat5, r),
}
out = pd.DataFrame(
    [(k, v, obs / v, lo / v, hi / v) for k, v in bench.items()],
    columns=["benchmark", "benchmark_pct", "ratio", "ci_lo", "ci_hi"],
).sort_values("benchmark_pct")
print(out.round(3).to_string(index=False))
print(f"\nDiagnosed Black-to-White prevalence ratio: {r:.6f}\n")

# ---- each trial against the counties hosting its own US sites ----
per = (cty.groupby("trial")[["black_alone_75plus", "white_alone_75plus"]].sum()
          .join(cty.groupby("trial").size().rename("counties"))
          .join(race.set_index("trial")[["black", "white"]]))
per["county_pct"] = share(per.black_alone_75plus, per.white_alone_75plus)
per["enrolled_pct"] = share(per.black, per.white)
per["ratio"] = per.enrolled_pct / per.county_pct
print(per[["counties", "county_pct", "enrolled_pct", "ratio"]].round(3)
         .sort_values("ratio").to_string())

# ---- ENDEAVOUR ----
e = race[race.trial == "ENDEAVOUR"].iloc[0]
e_bw = share(e.black, e.white)
e_all = 100.0 * e.black / e.race_reported_total
carrier = (inp["carrier_v142i_share_of_pathogenic_pct"] / 100.0
           * inp["carrier_black_share_of_v142i_pct"])
print(f"\nENDEAVOUR: {e.black} of {e.black + e.white} Black or White participants were Black "
      f"({e_bw:.1f}%), {e_bw / per.loc['ENDEAVOUR', 'county_pct']:.2f} times its county benchmark "
      f"and {e_bw / obs:.1f} times the pooled proportion.")
print(f"Among all {e.race_reported_total} participants, {e.black} were Black ({e_all:.1f}%).")
print(f"Exploratory carrier comparison: {inp['carrier_v142i_share_of_pathogenic_pct']}% x "
      f"{inp['carrier_black_share_of_v142i_pct']}% = {carrier:.1f}%; "
      f"{e_all:.1f}% / {carrier:.1f}% = {e_all / carrier:.2f}")

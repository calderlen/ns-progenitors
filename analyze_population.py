import numpy as np

from posydon.popsyn.synthetic_population import Population


pop = Population("1e+00_Zsun_population.h5")

history = pop.history[["event"]]
oneline = pop.oneline.select()


n_initial = len(oneline)

failed = oneline["FAILED"]
valid = ~failed

cc = history["event"].isin(["CC1", "CC2"])

n_failed = failed.sum()
n_cc_events = cc.sum()
n_cc_systems = history.loc[cc].index.nunique()

n_ccsn = (
    (oneline["S1_SN_type"] == "CCSN").sum()
    + (oneline["S2_SN_type"] == "CCSN").sum()
)

n_ecsn = (
    (oneline["S1_SN_type"] == "ECSN").sum()
    + (oneline["S2_SN_type"] == "ECSN").sum()
)

n_ns = (
    (oneline["S1_state_f"] == "NS").sum()
    + (oneline["S2_state_f"] == "NS").sum()
)

n_bh = (
    (oneline["S1_state_f"] == "BH").sum()
    + (oneline["S2_state_f"] == "BH").sum()
)


s1_ns = oneline["S1_state_f"] == "NS"
s2_ns = oneline["S2_state_f"] == "NS"

s1_ms = oneline["S1_state_f"] == "H-rich_Core_H_burning"
s2_ms = oneline["S2_state_f"] == "H-rich_Core_H_burning"

nsms = (s1_ns & s2_ms) | (s2_ns & s1_ms)

valid_nsms = valid & nsms

detached_nsms = valid_nsms & (oneline["state_f"] == "detached")
rlof_nsms = valid_nsms & (oneline["state_f"] == "initial_RLOF")

gaia_period_nsms = (
    detached_nsms
    & oneline["orbital_period_f"].notna()
    & (oneline["orbital_period_f"] < 1500)
)


print("Initial systems:", n_initial)
print("Failed systems:", n_failed)
print("Core-collapse events:", n_cc_events)
print("Systems reaching core collapse:", n_cc_systems)
print("CCSNe:", n_ccsn)
print("ECSNe:", n_ecsn)
print("NS remnants:", n_ns)
print("BH remnants:", n_bh)

print()
print("Valid NS+MS:", valid_nsms.sum())
print("NS+MS ending at RLOF:", rlof_nsms.sum())
print("Detached NS+MS:", detached_nsms.sum())
print("Detached NS+MS with P < 1500 d:", gaia_period_nsms.sum())

print()
print("Final states of valid NS+MS:")
print(oneline.loc[valid_nsms, "state_f"].value_counts(dropna=False))


survivors = oneline.loc[detached_nsms].copy()

survivors["M_NS"] = np.where(
    survivors["S1_state_f"] == "NS",
    survivors["S1_mass_f"],
    survivors["S2_mass_f"],
)

survivors["M_MS"] = np.where(
    survivors["S1_state_f"] == "NS",
    survivors["S2_mass_f"],
    survivors["S1_mass_f"],
)

cols = [
    "M_NS",
    "M_MS",
    "orbital_period_f",
    "eccentricity_f",
    "state_f",
    "event_f",
]

survivors = survivors[cols]

print()
print("Detached NS+MS summary:")
print(survivors.describe())

survivors.to_csv("data/nsms_survivors.csv")
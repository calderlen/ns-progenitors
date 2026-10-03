NS remnant-mass / Gaia project data package

Core first-pass inputs:
1. gaia2024_analysis_master.csv
2. ertl2020_progenitor_grid.csv
3. ertl2020_remnant_outputs_REQUIRED_TEMPLATE.csv  <-- values still need to be obtained from the Garching Ertl+2020 archive.

For a selection-corrected / publication-grade analysis also use:
- gaia2024_candidate_followup_universe.csv
- gaia2024_selection_cuts.csv
- gaia2024_orbit_fits.csv
- gaia2024_luminous_star_properties.csv

Validation:
- ertl2020_ns_mass_averages_validation.csv
  For the standard W18 setup, reconstruction should recover approximately median M_NS,g = 1.351 Msun and mean = 1.371 Msun.

Historical benchmark:
- pejcha2012_dns_sample_historical.csv

Not included because it is not required unless you re-fit the Gaia orbits from scratch:
- El-Badry+2024 Table 5 raw RV measurements
- Gaia DR3 full astrometric covariance inputs / epoch astrometry

All masses use solar-mass units unless column names specify otherwise.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astroquery.vizier import Vizier

# Download SPARC tables from CDS VizieR
print("Connecting to CDS VizieR to download SPARC catalog data.")
Vizier.ROW_LIMIT = -1
sparc_catalog_id = "J/AJ/152/157"

catalogs = Vizier.get_catalogs(sparc_catalog_id)

# Table 1 = Global properties, Table 2 = Rotation curves
global_table = catalogs[0].to_pandas()
rc_table = catalogs[1].to_pandas()

print("Catalog data successfully downloaded.")

# Filter for Dwarf Irregular Galaxies using SPARC column names
global_table["galaxy_name"] = global_table["Name"].astype(str)

# Exact Table 1 columns: "Type" (Hubble Type), "i" (Inclination), "Vflat" (Flat circular velocity)
type_filter = (global_table["Type"] >= 9) | (global_table["Vflat"] < 80.0)
inc_filter = (global_table["i"] >= 30.0) & (global_table["i"] <= 85.0)
valid_vflat = global_table["Vflat"] > 0

sample = global_table[type_filter & inc_filter & valid_vflat].copy()

# Compute physical masses and gas fractions
# L3.6 is in 10^9 L_sun, MHI is in 10^9 M_sun
sample["m_star"] = sample["L3.6"] * 1e9 * 0.50
sample["m_hi"] = sample["MHI"] * 1e9
sample["m_gas"] = 1.33 * sample["m_hi"]
sample["m_baryon"] = sample["m_star"] + sample["m_gas"]
sample["f_gas"] = sample["m_gas"] / sample["m_baryon"]

print("Extracted dwarf irregular count:", len(sample))

# Extract matching rotation curves from Table 2
rc_table["galaxy_name"] = rc_table["Name"].astype(str)
sample_names = sample["galaxy_name"].values
sample_rc = rc_table[rc_table["galaxy_name"].isin(sample_names)].copy()

# Save clean data tables to CSV
sample.to_csv("sparc_dirr_sample_summary.csv", index=False)
sample_rc.to_csv("sparc_dirr_rotation_curves.csv", index=False)

print("Saved sparc_dirr_sample_summary.csv and sparc_dirr_rotation_curves.csv")

# Plot Baryonic Mass vs Vflat
plt.figure(figsize=(10, 6))
plt.scatter(sample["Vflat"], sample["m_baryon"], color="navy", edgecolors="black", s=60)
plt.yscale("log")
plt.xlabel("Flat Circular Velocity Vflat [km/s]", fontsize=12)
plt.ylabel("Baryonic Mass Mbaryon [M_sun]", fontsize=12)
plt.title("SPARC Dwarf Irregular Sample: Baryonic Mass vs Vflat", fontsize=14, fontweight="bold")
plt.grid(True, linestyle="--", alpha=0.5)

plt.tight_layout()
plt.savefig("phase1_baryonic_tully_fisher.png", dpi=200)
plt.show()

print("Plot saved as phase1_baryonic_tully_fisher.png.")
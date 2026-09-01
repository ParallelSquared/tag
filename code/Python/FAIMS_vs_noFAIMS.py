#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 15:18:04 2025

@author: maddyyeh
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

'''

This figure shows the impact of FAIMS on the coelution of tag byproducts with tagged peptides. 
In the upper left corner, XIC of the tag byproducts are shown with the XIC of tagged peptides without FAIMS. 
In the bottom left corner, XIC of the tag byproducts are shown with the XIC of tagged peptides with FAIMS. 
In the upper right corner, peak area of the tag byproducts are compared across FAIMS and No FAIMS conditions. 
In the bottom right corner, MS1 intensities of intersected precursors (found in FAIMS and No FAIMS conditions) are compared
in a histogram. 

'''

# Define QualBrowser XIC paths
file_groups_FAIMS_free = {
    "FAIMS-free Tag Byproducts": [
        "/Volumes/Lab/MY/2025-03-05_T6_publication_fig/TaginFAIMS_vsnoFAIMS/QualBrowserXICHydrolyzedProd_FAIMS-free_MS1_d8.txt",
        "/Volumes/Lab/MY/2025-03-05_T6_publication_fig/TaginFAIMS_vsnoFAIMS/QualBrowserXICQuenchedProd_FAIMS-free_MS1_d8.txt"
    ],
    "FAIMS-free Labeled Peptides": [
        "/Volumes/Lab/MY/2025-03-05_T6_publication_fig/TaginFAIMS_vsnoFAIMS/QualBrowserXICLabeledPeptides_FAIMS-free_lowMR_d8_500-1000mz_MS1.txt"
    ]
}

file_groups_FAIMS_45 = {
    "FAIMS-45 Tag Byproducts": [
        "/Volumes/Lab/MY/2025-03-05_T6_publication_fig/TaginFAIMS_vsnoFAIMS/QualBrowserXICHydrolyzedProd_FAIMS-45_MS1_d8.txt",
        "/Volumes/Lab/MY/2025-03-05_T6_publication_fig/TaginFAIMS_vsnoFAIMS/QualBrowserXICQuenchedProd_FAIMS-45_MS1_d8.txt"
    ],
    "FAIMS-45 Labeled Peptides": [
        "/Volumes/Lab/MY/2025-03-05_T6_publication_fig/TaginFAIMS_vsnoFAIMS/QualBrowserXICLabeledPeptides_FAIMS-45_lowMR_d8_MS1.txt"
    ]
}

'''
Note for Tag Byproducts: expected masses of the quenched and hydrolyzed tag products were searched for at the MS1 level. 
These traces were exported as .csv and plotted. 

Note for Labeled Peptides: all peaks were extracted at the MS1 level. These traces were exported as .csv and plotted. 
'''

# Set colors for each condition to be preserved throughout the entire figure
colors = {
    "FAIMS-45 Tag Byproducts": "#1f77b4", # BLUE
    "FAIMS-free Tag Byproducts": "#d62728", # RED
    "FAIMS-45 Labeled Peptides": "#262626", # GRAY
    "FAIMS-free Labeled Peptides": "#262626"
}

# Define peak areas (these values were taken from QualBrowser)
peak_areas = {
    "Quenched Tag": {"FAIMS-45 Tag Byproducts": 460409, "FAIMS-free Tag Byproducts": 903671384}, 
    "Hydrolyzed Tag": {"FAIMS-45 Tag Byproducts": 0, "FAIMS-free Tag Byproducts": 91103635},
}

# Load evidence data for MS1 intensity histogram (MaxQuant evidence.txt) and filter by confidence interval and intensity > 0
evidence = pd.read_csv('/Volumes/Lab/MY/2025-03-05_T6_publication_fig/TaginFAIMS_vsnoFAIMS/combined/txt/evidence.txt', sep="\t")
evidence = evidence[evidence['Intensity'] > 0]
evidence = evidence[evidence['PEP'] <= 0.01]

# Create separate data frames for peptides identified in "No FAIMS" (or "FAIMS-free") and "FAIMS" conditions
MS1_intensity_FAIMS_free = evidence[evidence['Raw file'].str.contains("FAIMS-Free-Tag6-d8_5ng_regMR", na=False)]
MS1_intensity_FAIMS = evidence[evidence['Raw file'].str.contains("FAIMS-45-Tag6-d8_5ng_regMR", na=False)]

# Determine which peptides are shared between "No FAIMS" (or "FAIMS-free") and "FAIMS" conditions (SEQUENCE INTERSECTION)
shared_precursors = set(MS1_intensity_FAIMS_free['Modified sequence']) & set(MS1_intensity_FAIMS['Modified sequence'])

# Filter each dataframe to only contain peptides found in common between conditions
shared_df_FAIMS_free = MS1_intensity_FAIMS_free[MS1_intensity_FAIMS_free['Modified sequence'].isin(shared_precursors)]
shared_df_FAIMS = MS1_intensity_FAIMS[MS1_intensity_FAIMS['Modified sequence'].isin(shared_precursors)]


#%%

# Figure set up: use GridSpec to create a figure with multiple subplots
fig = plt.figure(figsize=(10, 8), dpi=300)
gs = GridSpec(2, 2, width_ratios=[1.5, 1], height_ratios=[1, 1], wspace=0.45, hspace=0.45)

################################################################################

# --- Time vs. Intensity Plot (Left side, top) ---

ax1 = fig.add_subplot(gs[0, 0]) # this defines the position of this subplot

# Rename legend elements
legend_label_mapping = { 
    "FAIMS-free Tag Byproducts": "Tag Byproducts",
    "FAIMS-free Labeled Peptides": "Labeled Peptides"
} 

plotted_labels = set()

for label, file_paths in file_groups_FAIMS_free.items(): 
    combined_df = None
    for file_path in file_paths:
        df = pd.read_csv(file_path, sep="\t", skiprows=3)
        df.columns = ["Time", "Intensity"]
        df["Time"] = pd.to_numeric(df["Time"], errors="coerce")
        df["Intensity"] = pd.to_numeric(df["Intensity"], errors="coerce")
        df["Intensity"] = df["Intensity"].where((df["Time"] >= 5) & (df["Time"] <= 40), 0)

        if combined_df is None:
            combined_df = df
        else:
            combined_df["Intensity"] += df["Intensity"]

    legend_label = legend_label_mapping.get(label, label)

    if legend_label not in plotted_labels:
        ax1.plot(combined_df["Time"], combined_df["Intensity"], label=legend_label, color=colors[label], linewidth=3, alpha=0.7)
        plotted_labels.add(legend_label)  # Mark this label as plotted
    else:
        ax1.plot(combined_df["Time"], combined_df["Intensity"], color=colors[label], linewidth=3, alpha=0.7)  # No label

ax1.set_xlabel("Time (min)", fontsize=16, fontname="Arial")
ax1.set_ylabel("Intensity", fontsize=16, fontname="Arial", color="black", labelpad=10)
ax1.set_ylim(-0.15e8,3.5e8)
ax1.tick_params(axis='y', labelcolor="black", labelsize=13)
ax1.tick_params(axis='x', labelcolor="black", labelsize=13)
ax1.set_xlim(5, 40)
#ax1.grid(True, linestyle="--", linewidth=0.5, alpha=0.7)
ax1.set_title("No FAIMS Chromatogram", fontsize=16, fontname="Arial")

# Annotate the leftmost peak as 'Quenched'
ax1.text(18, 2.65e8, "Quenched Tag", 
         fontsize=14, fontname="Arial", 
         ha="right", va="bottom", color="black")

ax1.text(18, 2.35e8, "342.1494 m/z", 
         fontsize=14, fontname="Arial", 
         ha="right", va="bottom", color="black")

# Annotate the rightmost peak as 'Hydrolyzed'
ax1.text(11.5, 0.85e8, "Hydrolyzed Tag", 
         fontsize=14, fontname="Arial", 
         ha="left", va="bottom", color="black")

ax1.text(12, 0.55e8, "327.1347 m/z", 
         fontsize=14, fontname="Arial", 
         ha="left", va="bottom", color="black")

ax1.legend(loc="upper right", fontsize=12, frameon=False)


################################################################################

# --- Time vs. Intensity Plot (Left side, bottom) ---

ax2 = fig.add_subplot(gs[1, 0]) # this defines the position of this subplot

legend_label_mapping = {
    "FAIMS-45 Tag Byproducts": "Tag Byproducts",
    "FAIMS-45 Labeled Peptides": "Labeled Peptides"
}

plotted_labels = set()

for label, file_paths in file_groups_FAIMS_45.items():
    combined_df = None
    for file_path in file_paths:
        df = pd.read_csv(file_path, sep="\t", skiprows=3)
        df.columns = ["Time", "Intensity"]
        df["Time"] = pd.to_numeric(df["Time"], errors="coerce")
        df["Intensity"] = pd.to_numeric(df["Intensity"], errors="coerce")
        df["Intensity"] = df["Intensity"].where((df["Time"] >= 5) & (df["Time"] <= 40), 0)

        if combined_df is None:
            combined_df = df
        else:
            combined_df["Intensity"] += df["Intensity"]
    
    legend_label = legend_label_mapping.get(label, label)

    if legend_label not in plotted_labels:
        if label == "FAIMS-45 Tag Byproducts":
            ax2.plot(combined_df["Time"], combined_df["Intensity"], label=legend_label, color=colors[label], linewidth=3.5, alpha=1, zorder=10)
        else: 
            ax2.plot(combined_df["Time"], combined_df["Intensity"], label=legend_label, color=colors[label], linewidth=2, alpha=0.7, zorder=5)
        plotted_labels.add(legend_label)  # Mark this legend label as plotted
    else:
        ax2.plot(combined_df["Time"], combined_df["Intensity"], color=colors[label], linewidth=2, alpha=0.7, zorder=5)  # No label

ax2.set_xlabel("Time (min)", fontsize=16, fontname="Arial")
ax2.set_ylabel("Intensity", fontsize=16, fontname="Arial", color="black", labelpad=10)
ax2.set_ylim(-0.075e7,1.75e7)
ax2.tick_params(axis='y', labelcolor="black", labelsize=13)
ax2.tick_params(axis='x', labelcolor="black", labelsize=13)
ax2.set_xlim(5, 40)
#ax2.grid(True, linestyle="--", linewidth=0.5, alpha=0.7)
ax2.set_title("FAIMS Chromatogram", fontsize=16, fontname="Arial")
ax2.legend(loc="upper right", fontsize=12, frameon=False)


################################################################################

# --- Bar Plot for Peak Areas (Top right) ---

ax3 = fig.add_subplot(gs[0, 1]) # this defines the position of this subplot

categories = list(peak_areas.keys())  
faims_free_values = [peak_areas[cat]["FAIMS-free Tag Byproducts"] for cat in categories]
faims_45_values = [peak_areas[cat]["FAIMS-45 Tag Byproducts"] for cat in categories]

bar_width = 0.4
x = np.arange(len(categories))

ax3.bar(x + bar_width/2, faims_free_values, bar_width, label="No FAIMS", color=colors["FAIMS-free Tag Byproducts"], alpha=0.7)
ax3.bar(x - bar_width/2, faims_45_values, bar_width, label="FAIMS-45", color=colors["FAIMS-45 Tag Byproducts"], alpha=0.7)
ax3.set_yscale("log")  # <-- Log scale for peak areas
ax3.set_ylim(-1e7,5e10)
ax3.set_xticks(x)
ax3.tick_params(axis='y', labelcolor="black", labelsize=13)
ax3.set_xticklabels(["Quenched\nTag", "Hydrolyzed\nTag"], fontsize=16, fontname="Arial", ha="center")
ax3.set_ylabel("Peak Area, log$_{10}$", fontsize=16, fontname="Arial", labelpad=7)
ax3.legend(fontsize=13, loc="upper right", frameon=False)
ax3.spines["top"].set_visible(False)
ax3.spines["right"].set_visible(False)


################################################################################

# --- Histogram for MS1 Intensities (Bottom right) ---

ax4 = fig.add_subplot(gs[1, 1]) # this defines the position of this subplot

ax4.hist(np.log10(shared_df_FAIMS_free["Intensity"]), alpha=0.7, bins=30, color=colors["FAIMS-free Tag Byproducts"], 
         label="FAIMS-free", orientation='horizontal', density=True)  # Set density=True
ax4.hist(np.log10(shared_df_FAIMS["Intensity"]), alpha=0.7, bins=30, color=colors["FAIMS-45 Tag Byproducts"], 
         label="FAIMS-45", orientation='horizontal', density=True)  # Set density=True
ax4.tick_params(axis='y', labelcolor="black", labelsize=13)
ax4.tick_params(axis='x', labelcolor="black", labelsize=13)
ax4.set_ylabel("MS1 Intensity Intersected\nPrecursors, log$_{10}$", 
               fontsize=16, fontname="Arial", labelpad=10)
ax4.set_xlabel("Density", fontsize=16, fontname="Arial")


################################################################################

plt.show()

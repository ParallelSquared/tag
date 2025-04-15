#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 16:50:13 2025

@author: maddyyeh
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap

'''

This figure shows a heatmap of tagged peptide ions with m/z on the x-axis and mobility on the y-axis. 
The data used to generate this figures is 200 pg tagged peptides on a timsTOF SCP. 
Black boxes are added to show which groups of peptides were selected for analysis (these are called diaPASEF windows).
Purple boxes are added to show ion mobility and m/z of tag byproducts relative to tagged peptides on the heat map. 
Tag byproduct ions are not selected for analysis due to their placement outside of the black windows. 

'''

# Set global font to Arial
plt.rcParams.update({'font.family': 'Arial'})

# Load heatmap data
heatmap_csv = pd.read_csv("/Volumes/Lab/MY/2025-03-28_timsTOF_T6_HeatMap/DIA_heatmap_200pg.tsv", sep="\t")
x = heatmap_csv["Mass"]
y = heatmap_csv["Mobility"]
z = heatmap_csv["Intensity"]
heatmap_df = heatmap_csv.pivot(index="Mobility", columns="Mass", values="Intensity")
heatmap_df = np.log1p(heatmap_df) # Apply log transformation to intensity values


# Define colormap
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", ["navy", "white", "red"])

# Create figure
fig, ax = plt.subplots(figsize=(10, 5))
c = ax.pcolormesh(heatmap_df.columns, heatmap_df.index, heatmap_df.values, cmap=custom_cmap, shading='auto')

# Add colorbar
cbar = plt.colorbar(c, ax=ax)
cbar.set_label('Intensity', fontsize=18)  # Increased fontsize for colorbar label
cbar.ax.tick_params(labelsize=14)  # Adjust label size for colorbar ticks

# Load window data exported from timsTOF Method Editor
window_df = pd.read_csv("/Volumes/Lab/MY/2025-03-28_timsTOF_T6_HeatMap/diaParameters.txt", sep=",")
window_df = window_df.dropna(subset=["Start IM [1/K0]", "End IM [1/K0]", "Start Mass [m/z]", "End Mass [m/z]"])

# Convert columns to numeric
window_df["Start IM [1/K0]"] = pd.to_numeric(window_df["Start IM [1/K0]"], errors="coerce")
window_df["End IM [1/K0]"] = pd.to_numeric(window_df["End IM [1/K0]"], errors="coerce")
window_df["Start Mass [m/z]"] = pd.to_numeric(window_df["Start Mass [m/z]"], errors="coerce")
window_df["End Mass [m/z]"] = pd.to_numeric(window_df["End Mass [m/z]"], errors="coerce")

# Filter PASEF windows
pasef_windows = window_df[window_df["#MS Type"] == "PASEF"]

# Overlay windows as rectangles
for _, row in pasef_windows.iterrows():
    start_im, end_im = row["Start IM [1/K0]"], row["End IM [1/K0]"]
    start_mass, end_mass = row["Start Mass [m/z]"], row["End Mass [m/z]"]

    # Create rectangle (x=start_mass, y=start_im, width, height)
    rect = patches.Rectangle(
        (start_mass, start_im),  # Bottom-left corner
        end_mass - start_mass,   # Width
        end_im - start_im,       # Height
        linewidth=1.5,
        edgecolor="black",     
        facecolor="none"         # Transparent inside
    )
    ax.add_patch(rect)

####################################

# Define the custom rectangle coordinates
start_mass_custom = 324
end_mass_custom = 360
start_im_custom = 0.83
end_im_custom = 0.94

# Create the rectangle for the custom window
rect_custom = patches.Rectangle(
    (start_mass_custom, start_im_custom),  # Bottom-left corner
    end_mass_custom - start_mass_custom,    # Width
    end_im_custom - start_im_custom,       # Height
    linewidth=2,
    edgecolor="#c701ff",                      
    facecolor="none"                       # Transparent inside
)

# Add the custom rectangle to the plot
ax.add_patch(rect_custom)

####################################

# Define the custom rectangle coordinates
start_mass_custom = 420
end_mass_custom = 433
start_im_custom = 0.92
end_im_custom = 0.98

# Create the rectangle for the custom window
rect_custom = patches.Rectangle(
    (start_mass_custom, start_im_custom),  # Bottom-left corner
    end_mass_custom - start_mass_custom,    # Width
    end_im_custom - start_im_custom,       # Height
    linewidth=2,
    edgecolor="#c701ff",                   
    facecolor="none"                       # Transparent inside
)

# Add the custom rectangle to the plot
ax.add_patch(rect_custom)

####################################

# Create a legend for the rectangles
legend_handles = [
    patches.Rectangle((0, 0), 1, 1, linewidth=1.5, edgecolor="black", facecolor="none", label="Analyzed peptide ions"),
    patches.Rectangle((0, 0), 1, 1, linewidth=1.5, edgecolor="#c701ff", facecolor="none", label="Tag byproducts")
]

# Legend
ax.legend(handles=legend_handles, loc='upper right', fontsize=16)

####################################

# Labels and title
ax.set_xlabel("Mass [m/z]", fontsize=18)  
ax.tick_params(axis='x', labelsize=16)  
ax.set_xlim(250, 1200) 
ax.set_ylabel("Mobility [1/K0]", fontsize=18)  
ax.tick_params(axis='y', labelsize=16) 


plt.show()

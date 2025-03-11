#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  4 15:23:48 2025

@author: maddyyeh
"""

import pandas as pd
import matplotlib.pyplot as plt

file_paths = [
    "/Volumes/Lab/MY/2025-03-05_T6_publication_fig/QualBrowserXICQuenchedProd.txt",
    "/Volumes/Lab/MY/2025-03-05_T6_publication_fig/QualBrowserXICHydrolyzedProd.txt", 
    "/Volumes/Lab/MY/2025-03-05_T6_publication_fig/QualBrowserXICLabeledPeptides.txt",
]

labels = ["Quenched Tag", "Hydrolyzed Tag", "Tagged Peptides"]
colors = ["#1f77b4", "#d62728", "#2ca02c"] 

plt.figure(figsize=(8, 6), dpi=300)  

for file_path, label, color in zip(file_paths, labels, colors):
    df = pd.read_csv(file_path, sep="\t", skiprows=3) # this skips metadata in first three lines
    df.columns = ["Time", "Intensity"]  
        
    df["Time"] = pd.to_numeric(df["Time"], errors="coerce")
    df["Intensity"] = pd.to_numeric(df["Intensity"], errors="coerce")
    df = df[(df["Time"] >= 10) & (df["Time"] <= 50)]

    plt.plot(df["Time"], df["Intensity"], label=label, color=color, linewidth=2, alpha=0.85)

plt.xlabel("Time (min)", fontsize=16, fontweight="bold", fontname="serif")
plt.ylabel("Intensity", fontsize=16, fontweight="bold", fontname="serif")
plt.xticks(fontsize=14, fontname="serif")
plt.yticks(fontsize=14, fontname="serif")

plt.xlim(10, 50)
plt.grid(True, linestyle="--", linewidth=0.5, alpha=0.7)

plt.title("Hydrolyzed and Quenched Tag Coelutes with Peptides", fontsize=18, fontweight="bold", fontname="serif")
plt.legend(fontsize=16, frameon=False)

plt.show()

# Qualbrowser data used is 2024-10-03_MY-KG-AS_T6-Clean-Up_Opt-2_C1-C18-SCX_100ng for labeled peptides
# and 2024-10-04_MY-KG-AS_T6_Clean-Up_Opt-2_C3-C18-new-SCX_100ng_lowMR


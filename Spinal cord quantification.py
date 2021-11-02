#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 11:59:15 2021

@author: Chia-Teng
"""

import pandas as pd
import matplotlib.pyplot as plt
import os

def read_Fl_files(file_folder, file_group, i):
    """
    Input: (where csv files are?)
        file_folder example: /Volumes/Teng/FL MS/Analysis/
        file_group example: 0dpt/nes_EW/
        i: sample #, int or str
    The folder should contains cvs files named #R/#G/#RB/#GB.csv,
    B stands for background.
    return a tuple of pandas DataFrame:
    (counts, R_BG_means, R_means_subtract_Bg, R_integrate,
             G_BG_means, G_means_subtract_Bg, G_integrate)
    Columns: X, Value 
    """
    # Read backgound csv
    R_BG = pd.read_csv(file_folder + file_group + str(i+1) + "RB.csv")
    G_BG = pd.read_csv(file_folder + file_group + str(i+1) + "GB.csv")
    
    # Average background intensity according to their X position
    R_BG_means = R_BG.groupby("X").mean().drop(columns = "Y").loc[0:2751]
    G_BG_means = G_BG.groupby("X").mean().drop(columns = "Y").loc[0:2751]
    
    # Read and process the mCherry intensity file
    Huc_data = pd.read_csv(file_folder + file_group + str(i+1) + "R.csv") # Read file
    counts = Huc_data.groupby("X").count().drop(columns = "Y").loc[0:2751] # Count pixel # (diameter) across X position
    R_means = Huc_data.groupby("X").mean().drop(columns = "Y").loc[0:2751] # Average mCherry intensity according to their X position
    R_means_subtract_Bg = R_means - R_BG_means # Subtract mean mCherry intensity to their corresponding backgound at the same X position
    R_integrate = R_means_subtract_Bg * counts # Mutiply background subtracted-mean mCherry intensity by their diameter (in pixel) at the same X position
    
    # Read and process the GFP intensity file, same process as mCherry
    GFP_data = pd.read_csv(file_folder + file_group + str(i+1) + "G.csv")
    G_means = GFP_data.groupby("X").mean().drop(columns = "Y").loc[0:2751]
    G_means_subtract_Bg = G_means - G_BG_means
    G_integrate = G_means_subtract_Bg * counts
    
    return (counts,
            R_BG_means, R_means_subtract_Bg, R_integrate,
            G_BG_means, G_means_subtract_Bg, G_integrate)

def ROI_means(means, X_range):
    """Calculate mean intensity of given ROI range (x1, x2)"""
    min_X, max_X = X_range
    sub_data = means.iloc[min_X: max_X]
    return sub_data.mean()

def ROI_diameter(counts, X_range):
    """mean diameter in ROI in μm"""
    min_X, max_X = X_range
    sub_data = counts.iloc[min_X:max_X]
    return sub_data.mean()*.908

def plot_means_vs_background(ID, title, 
                             means_subtract_Bg, BG_means, 
                             Y_min, Y_max, 
                             color = "#DC5C60", 
                             Fluorescence = "Fluorescence"):
    """
    Input:
        ID: assign a Figure ID in str. You can name by yourself. e.g. (dpf+treatment+genotype+str(i+1))
        Title: enter a str as the figure title. e.g. (dpt + "-" + treatment + "-" + Fluorescence)
        means_subtract_Bg, BG_means: pandas Dataframes
        Y_min, Y_max: Range of Y axis. RFP -2000, 1600 or GFP -7000, 7000
        color: color of Fluorescence. "#DC5C60" for red and "#70BE71" for green
    Output:
        Plot with Fluorescence mean and background (negative) intensity
    """
    plt.figure(ID)
    plt.fill_between(means_subtract_Bg.index * .908, 0, means_subtract_Bg.Value, color = color) # Fill real fluorescent signals above
    plt.fill_between(BG_means.index * .908, -BG_means.Value, color = "#D9D9D9") # Fill background signals, colored in gray below
    plt.xlim(0, 2500) # X axis range in μm
    plt.ylim(Y_min, Y_max)
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    plt.xlabel("position (μm)")
    plt.ylabel(Fluorescence + " mean intensity (AU)")
    plt.title(title)
    # plt.savefig(file_folder + "Figures/" + Fluorescence + "_mean" + dpt + "-" + treatment + "-" + genotype + Fluorescence + str(i+1) + ".pdf")

def conditionally_add_integrated_plot(ID, integrate, treatment, treatments):
    """
    Input:
        ID: assign a Figure ID in str. You can name by yourself. e.g. dpt + "HuC"
        integrate: A pandas Dataframe
    """
    plt.figure(ID)
    if treatment == treatments[0]: # Draw control
        plt.plot(integrate.index * .908, integrate.Value, '-', ms = 5, color= "#1397F1", alpha = 0.2, label = "control")
    elif treatment == treatments[1]: # Draw treatment
        plt.plot(integrate.index * .908, integrate.Value, '-', ms = 5, color= "#EE2B2A", alpha = 0.2, label = "treatment")
    else:
        raise ValueError("You may have more than one treatment group. Please optimize the code.")

def set_integrated_plot(ID, title, fluorescence, Y_max):
    plt.figure(ID)
    plt.xlabel("position (μm)")
    plt.ylabel(fluorescence + " integrated intensity (AU)")
    plt.title(title)
    plt.xlim(0, 2500)
    plt.ylim(0, Y_max)
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())

# This program only takes file folders arraged in a specific order and naming
# e.g. /Volumes/Teng/FL MS/Analysis/0dpt/nes_EW/

file_folder = "/Volumes/Teng/FL MS/Analysis/" # Change this to your folder directory
max_dpt = 4 # What is the max dpt you have?
treatments = ["EW", "MTZ"] # Enter your treatment as list. 
genotypes = ["nes", "sox2"]
max_sample = 11 # what is the largest sample size you have? Be free to enter a large number.
X_start, X_end = 500, 1500 # start/end X position, in μm
Plot_or_not = "Yes" # "Yes" to generate plots; otherwise "No"

# Create directory for plots
try:
    os.makedirs(file_folder + "/Figures")
except OSError:
    print ("Unable to create the directory for plots or it's already existed")
else:
    print ("Successfully created directories for plots")

# Create an empty pandas DataFrame
columns = ["genotype",
           "dpt",
           "treatment", 
           "mCherry", "GFP", 
           "SC_diameter", 
           "Weighted_mCherry",
           "Weighted_GFP"]
data = pd.DataFrame(columns = columns)

# Loop over the sample files
for dpt in [str(i) + "dpt" for i in range(max_dpt + 1)]:
    for treatment in treatments:
        for genotype in genotypes:
            file_group  = dpt + "/" + genotype + "_" + treatment + "/" # like "4dpt/sox2_MTZ/"
            for i in range(max_sample):

                try: # Loop until no more samples to analyze
                    
                    # Retrieve processed data
                    (counts,
                     R_BG_means, R_means_subtract_Bg, R_integrate,
                     G_BG_means, G_means_subtract_Bg, G_integrate) = read_Fl_files(file_folder, file_group, i)
                    
                    # Calculation based on ROI
                    X_range = (int((X_start+0.5)/.908), int((X_end+0.5)/.908)) # converted to pixel and int
                    mCherry_intensity = ROI_means(R_means_subtract_Bg, X_range)
                    GFP_intensity = ROI_means(G_means_subtract_Bg, X_range)
                    SC_diameter = ROI_diameter(counts, X_range)
                    
                    # append process/calculated data
                    sub_data = pd.DataFrame(columns = columns)
                    sub_data.loc[i] = [genotype,
                                       dpt,
                                       treatment,
                                       mCherry_intensity.Value,
                                       GFP_intensity.Value,
                                       SC_diameter.Value,
                                       (mCherry_intensity * SC_diameter).Value,
                                       (GFP_intensity * SC_diameter).Value]
                    data = data.append(sub_data)
                    
                    if Plot_or_not == "Yes":
                        
                        # plot mCherry mean intensity-background intensity figure, be free to commend these out if not needing it
                        plot_means_vs_background(dpt + treatment + str(i+1),
                                                 dpt + "_" + treatment + "_HuC",
                                                 R_means_subtract_Bg, R_BG_means,
                                                 -2000, 16000,
                                                 color = "#DC5C60",
                                                 Fluorescence = "mCherry")
                        plt.savefig(file_folder + "Figures/HuC_" + dpt + "_" + treatment + "_" + genotype + "_" + str(i+1) + ".pdf")
    
                        # plot GFP mean intensity-background intensity figure, be free to commend these out if not needing it
                        plot_means_vs_background(dpt + treatment + genotype + str(i+1),
                                                 dpt + "_" + treatment + "_" + genotype,
                                                 G_means_subtract_Bg, G_BG_means,
                                                 -7000, 7000,
                                                 color = "#70BE71",
                                                 Fluorescence = "GFP")
                        plt.savefig(file_folder + "Figures/" + genotype + "_GFP_" + dpt + "_" + treatment + "_" + str(i+1) + ".pdf")
    
                        # Plot and merge mCherry integrated intensity
                        conditionally_add_integrated_plot(dpt + "HuC", R_integrate, treatment, treatments)
                        # Plot and merge GFP integrated intensity
                        conditionally_add_integrated_plot(dpt + genotype, G_integrate, treatment, treatments)
                        # Plot and merge SC diameters
                        conditionally_add_integrated_plot(dpt + "SC diameter", counts*.908, treatment, treatments)
                    
                except FileNotFoundError: # when running out of samples in a group
                    
                    if Plot_or_not == "Yes":
                        # set parameters of mCherry integrated plots
                        set_integrated_plot(dpt + "HuC", dpt + "_" + "HuC", "mCherry", 2000000)
                        plt.savefig(file_folder + "Figures/" + "integrated_HuC_" + dpt + ".pdf")                 
                        # merge and add legends on GFP integrated plots       
                        set_integrated_plot(dpt + genotype, dpt + "_" + genotype, "GFP", 400000)
                        plt.savefig(file_folder + "Figures/" + "integrated_" + genotype + "_GFP_" + dpt + ".pdf")
                        # set parameters of mCherry integrated plots
                        set_integrated_plot(dpt + "SC diameter", dpt + " spinal cord diameter", "", 100)
                        plt.ylabel("Spinal cord diameter (μm)")
                        plt.ticklabel_format(axis="y", style = "plain", scilimits=(0,0))
                        plt.savefig(file_folder + "Figures/" + "SCdiameter_" + dpt + ".pdf")     
                    
                    print(str(i), "samples in", genotype, dpt, treatment, "have been completed.")
                    break


data.to_csv(file_folder + "data.csv", index = False)
# Bass Connections Fall 2022
# Author: Nikhil Gadiraju

# %% NOTES
# The purpose of this file is to read the output of the R statistical analysis file (provides corrected p-values for 3 treatment comparison groups)
# The p-values are then used to identify the top relevant regions for a given treatment comparison group and plotted via violin plots

# %% PACKAGE/MODULE IMPORTS
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# %% FUNCTION DEFINITIONS
# Define function to rename brain regions from the index.csv names
def structure_update(structure, hemisphere):
    struc_sep = " ".join(structure.split("_"))
    if hemisphere.lower() == 'left':
        new_struc = struc_sep + " (L)"
    elif hemisphere.lower() == 'right':
        new_struc = struc_sep + " (R)"
    elif hemisphere.lower() == 'bilateral':
        new_struc = struc_sep + " (L,R)"
    return new_struc

# Rename index.csv abbreviations by indicating brain location (right or left)
def abbreviation_update(abbrev, hemisphere):
    if hemisphere.lower() == 'left':
        new_abbrev = abbrev + "_L"
    elif hemisphere.lower() == 'right':
        new_abbrev = abbrev + "_R"
    else:
        new_abbrev = abbrev
    return new_abbrev


# %% DATA IMPORTS
# FILE: posthoc_standard.csv
df = pd.read_csv("/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22/Statistical Analysis/posthoc_standard.csv")
# FILE: index.csv
index = pd.read_csv("/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22/Reference Files/Absolute Files/index.csv")
# FILE: voxelvolumes.csv
df_vol = pd.read_csv("/Users/nikhilgadiraju/Box Sync/Home Folder nvg6/Sharing/Bass Connections/Processed Data/Mean Intensity & Voxel Volumes/voxelvolumes.csv")

# %% DATA PROCESSING
# Utilize the defined functions to create updated structures and abbreviation lists to eventually store through a dictionary
index_df = index.iloc[:,0:3]
structure_updated = list(map(lambda x,y: structure_update(x,y), index_df["Structure"], index_df["Hemisphere"]))
abbreviation_updated = list(map(lambda x,y: abbreviation_update(x,y), index_df["Abbreviation"], index_df["Hemisphere"]))
name_dict = {abbreviation_updated[i]: structure_updated[i] for i in range(len(abbreviation_updated))}

# Sort through posthoc_standard.csv to identify relevant regions for each comparison group based on corrected p-value
top_regs = []
top_regs_num = 5
comp_groups = ['ST','SW','TW']
for comparison in comp_groups:
    df_comp = df.sort_values(by=[comparison], ascending=True)
    top = pd.DataFrame({'Abbreviation': df_comp.iloc[:,0], '{} p-value'.format(comparison): df_comp.loc[:,comparison]}).head(top_regs_num)
    top_regs.append(top) # Append relevant regions based on comparison group to empty top_regs list

# Replace naming convention for different treatment groups to more accurately describe exercise conditions
df_vol = df_vol.replace({"Treatment": "wheel_only"}, "voluntary")
df_vol = df_vol.replace({"Treatment": "treadmill"}, "voluntary\n+ enforced")

# %% VIOLIN PLOTS
sns.set_theme()
sns.set(font_scale = 1.2)

title_groups = ['Sedentary vs. Voluntary + Enforced Exercise','Sedentary vs. Voluntary Exercise','Voluntary vs. Voluntary + Forced Exercise']
title_dict = {comp_groups[i]: title_groups[i] for i in range(len(comp_groups))} # Define dictionary to associate each comparison condition to the full comparison title

plot_num = 3
for listframe in range(len(comp_groups)): # This loop iterates through each comparison condition (within the comp_groups list); results in 3 output figures
    fig, axes = plt.subplots(1, plot_num, figsize=(15, 5))
    fig.suptitle('Most Significant Regions via Post-Hoc Analysis\n({})'.format(title_dict[comp_groups[listframe]]),
                 y=1.15, fontsize=16)
    fig.text(0.04, 0.5, 'Voxel Count', va='center', rotation='vertical')
    for val in range(plot_num): # This list iterates through each plot within a given figure (for the top significant regions as specified by plot_num)
        name = top_regs[listframe].iloc[val,0]
        reg = name_dict[name]
        sns.violinplot(ax = axes[val], x = "Treatment", y = name, data = df_vol.filter(['Treatment', name], axis=1)).set_ylabel("")
        axes[val].set_title(reg, fontsize = 14)
        axes[val].set_xlabel("")
    # Save figures in specified folder
    plt.savefig('{}/{}.png'.format('/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22/Statistical Analysis/Output Figures',
                                   comp_groups[listframe]))





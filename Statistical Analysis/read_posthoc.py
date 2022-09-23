# VoxVol_Treatment.py
# Author: Nikhil Gadiraju

# %% PACKAGE/MODULE IMPORTS
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# %% FUNCTION DEFINITIONS
def structure_update(structure, hemisphere):
    struc_sep = " ".join(structure.split("_"))
    if hemisphere.lower() == 'left':
        new_struc = struc_sep + " (L)"
    elif hemisphere.lower() == 'right':
        new_struc = struc_sep + " (R)"
    elif hemisphere.lower() == 'bilateral':
        new_struc = struc_sep + " (L,R)"
    return new_struc

def abbreviation_update(abbrev, hemisphere):
    if hemisphere.lower() == 'left':
        new_abbrev = abbrev + "_L"
    elif hemisphere.lower() == 'right':
        new_abbrev = abbrev + "_R"
    else:
        new_abbrev = abbrev
    return new_abbrev


# %% DATA IMPORTS
df = pd.read_csv("/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22/Statistical Analysis/posthoc_standard.csv")
index = pd.read_csv("/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22/Reference Files/Absolute Files/index.csv")
df_vol = pd.read_csv("/Users/nikhilgadiraju/Box Sync/Home Folder nvg6/Sharing/Bass Connections/Processed Data/Mean Intensity & Voxel Volumes/voxelvolumes.csv")

# %% DATA PROCESSING
index_df = index.iloc[:,0:3]
structure_updated = list(map(lambda x,y: structure_update(x,y), index_df["Structure"], index_df["Hemisphere"]))
abbreviation_updated = list(map(lambda x,y: abbreviation_update(x,y), index_df["Abbreviation"], index_df["Hemisphere"]))
name_dict = {abbreviation_updated[i]: structure_updated[i] for i in range(len(abbreviation_updated))}

top_regs = []
top_regs_num = 3
comp_groups = ['ST','SW','TW']
for comparison in comp_groups:
    df_comp = df.sort_values(by=[comparison], ascending=True)
    top = pd.DataFrame({'Abbreviation': df_comp.iloc[:,0], '{} p-value'.format(comparison): df_comp.loc[:,comparison]}).head(top_regs_num)
    top_regs.append(top)

df_vol = df_vol.replace({"Treatment": "wheel_only"}, "voluntary")
df_vol = df_vol.replace({"Treatment": "treadmill"}, "voluntary\n+ enforced")

# %% VIOLIN PLOTS
sns.set_theme()
sns.set(font_scale = 1.2)

title_groups = ['Sedentary vs. Voluntary + Enforced Exercise','Sedentary vs. Voluntary Exercise','Voluntary vs. Voluntary + Forced Exercise']
title_dict = {comp_groups[i]: title_groups[i] for i in range(len(comp_groups))}

# TODO: Update voxvol CSV so that reported figures have a Y-axis of proportinoal brain volumes rather than voxel count
# TODO: Consider how to repress matplotlib deprecation errors
for listframe in range(len(top_regs)):
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle('Most Significant Regions via Post-Hoc Analysis\n({})'.format(title_dict[comp_groups[listframe]]),
                 y=1.15, fontsize=16)
    fig.text(0.04, 0.5, 'Voxel Count', va='center', rotation='vertical')
    for val in range(top_regs_num):
        name = top_regs[listframe].iloc[val,0]
        reg = name_dict[name]
        sns.violinplot(ax = axes[val], x = "Treatment", y = name, data = df_vol.filter(['Treatment', name], axis=1)).set_ylabel("")
        axes[val].set_title(reg, fontsize = 14)
        axes[val].set_xlabel("")
    plt.savefig('{}/{}.png'.format('/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22/Statistical Analysis/Output Figures',
                                   comp_groups[listframe]))





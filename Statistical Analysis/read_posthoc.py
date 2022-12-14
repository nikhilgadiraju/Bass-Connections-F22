# Bass Connections Fall 2022
# Author: Nikhil Gadiraju

# %% NOTES
# The purpose of this file is to read the output of the R statistical analysis file (provides corrected p-values for 3 treatment comparison groups)
# The p-values are then used to identify the top relevant regions for a given treatment comparison group and plotted via violin plots

# %% PACKAGE/MODULE IMPORTS
import pandas as pd
import os

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
# FILE: posthoc_comparison_stats.csv
df = pd.read_csv("/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22/Reference Files/User-generated Files/posthoc_comparison_stats.csv")
# FILE: index.csv
index = pd.read_csv("/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22/Reference Files/Absolute Files/index.csv")
# FILE: voxelvolumes.csv
df_vol = pd.read_csv(
    "/Users/nikhilgadiraju/Box Sync/Home Folder nvg6/Sharing/Bass Connections/Processed Data/Mean Intensity & Voxel Volumes/voxelvolumes.csv")

# %% DATA PROCESSING
# Utilize the defined functions to create updated structures and abbreviation lists to eventually store through a dictionary
df = df.rename(columns={"ST Comparison Group Effect Size": 'ST',
               'SW Comparison Group Effect Size': 'SW', 'TW Comparison Group Effect Size': 'TW'})
df = df.rename(columns={"ST Comparison Group Pvalue": 'ST_p',
               'SW Comparison Group Pvalue': 'SW_p', 'TW Comparison Group Pvalue': 'TW_p'})
df = df.rename(columns={"ST Comparison Group Lower CI": 'ST_lci',
               'SW Comparison Group Lower CI': 'SW_lci', 'TW Comparison Group Lower CI': 'TW_lci'})
df = df.rename(columns={"ST Comparison Group Higher CI": 'ST_hci',
               'SW Comparison Group Higher CI': 'SW_hci', 'TW Comparison Group Higher CI': 'TW_hci'})

index_df = index.iloc[:, [0, 1, 2, 8]]
structure_updated = list(map(lambda x, y: structure_update(
    x, y), index_df["Structure"], index_df["Hemisphere"]))
abbreviation_updated = list(map(lambda x, y: abbreviation_update(
    x, y), index_df["Abbreviation"], index_df["Hemisphere"]))
name_dict = {abbreviation_updated[i]: structure_updated[i]
             for i in range(len(abbreviation_updated))}
index_dict = list(index_df["index"])
mean_dict = {"S": "sedentary", "W": "voluntary", "T": "voluntary + enforced"}

# Output Structure-Abbreviation CSV
os.chdir("/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22/Reference Files/User-generated Files/Internally Referenced")
pd.DataFrame(data={
    'Structure': structure_updated,
    'Abbreviation': abbreviation_updated,
    'Index': index_dict
}).to_csv("struc_abbrev.csv", encoding='utf-8', mode='w', index=False)

# Sort through posthoc_standard.csv to identify relevant regions for each comparison group based on corrected p-value
top_regs_pos = []
top_regs_neg = []

comp_groups = ['ST', 'SW', 'TW']
for comparison in comp_groups:
    df_iso = df.loc[:, ["Unnamed: 0", "{}_p".format(comparison), comparison,
                        "{}_lci".format(comparison),
                        "{}_hci".format(comparison),
                        "Mean sedentary group", "Mean voluntary group",
                        "Mean voluntary + enforced group",
                        "SD sedentary group", "SD voluntary group",
                        "SD voluntary + enforced group"]]
    df_sig = df_iso[df_iso["{}_p".format(comparison)] <= 0.05]

    df_pos = df_sig.sort_values(by=[comparison], ascending=False).query(
        "{} >= 0".format(comparison))
    top_pos = pd.DataFrame(
        {
            'Abbreviation': df_pos.iloc[:, 0],
            'Structure': [name_dict[key] for key in df_pos.iloc[:, 0]],
            'Mean Volume ({})'.format(mean_dict[comparison[0]].capitalize()):
            df_pos.loc[:, "Mean {} group".format(mean_dict[comparison[0]])],

            'Mean Volume ({})'.format(mean_dict[comparison[1]].capitalize()):
            df_pos.loc[:, "Mean {} group".format(mean_dict[comparison[1]])],

            'SD ({})'.format(mean_dict[comparison[0]].capitalize()):
            df_pos.loc[:, "SD {} group".format(mean_dict[comparison[0]])],

            'SD ({})'.format(mean_dict[comparison[1]].capitalize()):
            df_pos.loc[:, "SD {} group".format(mean_dict[comparison[1]])],

            'P-value': df_pos.loc[:, "{}_p".format(comparison)],
            'Effect Size': df_pos.loc[:, comparison],
            'Lower Confidence Interval':
                df_pos.loc[:, "{}_lci".format(comparison)],
            'Higher Confidence Interval':
                df_pos.loc[:, "{}_hci".format(comparison)]
        })
    # Append relevant regions based on comparison group to empty top_regs list
    top_regs_pos.append(top_pos)

    df_neg = df_sig.sort_values(by=[comparison], ascending=True).query(
        "{} < 0".format(comparison))
    top_neg = pd.DataFrame(
        {
            'Abbreviation': df_neg.iloc[:, 0],
            'Structure': [name_dict[key] for key in df_neg.iloc[:, 0]],
            'Mean Volume ({})'.format(mean_dict[comparison[0]].capitalize()):
            df_neg.loc[:, "Mean {} group".format(mean_dict[comparison[0]])],

            'Mean Volume ({})'.format(mean_dict[comparison[1]].capitalize()):
            df_neg.loc[:, "Mean {} group".format(mean_dict[comparison[1]])],

            'SD ({})'.format(mean_dict[comparison[0]].capitalize()):
            df_neg.loc[:, "SD {} group".format(mean_dict[comparison[0]])],

            'SD ({})'.format(mean_dict[comparison[1]].capitalize()):
            df_neg.loc[:, "SD {} group".format(mean_dict[comparison[1]])],

            'P-value': df_neg.loc[:, "{}_p".format(comparison)],
            'Effect Size': df_neg.loc[:, comparison],
            'Lower Confidence Interval':
                df_neg.loc[:, "{}_lci".format(comparison)],
            "Higher Confidence Interval":
                df_neg.loc[:, "{}_hci".format(comparison)]
        })
    top_regs_neg.append(top_neg)

# Replace naming convention for different treatment groups to more accurately describe exercise conditions
df_vol = df_vol.replace({"Treatment": "wheel_only"}, "voluntary")
df_vol = df_vol.replace({"Treatment": "treadmill"}, "voluntary + enforced")

# %% DATAFRAME OUTPUT
os.chdir("/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22/Reference Files/User-generated Files")
for ind in range(len(comp_groups)):
    top_regs_pos[ind].to_csv('{0}_regs/top_pos_regions_{0}.csv'.format(
        comp_groups[ind].lower()), encoding='utf-8', mode='w', index=False)
    top_regs_neg[ind].to_csv('{0}_regs/top_neg_regions_{0}.csv'.format(
        comp_groups[ind].lower()), encoding='utf-8', mode='w', index=False)

# Bass Connections Fall 2022
# Author: Nikhil Gadiraju

# %% NOTES
# The purpose of this file is to create a CSV that aligns a given brain dataset (provided its ID) to its treatment condition.
# The output of this program is referred to when conducting future statistical analysis

# %% PACKAGE/MODULE IMPORTS
import pandas as pd
import os

# %% FUNCTION DEFINITIONS
# Create function to convert Animal IDs from index.csv (Atlas) to those utilized in the output voxelvolumes.csv
def normFileNames(old_filename):
    if ":" in old_filename:
        old_filename_split = old_filename[:-2].split("_")
    else:
        old_filename_split = old_filename.split("_")
    if len(old_filename_split[1]) > 1:
        new_filename = "".join(old_filename_split)
    else:
        new_filename = "0".join(old_filename_split)
    return "A" + new_filename + "_T1"

# %% DATA IMPORTS
# FILE: voxelvolumes.csv
vox_vol_data = pd.read_csv("/Users/nikhilgadiraju/Box Sync/Home Folder nvg6/Sharing/Bass Connections/Processed Data/Mean Intensity & Voxel Volumes/voxelvolumes.csv")
# FILE: QCLAB_AD_mice011222.csv
treatment_data = pd.read_csv("/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22/Reference Files/Absolute Files/QCLAB_AD_mice011222.csv")

# %% DATA PROCESSING

id_treatment_df = treatment_data[["Animal", "Treatment"]] # Isolate Animal ID and Treatment columns
filtered_df = id_treatment_df[id_treatment_df["Treatment"].isin(['sedentary','wheel_only','treadmill'])] # Filter to include animal IDs that have the 3 specified treatments

id_updated = list(map(lambda x: normFileNames(x), filtered_df['Animal'])) # Create list of updated IDs
id_updated_df = pd.DataFrame({'ID': id_updated, 'Treatment': list(filtered_df.iloc[:,1])}) # Create new dataframe with column 1 = updated IDs, and column 2 = Treatment

# Save dataframe as CSV
os.chdir("/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22/Reference Files/User-generated Files")
id_updated_df.to_csv('ID_Treatment.csv', encoding='utf-8')
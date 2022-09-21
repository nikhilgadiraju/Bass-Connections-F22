# VoxVol_Treatment.py
# Author: Nikhil Gadiraju

# %% PACKAGE/MODULE IMPORTS
import pandas as pd

# %% DATA IMPORTS
# FILE: voxelvolumes.csv
vox_vol_data = pd.read_csv("/Users/nikhilgadiraju/Box Sync/Home Folder nvg6/Sharing/Bass Connections/Processed Data/Mean Intensity & Voxel Volumes/voxelvolumes.csv")
# FILE: QCLAB_AD_mice011222.csv
treatment_data = pd.read_csv("/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22/Reference Files/QCLAB_AD_mice011222.csv")

# %% DATA PROCESSING
id_treatment_df = treatment_data[["Animal", "Treatment"]] # Isolate Animal ID and Treatment columns
filtered_df = id_treatment_df[id_treatment_df["Treatment"].isin(['sedentary','wheel_only','treadmill'])]


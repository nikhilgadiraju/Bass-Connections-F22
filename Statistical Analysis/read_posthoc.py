# VoxVol_Treatment.py
# Author: Nikhil Gadiraju

# %% PACKAGE/MODULE IMPORTS
import pandas as pd

# %% DATA IMPORTS
df = pd.read_csv("/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22/Statistical Analysis/posthoc_standard.csv")
df_FDR = df.sort_values(by=['FDR'], ascending=True)
print(df_FDR)
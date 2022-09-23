# VoxVol_Treatment.py
# Author: Nikhil Gadiraju

# %% PACKAGE/MODULE IMPORTS
import pandas as pd

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

# %% DATA PROCESSING
index_df = index.iloc[:,0:3]
structure_updated = list(map(lambda x,y: structure_update(x,y), index_df["Structure"], index_df["Hemisphere"]))
abbreviation_updated = list(map(lambda x,y: abbreviation_update(x,y), index_df["Abbreviation"], index_df["Hemisphere"]))
struc_abbrev = pd.DataFrame({'Abbreviation': abbreviation_updated, 'Structure': structure_updated})
namedict = struc_abbrev.to_dict()

df_FDR = df.sort_values(by=['ST'], ascending=True)
print(df_FDR)
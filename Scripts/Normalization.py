# Bass Connections Fall 2022
# Author: Nikhil Gadiraju

# %% NOTES
# The purpose of this file is to process the raw collected mice data and output CSV's of the mean intensity values and voxel counts for each brain region within a number of brain datasets

# %% NOTES
# ** The only variables that the user should have to change is the file paths found in the below section (below imports) as well as those in the 'loopscript.sh' shell script

# %% Imports
from __future__ import print_function
import os
import numpy as np
import pandas as pd
from itertools import product
import statistics
import nibabel as nib
from VoxVol_Treatment import id_updated_df

# %% Defining file paths
# MAIN file paths; (fpath_main) this is where this python script is found as well as the folders indicated in the below file paths; (fpath_dat_main) This is where all input and output data is stored
fpath_main = "/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22"
fpath_dat_main = "/Users/nikhilgadiraju/Box Sync/Home Folder nvg6/Sharing/Bass Connections"

# INPUT file paths
# Location of all data (Volumes AND Labels)
fpath_dat = fpath_dat_main + "/Data/CVN_T1_Skullstripped_Labels"
fpath_norm_csv = fpath_main + "/Reference Files/Absolute Files/normVals.csv"
fpath_atlas_csv = fpath_main + "/Reference Files/Absolute Files"
fpath_id_treatment_csv = fpath_main + "/Reference Files/User-generated Files"

# OUTPUT file paths
fpath_norm = fpath_dat_main + "/Processed Data/Water Tube-Normalized Brain Volumes"
fpath_bfc_img = fpath_dat_main + \
    "/Processed Data/Bias Field Corrected (BFC) Brain Volumes"
fpath_bfield = fpath_dat_main + "/Processed Data/BFC Bias Fields"
fpath_reg_csv = fpath_dat_main + "/Processed Data/Regional Data CSVs"
fpath_meanint_voxvol_csv = fpath_dat_main + \
    "/Processed Data/Mean Intensity & Voxel Volumes"

# %% Read input files
# Save all brain volume and their associated label files' names into a list
file_names = sorted(next(os.walk(fpath_dat), (None, None, []))[2])

# Read Tube-Normalized CSV (includes normalized filenames and mRatio values)
df = pd.read_csv(fpath_norm_csv)

norm_filenames = list(df.iloc[:, df.columns.get_loc('filename')].copy())
mratio = list(df.iloc[:, df.columns.get_loc('mRatio')].copy())

# %% Defining data label dictionary
# Data label dictionary associates brain volume file name and its respective label set file name and mRatio (stores latter two values in array)
data_filename = []
lab_filename = []

for k in range(0, len(file_names)):
    if 'label' in file_names[k]:
        lab_filename.append(file_names[k])
    else:
        data_filename.append(file_names[k])

data_lab_dict = {}

for x in range(0, len(data_filename) - 1):
    data_lab_dict[data_filename[x]] = lab_filename[x]

# %% Deleting entries in data label dictionary that lack an mRatio value in the Tube-Normalized CSV
union = []
diff = []

for k in data_lab_dict:
    for y in range(0, len(norm_filenames)):
        if norm_filenames[y][0:12] == k[0:12]:
            data_lab_dict[k] = [data_lab_dict[k], mratio[y]]
            union.append(k)

for k in data_lab_dict:
    if str(k) not in union:
        diff.append(k)

for key in diff:
    del data_lab_dict[key]

# %% Importing and Defining Atlas legend dictionary
# Atlas dictionary associates brain region abbreviation and its associates label value
os.chdir(fpath_atlas_csv)  # Location of 'CHASSSYMM3AtlasLegends.csv'
df = pd.read_csv('index.csv')  # index.csv = CHASSSYMM3AtlasLegends.csv
abbreviation = df.iloc[:, df.columns.get_loc('Abbreviation')].copy()
hemisphere = df.iloc[:, df.columns.get_loc('Hemisphere')].copy()
index = df.iloc[:, df.columns.get_loc('index2')].copy()

abb_updated = []

# Rename region abbreviations to indicate left and right hemisphere
for k in range(0, len(abbreviation)):
    if hemisphere[k] == 'Left':
        str_updated = str(abbreviation[k]) + '_L'
        abb_updated.append(str_updated)
    elif hemisphere[k] == 'Right':
        str_updated = str(abbreviation[k]) + '_R'
        abb_updated.append(str_updated)

atlas = {}

for k in range(0, len(abb_updated) - 1):
    atlas[abb_updated[k]] = index[k]

# %% Begin main for loop
for x in data_lab_dict:
    example_filename = os.path.join(fpath_dat, x)
    label_filename = os.path.join(fpath_dat, data_lab_dict[x][0])

    # Whole-Brain Data processing
    img = nib.load(example_filename)
    lab = nib.load(label_filename)

    # Normalizing brain volume data
    data = img.get_fdata()
    data = np.multiply(data, np.multiply(
        np.ones(data.shape), data_lab_dict[x][1]))

    # Save normalized brain volumes to below folder
    newname = x[0:-7] + '_norm' + '.nii.gz'
    os.chdir(fpath_norm)  # Location of normalized brain volumes
    nib.save(nib.Nifti1Image(data, img.affine), newname)

# %% Bias field Correction
# Run in local system
os.chdir(fpath_main + "/Scripts")
os.system('bash loopscript.sh')

# %% Find Regional mean values
filenames = next(os.walk(fpath_bfc_img), (None, None, []))[2]

for x in filenames:
    # Load Bias Field Corrected (BFC) volumes
    lab_fname = x[:-22] + '.nii.gz'
    img = nib.load(os.path.join(fpath_bfc_img, x))  # Location of BFC volumes
    lab = nib.load(os.path.join(fpath_dat, data_lab_dict[lab_fname][0]))

    corrected_image = img.get_fdata()

    datf = []

    for k in range(0, len(abbreviation) - 1):
        reg_mask = np.array((np.ma.array(lab.get_fdata()) ==
                            atlas.get(abb_updated[k])) * 1)  # Region Mask
        reg_iso = np.multiply(np.array(corrected_image),
                              reg_mask)  # Isolated Region (image)
        # Mean intensity value
        mean_val = np.sum(reg_iso) / np.count_nonzero(reg_mask)

        # For whole-brain proportions, change 'np.count_nonzero(reg_mask)' to 'np.count_nonzero(reg_mask)/np.count_nonzero(corrected_image)' and re-run code
        val_arr = [abb_updated[k], atlas.get(abb_updated[k]), np.count_nonzero(
            reg_mask)/np.count_nonzero(corrected_image), mean_val]
        datf.append(val_arr)

    reg_df = pd.DataFrame(datf, columns=[
                          'Structure Abbreviation', 'Index', 'Voxel Number', 'Mean intensity'])

    # Save CSV indicating regional mean values for each brain
    newname = x[0:-22] + '_df' + '.csv'
    os.chdir(fpath_reg_csv)
    reg_df.to_csv(newname, encoding='utf-8')

    print(x[0:-7] + " DONE")

# %% Generate output CSVs
os.chdir(fpath_reg_csv)

output_mi = []
output_vv = []

for x in data_lab_dict:
    csv_name = x[0:-7] + '_df' + '.csv'

    df = pd.read_csv(csv_name)
    regs_abbrev = list(
        df.iloc[:, df.columns.get_loc('Structure Abbreviation')].copy())
    vox_vol = list(df.iloc[:, df.columns.get_loc('Voxel Number')].copy())
    mean_int = list(df.iloc[:, df.columns.get_loc('Mean intensity')].copy())

    # Remove duplicates region labels
    oc_set = set()
    res = []
    for idx, val in enumerate(regs_abbrev):
        if val not in oc_set:
            oc_set.add(val)
        else:
            res.append(idx)

    for k in sorted(res, reverse=True):
        del regs_abbrev[k]
        del vox_vol[k]
        del mean_int[k]

    rowarr_mi = [x, x[0:-7]]
    rowarr_vv = [x, x[0:-7]]

    for k in range(0, len(regs_abbrev)):
        rowarr_mi.append(mean_int[k])
        rowarr_vv.append(vox_vol[k])

    output_mi.append(rowarr_mi)
    output_vv.append(rowarr_vv)

columnsarr = ['Filename', 'ID'] + regs_abbrev

out_df_mi = pd.DataFrame(output_mi, columns=columnsarr)

out_df_vv = pd.DataFrame(output_vv, columns=columnsarr)

# Remove nans
reg_del = []

for columns in out_df_mi:
    temparr = list(out_df_mi[columns])
    if all(i != i for i in temparr):
        reg_del.append(columns)

out_df_mi = out_df_mi.drop(reg_del, axis=1)

out_df_vv = out_df_vv.drop(reg_del, axis=1)

# %% Regenerate z-score CSV

output_mi = []
partial_nan = []

for x in data_lab_dict:
    csv_name = x[0:-7] + '_df' + '.csv'

    df = pd.read_csv(csv_name)
    regs_abbrev = list(
        df.iloc[:, df.columns.get_loc('Structure Abbreviation')].copy())
    vox_vol = list(df.iloc[:, df.columns.get_loc('Voxel Number')].copy())
    mean_int = list(df.iloc[:, df.columns.get_loc('Mean intensity')].copy())

    # Remove duplicates region labels
    oc_set = set()
    res = []
    for idx, val in enumerate(regs_abbrev):
        if val not in oc_set:
            oc_set.add(val)
        else:
            res.append(idx)

    for k in sorted(res, reverse=True):
        del regs_abbrev[k]
        del vox_vol[k]
        del mean_int[k]

    # Remove nan columns (ALL NaN values)
    ind_rem = []

    for k in range(0, len(regs_abbrev)):
        if regs_abbrev[k] in reg_del:
            ind_rem.append(k)

    for k in sorted(ind_rem, reverse=True):
        del regs_abbrev[k]
        del vox_vol[k]
        del mean_int[k]

    rowarr_mi = [x, x[0:-7]]

    for k in range(0, len(regs_abbrev)):
        reg_vals = list(
            out_df_mi.iloc[:, out_df_mi.columns.get_loc(regs_abbrev[k])].copy())

        # Accounting for Partial NaN columns
        reg_mask = list(np.isnan(reg_vals))

        locs = []
        vals = []

        for x in range(0, len(reg_mask)):
            if reg_mask[x] == False:
                locs.append(x)
                vals.append(reg_vals[x])

        # Calculating Z-score
        stdev = statistics.pstdev(vals)
        mean = np.mean(vals)

        if np.isnan(mean_int[k]) == False:
            z_mi = (mean_int[k] - mean) / stdev
        else:
            z_mi = np.nan
            partial_nan.append(regs_abbrev[k])

        rowarr_mi.append(mean_int[k])
        rowarr_mi.append(z_mi)

    output_mi.append(rowarr_mi)

cols = product(regs_abbrev, ['Mean_int', 'Z-score'])
columnsarr = pd.MultiIndex.from_tuples(
    [('Filename', ''), ('ID', '')] + list(cols))

# Adding treatment column
id_treatment = id_updated_df
treatment_dict = id_treatment.set_index('Modified ID')["Treatment"].to_dict()
id_vv = pd.Series(out_df_vv.loc[:, "Modified ID"])


def treatment_match(input_id):
    if str(input_id) in list(treatment_dict.keys()):
        for id, treatment in treatment_dict.items():
            if str(input_id) == str(id):
                return treatment
    else:
        return None


treatment_map = list(map(lambda x: treatment_match(x), id_vv))
treatment_align_df = pd.DataFrame({'ID': id_vv, 'Treatment': treatment_map})
out_df_vv.insert(2, 'Treatment', list(treatment_align_df.iloc[:, 1]))
out_df_vv = out_df_vv.sort_values(by=['Treatment'], ascending=False)

os.chdir(fpath_meanint_voxvol_csv)
out_df_mi = pd.DataFrame(output_mi, columns=columnsarr)
out_df_mi.to_csv('meanintensities.csv', encoding='utf-8')

out_df_vv.to_csv('voxelvolumes.csv', encoding='utf-8')

# %% Print Nan Regions
print("Nan regions:")
for k in reg_del:
    print(k)

print("\nPartial Nan regions:")
for k in list(set(partial_nan)):
    print(k)

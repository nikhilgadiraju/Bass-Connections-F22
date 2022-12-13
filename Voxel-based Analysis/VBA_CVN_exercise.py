"""Voxel-Based Morphometry on Oasis dataset
========================================
This example uses Voxel-Based Morphometry (VBM) to study the relationship
between aging, sex and gray matter density.
The data come from the `OASIS <http://www.oasis-brains.org/>`_ project.
If you use it, you need to agree with the data usage agreement available
on the website.
It has been run through a standard VBM pipeline (using SPM8 and
NewSegment) to create VBM maps, which we study here.
VBM analysis of aging
---------------------
We run a standard GLM analysis to study the association between age
and gray matter density from the VBM data. We use only 100 subjects
from the OASIS dataset to limit the memory usage.
Note that more power would be obtained from using a larger sample of subjects.
THIS IS THE ONE YOURE WORKING ON 
"""
# Authors: Bertrand Thirion, <bertrand.thirion@inria.fr>, July 2018
#          Elvis Dhomatob, <elvis.dohmatob@inria.fr>, Apr. 2014
#          Virgile Fritsch, <virgile.fritsch@inria.fr>, Apr 2014
#          Gael Varoquaux, Apr 2014


from nilearn.image import load_img, math_img
from nilearn.glm import cluster_level_inference
from scipy.stats import norm
import matplotlib as mpl
from nilearn.glm import threshold_stats_img
from nilearn.glm.second_level import SecondLevelModel
from nilearn.plotting import plot_design_matrix
from nilearn.image import resample_to_img
from bunch import bunchify

############################################################################
# Load Oasis dataset
# ------------------
import numpy as np
import gzip
import os
import glob
import pandas as pd
import itertools
from nilearn.input_data import NiftiMasker
from nilearn.image import get_data
import matplotlib.pyplot as plt
import sklearn
from nilearn import plotting
from PIL import Image as img
import nibabel as nib

excel_path = '/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22/Reference Files/User-generated Files/Internally Referenced/ID_Treatment.xlsx'
mouse_database = pd.read_excel(excel_path, sheet_name='ID_Treatment').sort_values('Treatment')
mouse_database = mouse_database[mouse_database['N-number'] != 'Blank']
mouse_database_nnum = list(mouse_database.loc[:, "N-number"])

mouse_images = []
mouse_images_folder = ['/Users/nikhilgadiraju/Box Sync/Home Folder nvg6/Sharing/Bass Connections/Data/jac_smoothed']
mifl = list(itertools.chain.from_iterable(itertools.repeat(mouse_images_folder, len(mouse_database_nnum))))
fileextension1 = ['_jac_to_MDT.nii']
fileextension = list(itertools.chain.from_iterable(itertools.repeat(fileextension1, len(mouse_database_nnum))))


###############
#gray_matter_map_filenames = gm_maps_masked
# for j in range(0,len(np.unique(mouse_database['treatment']))):
#    a = mouse_database['treatment'] == sorted(np.unique(mouse_database['treatment']))[j] #0 is iron, 1 is LPS, 2 is LPS+Iron, 3 is Saline
#    mouse_database['TreatmentNum'][a] = j
# ['sedentary'0, 'treadmill'1, 'wheel_only'2]


# code the treatement as 0 1 0 -1 to compare iron saline and saline only
#mouse_database['TreatmentNum'] = mouse_database['treatment']
# for j in range(0,len(np.unique(mouse_database['Treatment']))):
mouse_database['Treatment_type'] = mouse_database['Treatment']
a = mouse_database['Treatment_type'] == "sedentary"
mouse_database['Treatment_type'][a] = 0
b = mouse_database['Treatment_type'] == "wheel_only"
mouse_database['Treatment_type'][b] = 1
c = mouse_database['Treatment_type'] == "treadmill"
mouse_database['Treatment_type'][c] = 2


# remove 0 (liuke sed if it is zero)
# mouse_database = mouse_database[(mouse_database['Treatment'] == "HFD")]


# mouse_database_1=mouse_database.dropna(axis=0)

# np.dropna(mouse_database)

###################


# Diet=mouse_database['Diet']
# mouse_database[ isin mouse_database['DWI']]
# str(Timepoints)
# mouse_names = list(mouse_database.DWI)
# index = [idx for idx in range(len(mouse_names))
#          if 'N' in str(mouse_names[idx])]
# mouse_names = [mouse_names[i] for i in index]
# #mouse_names = mouse_database.DWI
# mouse_names = [a[0:6] for a in mouse_names]

# mouse_database = mouse_database.iloc[np.array(index), :]


#mouse_names = mouse_names[4:80:3]
mouse_paths = []

#strings = ["%d" % Timepoints for Timepoints in Timepoints]
mouse_names = []

for num in mouse_database_nnum:
    mouse_names.append("/s1p5{}".format(num))

filepathmousenames = [''.join(z)
                      for z in zip(mifl, mouse_names, fileextension)]
#filepathtimepoints = [''.join(z) for z in zip(filepathmousenames, fileextension)]
#fullfilepaths = [''.join(z) for z in zip(filepathmousenames,fileextension)]
gray_matter_map_filenames = filepathmousenames


# remove sed
# mouse_database=mouse_database[mouse_database['TreatmentNum']!=0]


#mouse_database['GenotypeNum'] = mouse_database['Genotype']
# for k in range(0,len(np.unique(mouse_database['Genotype']))):
#    b = mouse_database['Genotype'] == sorted(np.unique(mouse_database['Genotype']))[k] #0 is 5xFAD, 1 is WT
#    mouse_database['GenotypeNum'][b] = k

#Diet = mouse_database['DietNum'].astype(float)

#treatment_subset = treatment[4:80:3]
#treatment = mouse_database[4]

###############################################################################
# Sex is encoded as 'M' or 'F'. Hence, we make it a binary variable.
#genotype = mouse_database['Genotype'] == b'F'
#genotype[102:504] = bool(123)
#genotype = mouse_database['GenotypeNum'].astype(float)
#genotype_subset = genotype[4:80:3]


#genotypestring = genotype.to_string()
#genotypesort = genotype.iloc[0:503:3]

###############################################################################
# Print basic information on the dataset.
print('First gray-matter anatomy image (3D) is located at: %s' %
      gray_matter_map_filenames[0])  # 3D data

###############################################################################
# Get a mask image: A mask of the  cortex of the ICBM template.
gm_mask = nib.load(
    '/Users/nikhilgadiraju/Box Sync/Home Folder nvg6/Sharing/Bass Connections/Data/median_images/MDT_mask_e3.nii')
bgimage = nib.load(
    '/Users/nikhilgadiraju/Box Sync/Home Folder nvg6/Sharing/Bass Connections/Data/median_images/MDT_fa.nii.gz')

###############################################################################
# Resample the images, since this mask has a different resolution.
mask_img = resample_to_img(
    gm_mask, gray_matter_map_filenames[0], interpolation='nearest')
# mask_img=gm_mask

#############################################################################
# Analyse data
# ------------
#
# First, we create an adequate design matrix with three columns: 'age',
# 'sex', 'intercept'.
n_subjects = len(mouse_database)  # more subjects requires more memory
intercept = np.ones(n_subjects)
#intercept_subset = intercept[4:80:3]


exercise = mouse_database['Treatment']
exercise[exercise == "sedentary"] = 0.0
exercise[exercise == "treadmill"] = 1.0
exercise[exercise == "wheel_only"] = -1.0
exercise = exercise.astype(float)
mouse_database['Treatment_num'] = exercise


design_matrix = pd.DataFrame(np.vstack((exercise, intercept)).T, columns=[
                             'exercise', 'intercept'])
#design_matrix = pd.DataFrame(np.vstack((treatment, intercept)).T, columns=['treatment', 'intercept'])

#############################################################################


# Let's plot the design matrix.

ax = plot_design_matrix(design_matrix)
ax.set_title('Second level design matrix', fontsize=12)
ax.set_ylabel('maps')

##########################################################################
# Next, we specify and fit the second-level model when loading the data and
# also smooth a little bit to improve statistical behavior.

#second_level_model = SecondLevelModel(mask_img=mask_img)
second_level_model = SecondLevelModel(smoothing_fwhm=.1, mask_img=mask_img)

second_level_model.fit(gray_matter_map_filenames, design_matrix=design_matrix)

##########################################################################
# Estimating the contrast is very simple. We can just provide the column
# name of the design matrix.
z_map = second_level_model.compute_contrast(
    second_level_contrast=[1, 0], output_type='z_score', )

###########################################################################
# We threshold the second level contrast at uncorrected p < 0.001 and plot it.


###
myalpha = 0.2


# size, threshold = threshold_stats_img(z_map, alpha=myalpha, height_control='fdr',
#                                      two_sided='false')


size, threshold = threshold_stats_img(
    z_map, alpha=myalpha, height_control='fdr')


print('The FDR=.05-corrected threshold is: %.3g' % threshold)
#nib.save(size, "size_sed_tread_0.05.nii.gz")


#condition_effect = np.hstack((  [0]* n_subjects  ,[1] * n_subjects, [0] * n_subjects , [- 1] * n_subjects))


#z_map = second_level_model.compute_contrast(second_level_contrast='treatment',output_type='z_score')
# z_map = second_level_model.compute_contrast(second_level_contrast=[1, 0],
#                                           output_type='z_score',second_level_stat_type='t' )


mpl.rcParams['figure.dpi'] = 1000

p_val = 0.05
p001_uncorrected = norm.isf(p_val)

proportion_true_discoveries_img = cluster_level_inference(
    z_map, threshold=p001_uncorrected, alpha=myalpha)

fig = plotting.plot_stat_map(
    size, threshold=p001_uncorrected, colorbar=True, display_mode='x',
    bg_img=bgimage, cut_coords=np.linspace(-5, 5, 12, endpoint=True))


#mask2 = math_img ('img < 0', img=size, threshold=threshold)

#from nilearn.image import load_img, math_img
#mask2 = math_img('img > 0', img=size)

# dataobj=z_map.dataobj
# dataobj[dataobj<0]=0
# mask2=z_map
# mask2.dataobj=dataobj

# fig = plotting.plot_stat_map(
#     dataobj, threshold=p001_uncorrected, colorbar=True, display_mode='x',
#     bg_img=bgimage, cut_coords=np.linspace(5, -5, 12, endpoint=True))


#fig.savefig('sed_vs_tread_0.05.png', dpi=1000)


myalpha = 1

size, threshold = threshold_stats_img(
    z_map, alpha=myalpha, height_control='fdr')
print('The FDR=.05-corrected threshold is: %.3g' % threshold)
#nib.save(size, "size_sed_tread_0.1.nii.gz")


#condition_effect = np.hstack((  [0]* n_subjects  ,[1] * n_subjects, [0] * n_subjects , [- 1] * n_subjects))


#z_map = second_level_model.compute_contrast(second_level_contrast='treatment',output_type='z_score')
# z_map = second_level_model.compute_contrast(second_level_contrast=[1, 0],
#                                           output_type='z_score',second_level_stat_type='t' )


mpl.rcParams['figure.dpi'] = 1000

p_val = 0.05
p001_uncorrected = norm.isf(p_val)

proportion_true_discoveries_img = cluster_level_inference(
    z_map, threshold=p001_uncorrected, alpha=myalpha)

fig = plotting.plot_stat_map(
    size, threshold=p001_uncorrected, colorbar=True, display_mode='x',
    bg_img=bgimage, cut_coords=np.linspace(5, -5, 12, endpoint=True))
#fig.savefig('sed_vs_tread_0.1.png', dpi=1000)


# filter fa mask for white matter analysis
# fa_mask_path = '/Users/nikhilgadiraju/Box Sync/Home Folder nvg6/Sharing/Bass Connections/Data/median_images/MDT_fa.nii.gz'

# fa_nifti = nib.load(fa_mask_path)

# mask = math_img('img > .3', img=fa_nifti)
# np.min(mask.get_data())
# np.max(mask.get_data())
# np.median(mask.get_data())

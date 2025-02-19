# Pipeline Steps

## Step 1: Parcellation
Parcellation of functional images (preprocessed with pipeline *preprocessed_rp_mo_csf_wm_nui*) registered to MNI using two atlases: *Schaefer* [(Github page)](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/MNI)[(Reference)](10.1093/cercor/bhx179), using 100 parcels, with direct correspondence to the intrinsic networks identified by Yeo et al. [(Reference)](10.1152/jn.00338.2011); *AAL116* with 116 parcels including subcortical and cerebellum regions [(Labels)](https://cran.r-project.org/web/packages/brainGraph/index.html) [(Reference)](https://doi.org/10.1006/nimg.2001.0978). 


We are currently using the function `parcellation.m`, with files:

- ../FC/files/atlases/Schaefer

- ../FC/files/atlases/AALAAL116_2mm_woSchaefer.nii.gz


```matlab
PARCELLATION

Parcellate functional image based on specified atlas.

[all_av_allvoxels, all_av_nonzero, percentage_included] = PARCELLATION(atlasdir, atlas, funcimg, plt)
```
### Inputs:

- `atlasdir`: string with the directory containing atlas files.
- `atlas`: string with the name of the atlas to use ('AAL90', 'AAL116', 'Desikan', 'HO', 'HOcort', 'HOsubcort', 'SchaeferH7100', 'SchaeferH17100', 'CDK50', 'CDK100', 'CDK500').
- `funcimg`: 4-D Functional image data.
- `plt`: Boolean flag indicating whether to plot the atlas (true or false).

### Outputs:

- `all_av_allvoxels`: Matrix containing mean intensity values for all voxels in each region (dim = #parcels x #volumes).
- `all_av_nonzero`: Matrix containing mean intensity values for non-zero voxels in each region (dim = #parcels x #volumes).
- `percentage_included`: Percentage of included voxels for each region (dim = #parcels x 1).

Parcellated data is organized by subject and session, stored in `data/parcellated_data`.

## Step 1a: Reorganization of parcels to get non-zero SchaeferSubCRB7100
Reorganization of Schaefer's parcellation to be grouped by network instead of hemisphere, using the `utils` script `reorganize_Schaefer_H2N.m`, which generates a modified labels file stored as `files/atlases/Schaefer/Schaefer2018_100Parcels_7Networks_order_new.txt` and a sorted index vector stored as `files/atlases/Schaefer/Schaefer/Schaefer_H2N_ind.mat`. The indices are used to reorganize the timecourses of the 100 parcels. Additionally, the timecourses of several subcortical and cerebellum parcels from AAL116 are appended to the reorganized files:
- hippocamppus - 37, 38 (L, R)
- amygdala - 41, 42 (L, R)
- caudate - 71, 72 (L, R)
- putamen - 73, 74 (L, R)
- pallidum - 75, 76 (L, R)
- thalamus - 77, 78 (L, R)
- cerebellum - 91-116 (9 L, R + vermis)

The same procedure was applied to the previously computed `percentage_included` files.

### Removal of parcels with non-zero voxels
Since the FOV did not allow to fully retrieve all the relevant structures for the whole sample (specially cerebellum, but also some frontal regions), the number of non-zero voxels contributing to each parcel was analyzed in `visualize_parcellation_nullROI`.

![Non-zero voxels](parcellation_nullROI_percentage_included_criteria.png) 

Figure details may be checked on `visualization.md`. Based on the results, the parcels from which more than 50% was zero for at least one subject (considering <mark style="background-color: lightyellow;">all subjects with available functional data</mark>). The indices of excluded regions (relative to *SchaeferSubCRB7100*) were stored in `files/atlases/parcellation_nullROI_parcel2exclude_SchaeferSubCRB7100_idx` and correspond to 4 regions, bilaterally:
- Crus II L (SchaeferSubCRB7100 - 115, AAL116 - 93)
- Crus II R (SchaeferSubCRB7100 - 115, AAL116 - 94)
- Cerebellum7b L (SchaeferSubCRB7100 - 123, AAL116 - 101)
- Cerebellum7b R (SchaeferSubCRB7100 - 124, AAL116 - 102)
- Cerebellum8 L (SchaeferSubCRB7100 - 125, AAL116 - 103)
- Cerebellum8 R (SchaeferSubCRB7100 - 126, AAL116 - 104)
- Cerebellum9 L (SchaeferSubCRB7100 - 127, AAL116 - 105)
- Cerebellum9 R (SchaeferSubCRB7100 - 128, AAL116 - 106)

All these parcels originally belonged to AAL116, no parcel from Schaefer's atlas was excluded. The final parcellated and reorganized files are stored in `data/parcellated_data/BOLD<pipeline>-parcellatedSchaeferSubCRB7100`.

## Step 2: Static functional connectivity (sFC) using Pearson Correlation
This MATLAB script performs functional connectivity (FC) analysis using Pearson correlation, using preprocessed and parcellated data for multiple subjects and sessions, calculating FC matrices, and extracting leading eigenvectors (to match the structure of the dFC results produced by [LEiDA](https://github.com/juanitacabral/LEiDA/blob/master/LEiDA/LEiDA.m)).

### Inputs:
- Parcellated data

### 1. Bandpass Filtering
- Applies a bandpass filter (2nd order Butterworth) to BOLD data, from 0.01 to 0.1 Hz.

### 2. FC Matrix Calculation
- Filters demeaned BOLD data for each seed (brain area).
- Computes Pearson correlation matrices (FC) for the filtered BOLD data for each subject and session, based on the configuration list.
- Extracts the leading eigenvector and calculates the variance explained.
- Logs information about the processing steps for each subject and session.
- Saves FC data, leading eigenvectors, variance explained, and other relevant information for the entire sample.

The FC matrices are stored as `FC_all`, in a #parcels x #parcels x #subjects/sessions matrix and the order of the subjects/sessions is given by the `filenames` variable. The number of parcels is stored as `N_areas`.

## Step 3: Differential Identifiability - template matrices
This MATLAB script creates templates for multilevel differential identifiability using the subject, session and group information for the MIG_N2Treat project. There are two types of templates, for differential identifiability and for grouping correlation values. The first are coded with binary values and the second are coded with integer numbers starting with 1. In both cases, irrelevant entries are coded as 0.5. Considering that each sample corresponds to one subject in one session, the template matrices are squared and each dimension corresponds to the total number of samples.

### Differential Identifiability
- Iwithin: correlation within a category of a certain level - coded as 1
- Iothers: correlation between different categories of a certain level - coded as 0
- (other parts of the matrix are coded as 0.5)

These templates were computed for:
- same session - `within_session`
- same subject - `within_subject`
- same group - `within_group`
- different groups, same menstrual cycle phase - `between_group_within_session`
- same menstrual cycle phase (regardless of the group) - `within_menstrual_session`

### Identifiability correlation
- within session - `wsession_corr`: 1 - preictal; 2 - ictal; 3 - postictal; 4 - interictal; 5 - premenstrual; 6 - midcycle
- within group between sessions - `wgroup_bsession_corr`: 1 - preictal-ictal; 2 - preictal-postictal; 3 - preictal-interictal; 4 - ictal-postictal; 5 - ictal-interictal; 6 - postictal-interictal; 7 - premenstrual-midcycle
- between group - `bgroup_corr`: 1 - preictal-premenstrual; 2 - ictal-premenstrual; 3 - postictal-premenstrual; 4 -  interictal-midcycle
- within group - `wgroup_corr`: 1 - migraine patients; 2 - healthy controls
- within subject between sessions - `wsubject_bsession_corr`: `within_subject.*wgroup_bsession_corr`

For all cases, other parts of the matrix are coded as 0.5.

The files are stored in `files/template_matrices`;

## Step 4: Differential Identifiability - computation

### Inputs:

- **pipeline**: String specifying the preprocessing pipeline used (e.g., 'preprocessed_rp_mo_csf_wm_nui').
- **atlas**: String specifying the brain atlas used for parcellation (e.g., 'SchaeferSubCRB7100').
- **networks**: Cell array containing the names of different functional networks.

### 1. Loading Pearson Correlation Data

- Loads precomputed Pearson correlation variables (**`FC_all`**) from a file based on the specified atlas and preprocessing pipeline.
- Extracts information such as the number of areas, subjects, filenames, and the correlation matrix.

### 2. FC Vectorization

- Converts the correlation matrices (**`FC_all`**) into vectorized lower triangular matrices (**`FCvec`**) without the diagonal for further analysis.
- Demeans the data by subtracting the mean.

### 3. Networks Mask

- Defines a mask (**`matrix_all`**) for different functional networks based on predefined areas.
- Creates a network template (**`matrix_FC_all`**) by summing the network masks.
- Vectorizes the lower triangular matrix of the network template (**`network_mask_vec`**).

### 4. PCA reconstruction + Multilevel Identifiability matrix

Performs PCA on the demeaned functional connectivity data (**`group_data`**) and reconstructs the functional connectivity matrices using different numbers of principal components (**`components`**).


- **PCA Matrix Calculation:**
    - Calculates the PCA matrix (**`PCA_matrix`**) using a subset of principal components and back the original mean to the PCA-reconstructed data (**`group_PCA_recon`**).
    - Stores the concatenated vectorized FC for each number of principal components (**`FCvec_PCA_recon`**).
    - Stores the FC matrices (**`FC_PCA_recon`**) for each subject after PCA reconstruction.
    - Calculates the explained variance (**`explained`**) for each number of components.
    - Stores principal component coefficients (**`coeffs`**), scores (**`score`**), and explained variance.
- **Explained Variance:**
    - Computes the cumulative explained variance (**`exp_var`**) up to the current number of principal components.
    - Stores the explained variance for each number of principal components (**`PCA_VE`**).
- **Identifiability Matrix for Networks:**
    - For each functional network or the whole brain (**`WB`**):
        - If not the last network (**`n <= nr_networks-1`**), extracts data corresponding to that network from the PCA-reconstructed group data.
        - If the last network (**`n == nr_networks`**), uses the entire reconstructed group data.
    - Calculates the Pearson correlation between vectorized FC matrices within the selected network/WB (**`corr_val`**).
    - Stores the identifiability matrix for each network/WB and number of principal components (**`corr_PCA_recon_networks`**).

### 5. Computing Differential Identifiability
Using the multilevel identifiability matrix (`corr_matrix`) and the desired identifiability template matrix as input, the function `diff_identifiability` provides the differential identifiability values: Iself/Iwithin, Iothers and Idiff.

The differential identifiability values for each template are stored in `data/results/DI/DI_values`. These are stored in a matrix with dimensions 3 x nr_networks x nr_components (considering 10 networks, including WB, and 68 components). 
- `within_session_DI`
- `within_subject_DI`
- `within_group_DI`
- `within_menstrual_session_DI`
- `between_group_within_session_DI`

### 6. Computing average correlation
Using the multilevel identifiability matrix (`corr_matrix`) and the desired correlation template matrix as input, the function `diff_identifiability` provides the average correlation values.

The average correlation values for each template are stored in `data/results/DI/DI_corr`. These are stored in a cell array with dimensions nr_networks x nr_components (considering 10 networks, including WB, and 68 components). Each element of the cell array will correspond to a matrix with a size that depends on the template.
- `corr_val_wgroup`
- `corr_val_wgroup_bsession`
- `corr_val_wsubject_bsession`
- `corr_val_wsession`
- `corr_val_bgroup`

## Step 5: Network-Based Statistic (NBS) Stats

### Inputs:

- nr_patients = 10;
- nr_controls = 14;
- nr_sessions_patients = 4;
- nr_sessions_controls = 2;

### 1. Create and store design matrices
Create NBS design matrices for one-sample t-test for patients and controls and two-sample t-test for between group comparison:
- patients - `design_patients.txt`
- controls - `design_controls.txt`
- patients vs controls - `design_patients_controls.txt`

The .txt files are stored in `files/nbs/design`.

### 2. Create and store .txt files with node coordinates and labels
Use the file coord_SchaeferSubCRB7100_130.txt from `files/atlases` to get the relevant information and split the coordinates and labels into different files.

The .txt files are stored in `files/nbs/nodes`.

### 3. Run NBS
- Create exchange blocks for patients, controls and patients vs controls

- Get FC data for the corresponding number of PCs and build matrices with concatenated FC data for each comparison
    - comparisons = {'preic-ic', 'preic-postic', 'preic-interic', 'postic-ic', 'ic-interic', 'postic-interic', 'prem-mid', 'preic-prem', 'ic-prem', 'postic-prem', 'interic-mid'};

- Run NBS comparison for each method, threshold, contrast type and comparison:
    - cluster threshold = [2, 3.1, 4, 5, 6]
    - methods = {'Extent','Intensity'};
    - contrast types = [1, -1], used to create contrasts for the corresponding comparison
    - Fixed parameters:   
        - alpha = 0.05
        - permutations = 5000
        - design_*.txt
        - coordinates = 'coord.txt'
        - labels = 'labels.txt'

- Results are stored in a structure array called `nbs`, with several fields: GLM, NBS, STATS and UI. The subfield `nbs.NBS.con_mat` contains a sparse array which describes the node pairs connected by edges considered significant. The .mat files are stored in `data/results/nbs-stats/stats-<nComp>PCs`.

## Step 6: Network-Based Statistic (NBS) Analysis

This script provides information about the significant differences found using NBS, by summarizing the number (or percentage) of significant edges for each network pair (either the whole network or split by hemisphere). 

### Inputs:
- components = [19, 68]
- thresholds = [3.1, 4];
- comparisons = {'preic-interic', 'ic-interic', 'postic-interic', 'prem-mid', 'ic-prem', 'postic-prem', 'preic-ic'};
- contrast_types = {'-1', '-1', '-1', '-1', '1', '1', '-1'};
- method = 'Extent';
- nr_networks = 9;
- areas = [1 13 1; 14 37 2; 38 52 3; 53 57 4; 58 69 5; 70 83 6; 84 100 7; 101 112 8; 113 130 9];

It first loads the NBS statistics for comparisons with significant results and compute the count, percentage and area percentage for significant edges. The analysis is performed for each combination of threshold and number of PCs and a structure array `sigmetrics` is used to store the computed metrics for all comparisons, including the fields:
- count_sig_FC: number of significant edges
- perc_sig_FC: percentage of significant edges considering the total number of significant edges
- percarea_sig_FC: percentage of significant edges considering the total number of possible edges for that area
- countH_sig_FC: number of significant edges, split by hemisphere
- percH_sig_FC: percentage of significant edges from the total of significant ones, split per hemisphere

The inputs that were used to run the analysis are also stored. 

Note: nbs_68PCs_postic-prem_contrast1_method-Extent_cluster-4 is not available, as no significant differences were found, so it was not used.

## Step 7: Network Metrics

### Inputs:

- **`nr_areas`**: Number of areas in the network.
- **`metric_names`**: Cell array containing names of local measures to compute (**`'nodedegree'`**, **`'clustercoef'`**, **`'betweencentr'`**).
- **`comparisons`**: Cell array specifying different comparisons in the network.
- **`contrast_types`**: Cell array indicating the type of contrast for each comparison.
- **`method`**: Method used for the analysis (e.g., **`'Extent'`**).
- **`dim`**: Dimensionality information for the analysis.
This script performs two main tasks:

### 1. **NBS Significant Edges:**
- Computes local (**`'nodedegree'`**, **`'clustercoef'`**, **`'betweencentr'`**).
- Organizes the results based on the dimensionality of the analysis (areas-subjects for nodal measures, subjects for global measures).
- Stores the results as a .mat file.

### 2. **Functional Connectivity after PCA Reconstruction:**
- Loads precomputed functional connectivity data from PCA reconstruction.
- Applies thresholding to retain the top 20% strongest connections.
- Computes both local (**`'nodedegree'`**, **`'clustercoef'`**, **`'betweencentr'`**) and global network measures (**`'charpathlength'`**, **`'globaleff'`**, **`'transitivity'`**).
- Organizes the results based on the dimensionality of the analysis (areas-subjects for nodal measures, subjects for global measures).
- Saves the results for further analysis.

The script iterates over different thresholds and components, computes local measures, and saves the results in a structure array `network`. For NBS, the results are stored in `data/results/network_metrics/metrics_nbs-<nComp>PCs` and for the functional connectivity matrices, they are stored in `data/results/network_metrics/metrics_nbs-<nComp>PCs`.

## Step 8: Clinical Features
This script is used to assess the relationship between migraine patients clinical features and the results from the analysis of identifiability metrics as well as FC fingerprints.

### Clinical data
The clinical data is loaded from .mat files (one per group) stored inside the folder `mign2treat_sample`. The information is relative to all patients and controls, so it has to be extracted specifically for the sample used for FC analysis (only the patients that completed the protocol and excluding control025 due to the truncation of the resting state data). 
Although the clinical variables concern only the patients, they are also plotted for controls, considering null values, for comparison.    

- Migraine Onset (years)
- Migraine Disease Duration (years)
- Attack Frequency (number per month)
- Migraine Duration (hours)
- Migraine Pain Intensity (score)

### Association with identifiability metrics
- within subject, between sessions avg - 1 value per subject
- Mig vs HC - 4 values per subject (each session)
- Mig vs HC - 1 value per subject (avg)

### Association with average FC for significant edges using NBS 
- mask ict vs pm - 1 value per subject
- mask post vs pm - 1 value per subject




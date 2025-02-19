# Visualization Guide

## Introduction
Brief introduction to visualization scripts and their purpose.

## Available Scripts
### 1. `visualize_qc_carpetplot`: 

**Figures Path**: *../FC/figures/carpet_plots* 

The given MATLAB code is part of the data visualization for quality control.  

### 2. `visualize_parcellation_nullROI.m`:

**Figures Path**: *../FC/figures/parcellation_nullROI* 

The script assesses the percentage of inclusion for each parcel and identifying parcels that do not meet specific criteria across all subjects.

- Two figures are generated to visualize the percentage of inclusion for each parcel across data samples.
- The first figure presents a scatter plot of percentage inclusion for each parcel.
- The second figure is an imagesc plot providing an overview of inclusion percentages for all parcels across data samples.

#### Inclusion Criteria Analysis

- Calculates the number of data samples where inclusion percentages meet specific criteria (100%, >90%, >50%).
- Generates bar plots to visualize the distribution of inclusion percentages across parcels.
- Determines parcels that do not meet the inclusion criteria for all subjects.
- Generates a list of parcel indices to exclude and stores it in a file for further use.

**Figures Path**: *../FC/figures/parcellation_nullROI* 

### 3. `visualize_DI_templates.m`

**Figures Path**: *../FC/figures/DI_templates* 

The script produces figures of template matrices coded with integer numbers and the value 0.5, which corresponds to entries that should be ignored (upper triangular matrix, due to symmetry).

- `DI_templates_matrices_Idiff.png`: Iself values coded as 1 and Iothers values coded as 0, for several identifiability divisions: group, subject, session, cycle and group-cycle (binary values).

![DI_templates_matrices_Idiff](DI_templates_matrices_Idiff.png) 

- `DI_templates_matrices_correlation.png`: Values coded according to different divisions: within-group, within-group between-sessions, within-subject between-sessions, within-session and between groups.

![DI_templates_matrices_correlation](DI_templates_matrices_correlation.png) 

### 5. `visualize_DI_FC_PCA.m`

**Figures Path**: *../FC/figures/DI_FC_PCA* 

The script produces figures of FC matrices for the whole sample or separated by patients and controls (raw and PCA reconstructed). Results are presented for the whole brain and several networks.

### 6. `visualize_DI_Imatrix_PCA.m`

**Figures Path**: *../FC/figures/DI_Imatrix_PCA* 

The script produces figures of identifiability matrices (raw and PCA reconstructed). Results are presented for the whole brain and several networks.

- `DI_FC-Imatrix_PCA_networks-WB.png`: identifiability matrix as a function of the number of PCs used for data reconstruction. The identifiability matrices were computed from individual FC data concerning each one of the networks and also the whole brain.

### 7. `visualize_DI_values.m`

**Figures Path**: *../FC/figures/DI_values* 

The script produces figures of raw and PCA reconstructed identifiability matrices, and Idiff, Iself and Iothers as a function of the number of principal components of the reconstruction.  Results are presented for the whole brain and several networks.

- `DI_values_Iself-Iothers_network.png`:
Idiff and Iself vs Iothers, as a function of the number of PCs.

- `DI_values_histogram-gif_network.gif`:
Original vs PCA reconstructed matrix, Idiff and Correlation distribution (histogram Iself vs Iothers) as a function of the number of PCs, for group, subject, session, cycle and group-cycle.

- `DI_values_BrainFunctionWorkshop-panelA_Idiff.png`:
Figure presented at the Brain Function Workshop showing the Original identifiability matrix and Idiff values as a function of the number of PCs, for within-subject, within-session and within-group (Panel A).

### 8. `visualize_DI_correlation.m`

**Figures Path**: *../FC/figures/DI_correlation* 

The script produces figures of correlation values as a function of the number of principal components of the reconstruction, for several data subdivisions. Results are presented for the whole brain and several networks.

- `DI_correlation_boxplots-gif_network.gif`: Original vs PCA reconstructed matrix, Idiff (for group, subject, session, cycle and group-cycle) and Correlation boxplots as a function of the number of PCs, for within-group, within-group between-sessions, within-subject between-sessions, within-session and between groups. 

- `DI_correlation_BrainFunctionWorkshop-panelB_Correlation.png`: Figure presented at the Brain Function Workshop showing correlation values (points selected from the identifiability matrix) for within-subject, within-session and within-group (Panel B). For within-subject, in the case of MIG, the values correspond to all possible session combinations.

### Linked Functions
-  `plot_bp_pnts.m`: Plot scatter over boxplot to show individual points

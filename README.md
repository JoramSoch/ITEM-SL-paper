# ITEM-SL-paper

**MATAB Code for the ITEM-SL Paper, submitted to Imaging Neuroscience**

This code belongs to the paper on "Searchlight-based trial-wise fMRI decoding in the presence of trial-by-trial correlations" by Joram Soch, publicly available from *bioRxiv* and currently under review in *Imaging Neuroscience*. It consists of two sub-folders, "Simulation" (Section 3, Figures 3/4/7 in the paper) and "Application" (Section 4, Figures 5/6 in the paper).

- Preprint: https://www.biorxiv.org/content/10.1101/2023.12.05.570090
- Data: https://openneuro.org/datasets/ds002013 (Version 1.0.3)
- Code: https://github.com/JoramSoch/ITEM-SL-paper
- Toolbox: https://github.com/JoramSoch/ITEM


### Requirements

This code was developed and run using the following software:
- [MATLAB R2021a](https://de.mathworks.com/help/matlab/release-notes-R2021a.html) (Version 9.10)
- [MATLAB Statistics Toolbox](https://de.mathworks.com/products/statistics.html) (Version 12.1)
- [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) (Revision 7771 as of 13/01/2020)
- [ITEM Toolbox](https://github.com/JoramSoch/ITEM) (as on GitHub)
- [SPM Helpers](https://github.com/JoramSoch/spm_helper) (as on GitHub)
- [GLMsingle Toolbox](https://github.com/cvnlab/GLMsingle) (Version 1.2 or later)
- [fracridge Scripts](https://github.com/nrdg/fracridge) (Version 2.0 or later)


### Simulation

For re-running analyses of the simulation study, you need to perform the following steps:
1. Clone this repository into some folder on your computer.
2. Run the scripts `Simulation_A.m` and `Simulation_B.m` located in the sub-folder "Simulation". This performs the simulations as reported in Section 3 and Appendix A of the paper and should create results files called `Simulation_A.mat` and `Simulation_B.mat`. They also generate a number of figures which are only for diagnosis and can be closed.
3. Run the script `Figure_SimRes.m` to reproduce Figure 3 from the paper.
4. Run the script `Figure_SimResTrl.m` to reproduce Figure 4 from the paper.
5. Run the script `Figure_ParaMod.m` to reproduce Figure 7 from the paper.

Note: The directory contains random number generators for the uniform distribution (`MD_unirnd.m`) and the matrix-normal distribution (`MD_matnrnd.m`) which only rely on the built-in MATLAB functions (`rand`, `randn`). Therefore, MATLAB's Statistics and Machine Learning Toolbox is only required for SVM training (`svmtrain`) and testing (`svmpredict`) based on LS-A and LS-S estimates (cf. Step 3c in [Simulation A](https://github.com/JoramSoch/ITEM-SL-paper/blob/main/Simulation/Simulation_A.m#L297-L343) and [Simulation B](https://github.com/JoramSoch/ITEM-SL-paper/blob/main/Simulation/Simulation_B.m#L314-L367)).

To investigate a null scenario without condition decodability, set `r = 0` in [line 39](https://github.com/JoramSoch/ITEM-SL-paper/blob/main/Simulation/Simulation_A.m#L39) of `Simulation_A.m`. This will set the proportion of voxels with information to zero, such that there are no multivariate differences between the two experimental conditions. In this case, decoding accuracies should be symmetrically distributed around the chance level (see Section 3.2 and cf. Figure 4).


### Application

For re-running analyses of the empirical data, you need to perform the following steps:
1. Create a folder on your computer that hereafter is referred to as the "study directory". Within the study directory, create sub-folders called "data" and "stats".
2. Clone the GitHub repository [OpenNeuroDatasets/ds002013](https://github.com/OpenNeuroDatasets/ds002013) into a folder on your computer that hereafter is referred to as the "tools directory".
3. Enter the study directory from step 1 into [line 11](https://github.com/JoramSoch/ITEM-SL-paper/blob/main/Application/project_directories.m#L11) of `project directories.m`.
4. Enter the tools directory from step 2 into [line 14](https://github.com/JoramSoch/ITEM-SL-paper/blob/main/Application/download_dataset.m#L14) of `download_dataset.m`.
5. Now you are ready to run the main analysis script `analyses_ITEM_SL.m`. Note that this script is divided into sub-sections which can be run step-wise by commenting anything else ([Step 1](https://github.com/JoramSoch/ITEM-SL-paper/blob/main/Application/analyses_ITEM_SL.m#L41-L70) should always be uncommented). Ideally, run the code step by step to ensure that each stage of the analysis succeeds.

When data analysis has finished, the following directories should be in the "stats" sub-folder of the study directory:
- `glms-item_glm-full_ITEM_sects-all_SL-6mm_wcvCC/`
- `glms-item_glm-full_LS_A_sects-all_SL-6mm_wcvCC/`
- `glms-item_glm-full_LS_S_sects-all_SL-6mm_wcvCC/`
- `glms-item_glm-full_GLM_single_sects-all_SL-6mm_wcvCC/`

These directories represent group-level SPM analyses of voxel-wise cross-validated correlation coefficients (`cvCC`) from searchlight-based decoding of all sectors in the visual field (`sects-all`) using a searchlight radius of 6 mm (`SL-6mm`) based on either ITEM, LS-A, LS-S or GLMsingle (`ITEM`, `LS_A`, `LS_S`, `GLM_single`). In each directory, contrast images 60 to 79 (e.g. `con_0060_L_FWE_0.05_0.nii` or `con_0064_E1_unc_0.001_10.nii`) are those which are presented in the paper (see Figures 5B and 6).


### Bonus: Graphical Abstract

<img src="https://raw.githubusercontent.com/JoramSoch/ITEM-SL-paper/main/Figure_GA.png" alt="Graphical Abstract" width=1000>

This graphical abstract illustrates the core idea of the paper: When multiplying the trial-wise design matrix with itself, weighted by the (inverse of the) scan-by-scan covariance matrix, this results in the (inverse of the) trial-by-trial covariance matrix which describes the distribution of the trial-wise parameter estimates. Click [here](https://github.com/JoramSoch/ITEM-paper-SL/blob/main/Figure_GA.pdf) for a PDF version.
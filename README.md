# g_ci_spm
Matlab function for the estimation of effect size g and its confidence interval from SPM t maps

Use by es_ci_spm(t_map,con,X,mask_img,confLevel,out_name) with

t_map:      string with the name of a SPM results nifti file with t values (usually something like 'spmT_000X.nii');
con:        contrast vector used to estimate the t map;
X:          filtered and pre-whitened SPM design matrix (This should be SPM.xX.xKXs.X in the respective SPM.mat file. Please note that SPM is mean centering covariates); 
mask_img:   string with the name of the SPM mask file for the analysis (usually 'mask.nii');
confLevel:  confidence level of the estimated condfidence interval (usually something like .90 or .95);
out_name:   string with a prefix to add to the name of the results nifti files 

The results are saved as three nifti files [out_name '_g.nii'] containing the g values map, and [out_name '_g_ci_l.nii'] and [out_name '_g_ci_u.nii'] containing the lower and upper CI limit maps respectively.

Please note that the calculations are computationally expensive and can take up to several hours.

This function is dependent on SPM 12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) and the Measures of Effect Size (MES) toolbox (https://github.com/hhentschke/measures-of-effect-size-toolbox).

Examples.zip contains the code for the simulations in Figure 1 and the second level maps and necessary information to run g_ci_spm for the data presented in Figures 2-5 of the manuscript.

If this function is of help for your own work, please cite Gerchen, M.F., Kirsch, P., & Feld, G.B. (2021) Brain-Wide Inferiority and Equivalence Tests in fMRI Group Analyses: Selected Applications. Human Brain Mapping. DOI: 10.1002/hbm.25664

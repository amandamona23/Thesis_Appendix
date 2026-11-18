# Thesis_Appendix
ALL code used for Amanda Awad's Thesis Project
# White matter ageing in the fornix and uncinate fasciculus: analysis code

Code for an honours thesis on age-related white matter change in the fornix (FX) and uncinate fasciculus (UF), using DTI, fixel-based analysis (FBA), hippocampal-prefrontal structural connectivity and lifestyle measures. Data are from HCP-Aging and AABC Release 2 (80 participants, 160 sessions: first and last diffusion visit).

**No imaging or participant data are included.** The HCP-Aging/AABC Data Use Agreement does not allow data to be shared outside the institution. Scripts are provided for transparency and reproducibility only. Paths, usernames and session lists must be edited before use.

## How to run

1. Edit paths and parameters in `pipeline_config_v2.sh` (every step script sources it).
2. Run steps in the order below. Steps take session IDs (`HCAxxxxxxx_Vx_MR`) on the command line or with `--list sessions.txt`.
3. Run `make_dataset.py`, then the R scripts.

Software versions are listed in Table A2 of the thesis appendix (MRtrix3 3.0.4, FreeSurfer 7.1.1).

## Table A1. Script index

### Configuration

| Step | Script | What it does | Inputs | Outputs |
|---|---|---|---|---|
| Config | `pipeline_config_v2.sh` | Shared paths, thread limits and parameters (upsampled voxel size, DTI shells, streamline counts, tract list) | None | Variables sourced by all steps |

### Preprocessing

| Step | Script | What it does | Inputs | Outputs |
|---|---|---|---|---|
| 1 | `step1_per_subject_preproc.sh` | Denoising (`dwidenoise`), bias correction (`dwibiascorrect ants`), per-session response functions (`dwi2response dhollander`, native resolution), upsampling and mask regridding (`mrgrid`); FreeSurfer `recon-all` | HCP-preprocessed `T1w/Diffusion/data.nii.gz`, `bvals`, `bvecs`, `nodif_brain_mask.nii.gz`; `T1w_acpc_dc_restore` (and T2w) | `dwi_unbiased.mif`, `dwi_upsampled.mif`, `mask_upsampled.mif`, `response/{wm,gm,csf}.txt`; FreeSurfer subject folder |
| 2 | `step2_group_responsemean.sh` | Averages response functions across sessions | Per-session `response/*.txt` | `group_response/group_{wm,gm,csf}.txt` |
| 3a | `step3a_fod_normalise_v2.sh` | Multi-shell multi-tissue CSD (`dwi2fod msmt_csd`) with group responses; intensity normalisation (`mtnormalise`) | `dwi_upsampled.mif`, `mask_upsampled.mif`, group responses | `wmfod_norm.mif`, `gmfod_norm.mif`, `csffod_norm.mif` |
| 3c | `step3c_dti_tensor_v2.sh` | Tensor fit on the low-b shell at native resolution; FA, MD, RD, AD maps | `dwi_unbiased.mif`, mask | `dti/{fa,md,rd,ad}.mif` |
| 3c QC | `step3c_qc_check.sh` | Checks FOD and DTI outputs exist and are plausible | Step 3a and 3c outputs | Printed QC summary |

### Templates and fixel-based analysis

| Step | Script | What it does | Inputs | Outputs |
|---|---|---|---|---|
| 4a | `step4a_intrasubject_template_v3.sh` | Within-subject halfway FOD template from the two visits | `wmfod_norm.mif` (both visits) | `subject_templates/<SUBJECT>/halfway_template.mif`, `halfway_mask.mif`, `warps/<SESSION>.mif` |
| 4b | `step4b_group_template_v2.sh` | Population FOD template from 40 participants (10 per decade; list in `template40_subjects.txt`) | Halfway templates of the 40 participants | `group/wmfod_template.mif`, `group/template_mask.mif` |
| 4c | `step4c_halfway_to_group.sh` | Registers all 80 halfway templates to the population template; builds the analysis mask | Halfway templates, `wmfod_template.mif` | `halfway2group_warp.mif`, `group2halfway_warp.mif` |
| 5 | `step5_fixels_v3.sh` | Template fixel mask (`fod2fixel`, threshold 0.06); composes session warps; `fixelcorrespondence`; FD, FC, log FC and FDC per session | `wmfod_norm.mif`, step 4a and 4c warps, `wmfod_template.mif` | `fixel_mask/`, `session_warps/<SESSION>.mif`, `fd/`, `fc/`, `log_fc/`, `fdc/` (and `_fmls0p06` versions) |
| 6 | `step6_tracts_v3.sh` | TractSeg on the template (`sh2peaks`, `tract_segmentation`, `TOM`, `Tracking`, 5,000 streamlines per tract); tract voxel masks; tract fixel masks (`tck2fixel`, streamline cutoffs 1, 3, 10) | `wmfod_template.mif`, `template_mask.mif`, `fixel_mask/` | `tracts/tractseg/`, `tracts/voxel_masks/`, `tracts/fixels/<TRACT>/<TRACT>_ge<N>.mif`, `tract_summary.tsv` |
| 6 (tables) | `step6_tracts_v3.sh --zero-table`, `make_fdc_table.sh` | Per-session tract mean FD, log FC and FDC | Step 5 outputs, tract fixel masks | `tracts/zeros_by_session.tsv`, `tracts/fdc_by_session.tsv` |
| 7 | `[script name to add]` | Whole-brain template tractogram (20M, SIFT to 2M), fixel-fixel connectivity matrix, connectivity-based smoothing of FD, log FC, FDC | `wmfod_template.mif`, `fixel_mask/`, step 5 outputs | `step7/` (analysis mask, `tracks_sift.tck`, `matrix/`), `fd_smooth/`, `log_fc_smooth/`, `fdc_smooth/` |

### DTI tract values

| Step | Script | What it does | Inputs | Outputs |
|---|---|---|---|---|
| 8 | `step8_dti_tracts_v3.sh` | Warps DTI maps to template space with the composed warps; tract-mean FA, MD, RD, AD in the tract voxel masks | `dti/*.mif`, `session_warps/`, `tracts/fixels/` (via `fixel2voxel`) | `dti_tracts/masks/`, `dti_tracts/dti_by_session.tsv` |
| 8b | `step8b_along_tract_v3.sh` | Along-tract profiles (explored, not used in the thesis) | Step 6 and 8 outputs | `alongtract/` |

### Fixel-wise statistics

| Step | Script | What it does | Inputs | Outputs |
|---|---|---|---|---|
| 9 | `step9_make_change.sh` | Annual rate-of-change images per participant | Smoothed FD, log FC, FDC for both visits; years between visits | `change/{fd,log_fc,fdc}_rate/` |
| 9 | `step9_fixelcfestats.sh` | `fixelcfestats` (CFE, 5,000 permutations) within each tract mask; cross-sectional and rate-of-change models | Smoothed fixel data, `matrix/`, design and contrast files, tract masks | `stats/<model>_<metric>_<tract>_ge3/` |
| 9 | `step9_summary.sh` | Significant fixels per hypothesis (FWE p < 0.05), percentage of tract, effect size; optional streamline figure files | `stats/` results | `stats/step9_summary.tsv`, `stats/figures/*.tsf` |

### Connectome

| Step | Script | What it does | Inputs | Outputs |
|---|---|---|---|---|
| 10 | `step10_connectome.sh` | Per session, in native space: `5ttgen hsvs`, `5tt2gmwmi`, `labelconvert` (84 Desikan-Killiany nodes), `tckgen` (iFOD2, ACT, 10M streamlines), `tcksift2`, `tck2connectome`, hippocampus-prefrontal edge streamlines, QC | `wmfod_norm.mif`, FreeSurfer subject folder | `connectome/<SESSION>/` (`5tt.mif.gz`, `gmwmi.mif.gz`, `nodes.mif`, `connectome.csv`, `sift_mu.txt`, `edges/`, `tracks_200k.tck`, `qc.tsv`), `connectome_qc.tsv` |
| 10 | `step10_extract.py` | Fibre bundle capacity (SIFT2 weights x mu) for hippocampal-prefrontal edges, left and right | `connectome/<SESSION>/connectome.csv`, `sift_mu.txt` | `connectome_edges.tsv`, `connectome_reliability.tsv` |

### Covariates and analysis tables

| Step | Script | What it does | Inputs | Outputs |
|---|---|---|---|---|
| ICV | `overnight_icv.sh` | Intracranial and ventricle volumes with SAMSEG; Talairach/eTIV comparison | FreeSurfer subject folders | `icv_samseg.tsv`, `icv_comparison.tsv` |
| Covariates | `build_covariates.py` | Session and subject covariate tables and CFE design matrices (age, sex, site, education, ICV, ventricles) | AABC non-imaging data, ICV tables, session list | `covariates/covariates_by_session.tsv`, `covariates_by_subject.tsv`, `design/` |
| Dataset | `make_dataset.py` | Merges covariates, DTI, FBA and connectome tract values | Covariate tables, `dti_by_session.tsv`, `zeros_by_session.tsv`, `fdc_by_session.tsv`, `connectome_edges.tsv` | `thesis_by_session.csv` (160 rows), `thesis_by_subject.csv` (80 rows), `data_dictionary.csv` |

### Statistics and figures (R)

All R scripts read `thesis_by_session.csv` (or `thesis_by_subject.csv`) and write to `thesis_results/`.

| Step | Script | What it does | Outputs |
|---|---|---|---|
| Stats | `eda_assumptions.R` | Distributions, outliers, collinearity and model assumption checks | Diagnostic plots |
| Stats | `thesis_final.R` | Main linear mixed-effects models (`age_dec * time + sex + educ_c + icv_c + site` with a random intercept per subject), FDR correction, cohort table, main figures | `tables/*.csv`, `figures/*` |
| Stats | `thesis_robustness.R` | Influence, cluster-robust SEs, unequal variance, change-score, site x time and visit-offset checks | `tableS_robustness*.csv` |
| Stats | `thesis_ventricle_adjusted.R` | Models adjusted for lateral ventricle volume | Ventricle-adjusted tables |
| Stats | `thesis_connectome.R` | Hippocampal-prefrontal connectivity models, left/right and sensitivity analyses | Connectome tables and figures |
| Stats | `thesis_lifestyle.R` | Primary and secondary cognition and sleep models (RAVLT, Trails B, PSQI) | Lifestyle tables and figures |
| Stats | `thesis_lifestyle_exploratory.R` | Post hoc exploratory sleep analysis | Exploratory tables |
| Stats | `quick_stats.py` | First-pass tract-level mixed models (Python, statsmodels) | `results/quick_stats/` |
| Figures | `thesis_figures.R`, `thesis_scatter_figures.R`, `thesis_scatter_figures_main.R`, `abstract_figures.R` | Thesis and abstract figures | `figures/` |
| Tables | `make_word_tables.R`, `make_main_word_table.R`, `make_vent_word_table.R` | Word versions of the results tables | `tables/*.docx` |

### Utilities (not part of the analysis)

`backup_to_rds.sh`, `nas_to_rds.sh` (archiving), `view_cfe.sh`, `view_step9_summary.R` (result inspection), `run4a_lanes.sh` (parallel launcher for step 4a).

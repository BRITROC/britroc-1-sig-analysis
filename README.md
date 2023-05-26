# britroc-1-sig-analysis

Compositional analyses for the data presented in the publication [The copy number and mutational landscape of recurrent ovarian high-grade serous carcinoma](https://www.medrxiv.org/content/10.1101/2022.10.21.22280992v1) (Smith & Bradley _et al._ 2023 - Accepted _Nature Communications_).

## Compositional CN signatures
### setup
Clone this git repository and `cd` into the directory
```
git clone https://github.com/BRITROC/britroc-1-sig-analysis.git
cd britroc-1-sig-analysis
```
#### Install the conda environment


#### Activate the newly installed environment
```
conda activate britroc-fitting
```

### Run analysis
```
#./run_analysis_pipeline.sh
```

## Organisation of files
### Main files
```
.
├── README.md
├── code
│   ├── exploratory_analyses
│   │   ├── analysis_exposures.R
│   │   └── signatures_by_location.pdf
│   ├── models
│   │   ├── run_bernoulli_zeros
│   │   ├── run_partialILR
│   │   └── tmb_RE
│   ├── models_other_regressions
│   │   └── LN_modelling_zeros_FE_partialILR_TMB_TP53.R
│   └── signature_general
├── data
├── out
│   ├── inference
│   │   ├── correlated_binom_0_63ddb8a4-d766-4b28-bb2d-3d4b0c2ebf62.Data
│   │   ├── correlated_binom_0_63ddb8a4-d766-4b28-bb2d-3d4b0c2ebf62.RDS
│   │   └── partialILR
│   │       ├── partialILRcor_nlminb_allsigs.RDS
│   │       ├── partialILRnocor_nlminb_allsigs.RDS
│   │       ├── partialILRnocor_nlminb_allsigs_TMBdata.RDS
│   │       ├── repeat_full_with_different_initial_obs.RDS
│   │       ├── res_nlminb.RDS
│   │       ├── res_nlminb_nocoroutsidesd_only_matched.RDS
│   │       ├── res_nlminb_nocoroutsidesd_only_matched_allsigs.RDS
│   │       ├── res_nlminb_outsidesd_only_matched.RDS
│   │       └── res_nlminb_outsidesd_only_matched_allsigs.RDS
├── results
│   ├── exploratory
│   │   ├── bleeding_signature_s5
│   ├── partialILRmodelling_FE
│   ├── partialILRmodelling_FE_other_regression
│   ├── partialILRmodelling_ME
│   │   └── ternary_all.pdf
│   └── zeros_modelling
└── text
    └── morrill_britroc_latex
        ├── morrill_britroc_latex.tex
        ├── morrill_britroc_latex_content.tex
        └── morrill_britroc_latex_content2.tex
```

### Models
The models are implemented using Template Model Builder (TMB) and are found in `code/models/tmb_RE/`. They are run from R (see the files under `code/models/`).

### Author(s)
- [@lm687](https://github.com/lm687)
- [@phil9s](https://github.com/Phil9S)

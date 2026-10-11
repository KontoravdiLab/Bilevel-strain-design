# Bilevel-strain-design
Bilevel strain design using a genome-scale model of CHO cells
# Data and code deposit

Models, constraints, optimisation code, scenario tables and analysis code for the
bilevel strain-design study of IgG production in CHO-K1. `notebooks/01_reproduce.ipynb`
reproduces every computational result in the paper, with outputs stored. The
BCAT2-restored and iCHO3K transfer model variants are built there at run time
rather than deposited as files.

## Contents

    models/          published CHO-K1 model, the modified model used (iCHO_K1_used),
                     iCHO1766 and iCHO3K
    constraints/     flux bounds for all 3952 reactions, one file per dataset
    preprocessing/   ccFBA, ncFBA and FVA code (MATLAB) that produced the bounds
    optimisation/    bilevel.gms and its inputs; regulation-removal scripts; RKO_Results.m
    scenarios.xlsx   93 initial, 89 refined and 20 high-producing scenarios
    notebooks/       01_reproduce.ipynb, eight steps, outputs stored
    data/            fed-batch measurements
    figures/         code that regenerates Figures 1, 4 and S3

In `scenarios.xlsx`, `scenario_id` is the numbering used in the paper; `id_in_93`
is the numbering of the optimisation logs. The two agree up to 67 and diverge
thereafter, because four duplicate scenarios (93-ids 68, 74, 87 and 89) were
removed at the refinement step.

## Datasets

| File | Dataset |
|---|---|
| `FVA_EX_SUE_K1_NoEPO.mat` | D1 (CHO-T, FeedC) |
| `FVA_EX_SAR_FeedC_K1_NoEPO.mat` | D2 (CHO-GS46, FeedC) |
| `FVA_EX_SAR_FeedAll_K1_NoEPO.mat` | D3 (CHO-GS46, FeedAll) |
| `FVA_EX_SAR_FeedAll40_K1_NoEPO.mat` | D4 (CHO-GS46, FeedAll + 40%) |

## Running

- `notebooks/01_reproduce.ipynb` and `figures/Make_Figures.ipynb`: Python, run in place.
- `optimisation/bilevel.gms`: GAMS with a CPLEX licence, run from `optimisation/`.
  It is a representative script; set the run settings at the top of the file.
  Read the header before running: it records two respects in which the deposited
  script differs from the scripts that produced the published scenarios (the
  corrected reversibility constraint, and the D2-specific run settings).
- `preprocessing/`: MATLAB; needs the source concentration data (see below).
  `PrintToGams.m` takes the directory of FVA bound files as its fourth argument
  and defaults to `../constraints`.

`optimisation/inputs/` holds one folder per dataset. `SarantosFeedAll40` (D4) has
no `ModelDrainRXN.txt`; the file is not read by `bilevel.gms` and its absence does
not affect any run.

For the Python code, use the `glpk_exact` solver, not default GLPK: with the design
applied, default GLPK returns solutions that violate the IgG demand bound.

## Not included

The concentration data behind D1-D4 come from the publications cited in the paper.
Western blot densitometry values are not available.

## Software

MATLAB R2019b, COBRA Toolbox v3.0, GAMS 28.2.0, CPLEX 12.9.0.0.
Python 3.12, cobrapy 0.32.1, optlang 1.9.1, scipy, openpyxl, matplotlib.

## Licence and sources

The CHO-K1 model and iCHO1766: Hefzi et al. (2016), Cell Systems 3, 434-443,
doi:10.1016/j.cels.2016.10.020; reuse subject to the publisher's terms.

iCHO3K: Di Giusto, Choi, Antonakoudis et al. (2026), Cell Systems 17, 101574,
doi:10.1016/j.cels.2026.101574, https://github.com/LewisLabUCSD/iCHO3K, CC BY 4.0.
`iCHO3K.mat` is unmodified; the transfer analysis closes LEUTA, LEUTAm,
EX_pheleu_e, HMR_5166 and r1402 at run time in the notebook.

`carbonConstraints.m` and `findCoAs.m` are edited from M. Lularevic's ccFBA code.
Everything else is original to this work.

## Figures

`figures/Make_Figures.ipynb` generates Figures 1, 4 and S3 from the fed-batch data.
Figures 2 (an Escher flux map) and 3 (Western blots) are not generated from code;
the uncropped membranes behind Figure 3 are in Supplementary Figure S2.

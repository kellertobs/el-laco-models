# el-laco-models
Materials to reproduce thermodynamic and mechanical modelling in the manuscript "Extrusion of iron-rich melt formed magnetite-apatite deposits on El Laco volcano" by Keller et al., under review.

## Thermodynamic model: `/alphaMELTS`
Results were produced with alphaMELTS v.1.9 (available at https://magmasource.caltech.edu/alphamelts/). All models including liquid immiscibility were run as isobaric crystallisation sequences with the general parameters listed in "isobaric_immisc_xtaln.txt". The alternative model with only one liquid phase allowed was run with the parameters in "isobaric_oneliq_xtaln.txt".

All results for El Laco andesite compositions are found in the folder `alphaMELTS/LACO`, using the same abbreviations to identify model compositions as in the main text and figures, as well as extended data figures and tables. Reproduced experiments from Hou et al., Nature Comms, 2018, are in folder `/alphaMELTS/HOU+2018_REPRODUCED` with individual model runs labelled using the same abbreviations as in the original paper.


## Volcano deformation model: `/VolcanoDeformation`
The source code for the finite-element volcano deformation is located in `/VolcanoDeformation/src`. The code was built and tested with Matlab R2020a. User scripts to reproduce results in the main text and figures, as well as extended data figures, are provided in `/VolcanoDeformation/usr`. To run the code, open and execute a user script (e.g., `run_LACO.m`) in Matlab. The output will be saved to `/VolcanoDeformation/out`

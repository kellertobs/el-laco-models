# el-laco-models
Materials to reproduce thermodynamic and mechanical models and scaling analyses in the manuscript "A genetic model of the magnetite-apatite deposits on El1
Laco volcano by extrusion of iron-rich melt" by Keller et al., in review __Nature Comms__.

## Thermodynamic model: `/alphaMELTS`
Results were produced with alphaMELTS v.1.9 (available at https://magmasource.caltech.edu/alphamelts/). All models including liquid immiscibility were run as isobaric crystallisation sequences with the general parameters listed in "isobaric_immisc_xtaln.txt". The alternative model with only one liquid phase allowed was run with the parameters in "isobaric_oneliq_xtaln.txt".

All results for El Laco andesite compositions are found in the folder `alphaMELTS/LACO`, using the same abbreviations to identify model compositions as in the main text and figures, as well as extended data figures and tables. Reproduced experiments from Hou et al., Nature Comms, 2018, are in folder `/alphaMELTS/HOU+2018_REPRODUCED` with individual model runs labelled using the same abbreviations as in the original paper.

The file `alphaMELTS/fit_remix.m` contains a Matlab routine used to reconstruct an averaged El Laco host andesite composition by mixing analysed compositional end-members including phenocryst phases in host andesite and immiscible melts analysed in unmixed melt inclusions. To reproduce the fitting by least-squares run the script in Matlab.

## Fe-rich melt separation scaling analysis: `/PhaseSeparation`
Provided are Matlab scripts that reproduce the immiscible melt viscosity model fitted to viscometry measurements (`/PhaseSeparation/fit_viscosity_model.m`), as well as the three-phase transport coefficient calibration (`/PhaseSeparation/ThreePhaseCalibration.m`) and scaling analysis (`/PhaseSeparation/ThreePhaseScalingAnalysis.m`) used to constrain rates and timing of Fe-rich melt segregating from a variably crystalline parent magma/mush. The latter two depend on the `ternplot` library, an open-source plotting tool for ternary plots (https://mathworks.com/matlabcentral/fileexchange/2299-alchemyst-ternplot), which is included here (`/PhaseSeparation/ternplot`). To reproduce the models, run the scripts in Matlab. 

## Volcano deformation model: `/VolcanoDeformation`
The source code for the finite-element volcano deformation is located in `/VolcanoDeformation/src`. The code was built and tested with Matlab R2020a. User scripts to reproduce results in the main text and figures, as well as extended data figures, are provided in `/VolcanoDeformation/usr`. To run the code, open and execute a user scripts (e.g., `run_LACO.m`) in Matlab. The output will be saved to `/VolcanoDeformation/out`

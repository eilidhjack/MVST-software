# MVST-software
Software and data to reproduce analysis from "Estimating the changing nature of Scotland's health inequalities using a multivariate spatio-temporal model".

Description of files:

"Ten Years Three Disease - IG.csv" - Contains the disease and covariate data. 

"Model.run.R" - R file which will load required libraries, read in the data, source the necessary files and run the model.

"MVLerouxAR.R" - Main R function.

"Leroux.phi.cpp" - cpp file which is called from main function to sample the random effects.

"AR1.HB.cpp" - cpp file which is called from main function to sample the HB effects.

"kronecker.cpp" - cpp file which calculates the kronecker product which is needed when sampling rho.

"W.RData" - Neighbourhood matrix, W, which has been altered to connect the islands to the mainland.


# SIR animated maps
Animations of the SIR maps for each disease across time are available here to download.

Description of files:

"SuppCereSIR.pdf" - Standardised Incidence Ratios (SIR) for cerebrovascular disease for each IG in Scotland from
2003 to 2012.

"SuppCHDSIR.pdf" - Standardised Incidence Ratios (SIR) for coronary heart disease for each IG in Scotland from
2003 to 2012.

"SuppRespSIR.pdf" - Standardised Incidence Ratios (SIR) for respiratory disease for each IG in Scotland from
2003 to 2012.

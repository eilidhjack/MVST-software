# MVST-software
Software and data to reproduce analysis from "Estimating the changing nature of Scotland's health inequalities using a multivariate spatio-temporal model"

Description of files:
"Ten Years Three Disease - IG.csv" - Contains the disease and covariate data. 
"Model.run.R" - R file which will load required libraries, read in the data, source the necessary files and run the model.
"MVLerouxAR.R" - Main R function.
"Leroux.phi.cpp" - cpp file which is called from main function to sample the random effects.
"AR1.HB.cpp" - cpp file which is called from main function to sample the HB effects.
"kronecker.cpp" - cpp file which calculates the kronecker product which is needed when sampling rho.

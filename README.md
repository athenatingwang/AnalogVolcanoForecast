# AnalogVolcanoForecast

The computational code and data are for the paper: Wang, T., Bebbington, M., Cronin, S., & Carman, J. (2022). Forecasting eruptions at poorly known volcanoes using analogs and multivariate renewal processes. Geophysical Research Letters, 49, e2021GL096715. https://doi.org/10.1029/2021GL096715.

All .R and .jags files can be downloaded and saved in the same folder. Instructions are given if they depend on code in other sections, in which case one simply needs to run the R code in that section first.

### GVPmod0WeibsepZuse.jags 
This file contains the JAGS "model" for the Weibull model in equation (1) in the paper Wang et al. (2022).

### GVPmod1WeibsepZuse.jags
This file contains the JAGS "model" for model M1 in equation (2) in the paper Wang et al. (2022).

### GVPmod2WeibsepZuse.jags
This file contains the JAGS "model" for model M2 in equation (3) in the paper Wang et al. (2022).

### GVPmod3WeibsepZuse.jags
This file contains the JAGS "model" for model M3 in equation (4) in the paper Wang et al. (2022).

### ReadGVPdata.R
This file reads all the data files for the empirical analogues saved in the folder "HoloceneRecords". These data files are obtained from the Smithsonian Institution's GVP catalogue.

##
## When fitting models, make sure to save all .R and .jags files and the two data folders "HoloceneRecords" and "QuaternaryRecords" in the same main folder.

## Fitting models M1, M2, and M3 to the empirical analogues 
### GVPEmpAnalog.R
Run the file GVPEmpAnalog.R, and the MCMC samples of the posterior distributions will be saved in .image files.

## Fitting models M1, M2, and M3 to the statistical analogues with Holocene records, and carry out residual analysis 
### GVPStatsAnalog.R
Run the file GVPStatsAnalog.R, and the MCMC samples of the posterior distributions will be saved in .image files. The residual analysis for each model will be saved as .eps file. If a different format is preferred, such as a .pdf file, just change the "postscript" command to "pdf" and change ".eps" to ".pdf" in the R code.

## Fitting models M1, M2, and M3 to the statistical analogues with Quaternary records (VEI 4+), and carry out residual analysis 
### QuatStatsAnalogVEI4+.R
Run the file QuatStatsAnalogVEI4+.R, and the MCMC samples of the posterior distributions will be saved in .image files. The residual analysis for each model will be saved as .eps file. If a different format is preferred, such as a .pdf file, just change the "postscript" command to "pdf" and change ".eps" to ".pdf" in the R code.

## Fitting models M1, M2, and M3 to the statistical analogues with Quaternary records (VEI 3+), and carry out residual analysis 
### QuatStatsAnalogVEI3+.R
Run the file QuatStatsAnalogVEI3+.R, and the MCMC samples of the posterior distributions will be saved in .image files. The residual analysis for each model will be saved as .eps file. If a different format is preferred, such as a .pdf file, just change the "postscript" command to "pdf" and change ".eps" to ".pdf" in the R code.

## Fitting model M0, a simple Weibull renewal process to the Tongariro GVP Holocene record only 
### GVP-TgrOnly-VEI3+.R
Run the file GVP-TgrOnly-VEI3+.R, and the MCMC samples of the posterior distributions will be saved in .image files. It also provides forecasts from using Holocene record only from Tongariro volcano.

## Fitting model M0, a simple Weibull renewal process to the Tongariro Quaternary record only 
### Quat-TgrOnly-VEI3+.R 
Run the file Quat-TgrOnly-VEI3+.R, and the MCMC samples of the posterior distributions will be saved in .image files. It also provides forecasts from using Quaternary record only from Tongariro volcano.


## Plotting the VEI, posterior distributions, and forecasts used in the paper 
### GRLplots.R
All the results files used in this code are saved in the Results folder.





